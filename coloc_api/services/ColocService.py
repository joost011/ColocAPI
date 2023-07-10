import subprocess, json, itertools, re, os, time
import pandas as pd
import numpy as np
import multiprocessing as mp
from django.conf import settings
from django_rq import job
from ..models import ColocAnalysis, ColocAnalysisStatus

class ColocService:
    '''Class that is responsible for performing colocalization analysis on an input file'''

    def __init__(self):
        self.manager = mp.Manager()
        self.status_message = self.manager.Value('c', '')
        self.status_order = self.manager.Value('i', 1)
        self.status_update = True

    def start_coloc(self, file_name): 
        '''Main function that initializes the colocalization process.
        It has the @job decorator from Redis so it will be run in a different process on the Redis server'''
        uuid = file_name.split('.')[0]
        status_update_interval = 0.25 # Seconds
        status_update_process = mp.Process(target=self.schedule_status_update, args=[uuid, status_update_interval])
        status_update_process.start()

        # Start analysis
        self.process_file(file_name)

        with open(settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json', 'r') as f:
            output = json.load(f)

        return output
    

    def schedule_status_update(self, uuid, interval):
        '''Function that starts the loop in which the process status is updated'''
        while self.status_update:
            self.update_status(uuid)
            time.sleep(interval)


    def update_status(self, uuid):
        '''Function that updates the status of a coloc process in the database'''
        if self.status_update != '':
            coloc_analysis_object = ColocAnalysis.objects.get(uuid=uuid)
            coloc_analysis_status_object, created = ColocAnalysisStatus.objects.update_or_create(coloc_analysis=coloc_analysis_object, status_order=self.status_order.value, defaults={
                'coloc_analysis': coloc_analysis_object, 
                'status_message': self.status_message.value, 
                'status_order': self.status_order.value
                })

            coloc_analysis_status_object.status_message = self.status_message.value
            coloc_analysis_status_object.status_order = self.status_order.value
            coloc_analysis_status_object.save()

            if self.status_message.value.lower() == 'done':
                coloc_analysis_object.finished = True
                coloc_analysis_object.save()


    def process_file(self, file_name):
        '''Function that starts processing the input file'''
        in_file_path = settings.STORAGE['IN_FILES_PATH'] / file_name
        uuid = file_name.split('.')[0]
        window_dict = self.get_windows(in_file_path)
        genes = self.get_genes(window_dict)
        gene_ids, gene_positions_dict = self.get_gene_info(genes)

        # Create processed and output directory with uuid
        os.mkdir(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid)
        os.mkdir(settings.STORAGE['OUT_FILES_PATH'] / uuid)

        # Create int value that can be shared in separate processes
        manager = mp.Manager()
        num_coloc = manager.Value('i', 0)

        gwas_snps_dict = self.obtain_gwas_snps(in_file_path, gene_ids, gene_positions_dict)

        # Analyse every gene in a separate process
        with mp.Pool(mp.cpu_count() - 2) as pool:
            args_list = [(uuid, gene_ids[i], gene_positions_dict[gene_ids[i]], gwas_snps_dict[gene_ids[i]], num_coloc) for i in range(len(gene_ids))]
            pool.starmap(self.process_gene, args_list)

        self.status_message.value = 'Preparing output'
        self.status_order.value = 5
        self.create_output_file(uuid)


    def process_gene(self, uuid, gene_id, gene_position, gwas_snps, num_coloc):
        '''Function that starts processing a gene'''
        if len(gwas_snps) > 0:
            eqtl_snps = self.obtain_eqtl_snps(gene_id)

            if len(eqtl_snps) > 0:
                if self.prepare_coloc_data(gene_id, eqtl_snps, gwas_snps, gene_position, uuid):
                    self.perform_coloc_all_genes(gene_id, uuid)

        num_coloc.value += 1
        self.status_message.value = f'{num_coloc.value} genes processed'
        self.status_order.value = 4


    def merge_windows(self, arr):
            '''Function that clusters an array of windows'''
            arr.sort(key=lambda x: x[0])
            index = 0
        
            for i in range(1, len(arr)):

                if (arr[index][1] >= arr[i][0]):
                    arr[index][1] = max(arr[index][1], arr[i][1])
                else:
                    index = index + 1
                    arr[index] = arr[i]
        
            res = []
            for i in range(index+1):
                res.append(arr[i])

            return res


    def get_windows(self, file_path):    
        '''Function that gets windows around significant snps from GWAS file'''
        window_dict = {}
        num_significant_snps = 0

        # Loop through GWAS in chunks of 1.000.000
        for df in pd.read_csv(file_path, sep='\t', comment='#', chunksize=1000000):
            
            # Identify significant snps (cutoff of 5 × 10−8)
            df = df.loc[df['p_value'] < 5e-8]

            # Update analysis total number of significant snps found
            num_significant_snps += len(df)

            # For each significant snp, create start and end window columns
            df['window_start'] = df['base_pair_location'] - 10000
            df['window_end'] = df['base_pair_location'] + 10000

            # Group the data by chromosome number
            grouped = df.groupby('chromosome')

            # Loop throug all the chromosomes in the grouped dataset
            for k, v in grouped:

                start_list = v['window_start'].tolist()
                end_list = v['window_end'].tolist()
                windows = [[start_list[i], end_list[i]] for i in range(len(start_list))]

                # Merge/cluster the windows of significant snps so genes don't get colocalized twice
                merged_windows = self.merge_windows(windows)

                # Add the merged windows to a dictionary where the windows of all chunks are saved
                try:
                    window_dict[k] += merged_windows
                except:
                    window_dict[k] = merged_windows


            self.status_message.value = f'{num_significant_snps} significant snps found on {len(window_dict)} chromosomes'
            self.status_order.value = 1

        return window_dict
        
        
    def get_genes(self, window_dict):
        '''Function that retrieves genes from the Gencode file based on the windows around significant snps'''

        # Open de gencode file (file with position data of snps)
        gene_store = pd.HDFStore(settings.STORAGE['STATIC_FILES']['GENE_POSITIONS'])
        genes = []
        num_genes = 0

        # Retrieve all genes that are within the windows of significant snps
        for chromosome, windows in window_dict.items():
            for window in windows:
                c = f'chr{chromosome}'
                data = gene_store.select('genes', where='chromosome == c and ((start >= window[0] and start <= window[1]) or (end >= window[0] and end <= window[1]))')
                num_genes += len(data)

                self.status_message.value = f'{num_genes} genes found in windows of significant snps'
                self.status_order.value = 2
                genes.append(data)
                
        gene_store.close()
            
        return genes


    def get_gene_info(self, genes):
        '''Function that takes a list of dataframes with gene info as input and converts it to a list'''

        # Create a list with all gene ids
        gene_ids = [df['gene_id'].to_list() for df in genes]
        gene_ids = list(itertools.chain(*gene_ids))
        gene_ids = [gene_id.split('.')[0] for gene_id in gene_ids]

        # Create a list with all gene names
        gene_names = [df['gene_name'].to_list() for df in genes]
        gene_names = list(itertools.chain(*gene_names))
        gene_names = [gene_name.split('.')[0] for gene_name in gene_names]

        # Create a list with all start positions
        gene_start_positions = [df['start'].to_list() for df in genes]
        gene_start_positions = list(itertools.chain(*gene_start_positions))

        # Create a list with all end positions
        gene_end_positions = [df['end'].to_list() for df in genes]
        gene_end_positions = list(itertools.chain(*gene_end_positions))

        # Create a list with all chromosome numbers of the genes
        # This should be altered, because the sex chromosomes are missing now
        gene_chromosomes = [df['chromosome'].to_list() for df in genes]
        gene_chromosomes = list(itertools.chain(*gene_chromosomes))
        gene_chromosomes = [re.sub('\D', '', chromosome) if re.sub('\D', '', chromosome) != '' else chromosome[-1] for chromosome in gene_chromosomes]

        # Combine the above created lists
        gene_positions = [(gene_start_positions[i], gene_end_positions[i], gene_chromosomes[i], gene_names[i]) for i in range(len(gene_start_positions))]

        # Put the items of the combined list to a dict with the corresponding gene as key
        gene_positions_dict = {gene_ids[i]: gene_positions[i] for i in range(len(gene_ids))}
        
        return gene_ids, gene_positions_dict


    def obtain_gwas_snps(self, file_path, gene_ids, gene_positions_dict):
        '''Function that extracts SNPs from GWAS file for every locus'''
        gwas_snps = {}

        # Loop throug GWAS again, now that we know the location of the genes
        for df in pd.read_csv(file_path, sep='\t', comment='#', chunksize=1000000):

            # Process data
            df = df[['variant_id', 'chromosome', 'base_pair_location', 'p_value', 'beta', 'standard_error']]
            df['varbeta'] = df['standard_error']**2
            df.rename(columns={'variant_id': 'snp'}, inplace=True)
            df.loc[df['p_value'] == 0, 'p_value'] = 2.22e-308
            df['logp'] = -np.log10(df['p_value'])
            
            # Loop through all genes
            for gene in gene_ids:

                gene_position = gene_positions_dict[gene]
                    
                try:
                    chromosome = int(gene_position[2])
                except:
                    chromosome = gene_position[2]

                # Select all snps from GWAS data for every gene
                data = df.loc[
                        (df['base_pair_location'] >= int(gene_position[0]) - 100_000) & 
                        (df['base_pair_location'] <= int(gene_position[1]) + 100_000) &
                        (df['chromosome'] == chromosome)
                    ]
                
                # Add the data to a dict, with gene ids as key
                if gene in gwas_snps:
                    gwas_snps[gene] = pd.concat([gwas_snps[gene], data])
                else:
                    gwas_snps[gene] = data

                # Remove duplicate snps
                gwas_snps[gene].drop_duplicates(subset='snp', inplace=True)

                total_snp_extracted = sum([len(v) for k, v in gwas_snps.items()])
                self.status_message.value = f'{total_snp_extracted} SNPs extracted from GWAS data'
                self.status_order.value = 3

        return gwas_snps


    def obtain_eqtl_snps(self, gene_id):
        '''Function that retrieves all SNPs for a gene from the eQTLGen dataset'''
        data = pd.read_hdf(settings.STORAGE['STATIC_FILES']['EQTLGEN'], key='eqtls', where='Gene == gene_id')    
        return data


    def prepare_coloc_data(self, gene_id, eqtl_snps, gwas_snps, gene_positions, uuid):
        '''Function that prepares the datasets for coloc analysis'''

        # Store total number of GWAS and eQTL snps in variables
        total_gwas_snps = len(gwas_snps)
        total_eqtl_snps = len(eqtl_snps)

        # Get intersecting snps
        eqtl_snps_list = eqtl_snps['snp'].tolist()
        gwas_snps = gwas_snps.loc[gwas_snps['snp'].isin(eqtl_snps_list)]
        gwas_snps_list = gwas_snps['snp'].tolist()
        eqtl_snps = eqtl_snps.loc[eqtl_snps['snp'].isin(gwas_snps_list)]

        # Store total number of intersecting snps in variable
        total_intersecting_snps = len(gwas_snps)

        # If the length of the data in the GWAS dict and the eQTLGEN dict is not 0, write them to a JSON file, which can later be loaded in the Coloc R script
        if len(gwas_snps) > 0 and len(gwas_snps) > 0:

            meta_data = {
                'total_gwas_snps': total_gwas_snps,
                'total_eqtl_snps': total_eqtl_snps,
                'total_intersecting_snps': total_intersecting_snps,
                'chromosome': gene_positions[2],
                'start_position': gene_positions[0],
                'stop_position': gene_positions[1],
                'gene_name': gene_positions[3],
            }

            d = {
                'meta_data': meta_data,
                'gwas': gwas_snps.to_dict(orient='list'), 
                'eqtls': eqtl_snps.to_dict(orient='list')
                    }

            gene_file_name = gene_id + '.json'
            with open(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid / gene_file_name, 'w') as out_file:
                out_file.write(json.dumps(d))

            return True
        
        return False


    def perform_coloc_all_genes(gene_id, uuid):
        '''Function that calls the coloc R script in a subprocess'''
        subprocess.call(['Rscript', settings.R['COLOC_PATH'], gene_id, uuid, '--no-save'])


    def create_output_file(self, uuid):
        '''Function that creates an output file of the analysis'''
        directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        res_dict = {
            'genes': {},
        }

        for filename in os.listdir(directory):
            f = os.path.join(directory, filename)
            gene = f.split('/').pop().split('.')[0]

            # For every colocalizing gene, add gene object to output
            with open(f, 'r') as input:
                data = json.load(input)
                res_dict['genes'][gene] = data['datasets']

                # Add top snp to output
                pph4_list = res_dict['genes'][gene]['coloc']['SNP.PP.H4']
                max_pph4 = max(pph4_list)
                max_index = pph4_list.index(max_pph4)
                top_snp = res_dict['genes'][gene]['coloc']['snp'][max_index]
                res_dict['genes'][gene]['meta_data']['top_snp'] = top_snp

        # Add additional data to output
        res_dict['total_tested_genes'] = len(os.listdir(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid))
        res_dict['total_colocalizing_genes'] = len(os.listdir(directory))
        res_dict['total_gwas_snps'] = sum([v['meta_data']['total_gwas_snps'] for k, v in res_dict['genes'].items()])
        res_dict['total_eqtl_snps'] = sum([v['meta_data']['total_eqtl_snps'] for k, v in res_dict['genes'].items()])

        # Write output to file
        with open(settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json', 'w') as output:
            self.status_message.value = 'Done'
            self.status_order.value = 6
            self.status_update = False
            json.dump(res_dict, output)
        