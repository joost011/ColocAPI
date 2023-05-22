import subprocess
import pandas as pd
import json
import itertools
import numpy as np
import re
import os
import sqlite3
from django.conf import settings


def coloc(file_name):    
    process_file(file_name)
    uuid = file_name.split('.')[0]

    with open(settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json', 'r') as f:
        output = json.load(f)

    return output


def process_file(file_name):
    file_path = settings.STORAGE['IN_FILES_PATH'] / file_name
    uuid = file_name.split('.')[0]
    window_dict = get_windows(file_path, uuid)
    genes = get_genes(window_dict, uuid)
    gene_ids, gene_positions_dict = get_gene_info(genes)

    os.mkdir(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid)
    os.mkdir(settings.STORAGE['OUT_FILES_PATH'] / uuid)

    update_status(uuid, f'A total of {len(gene_ids)} genes is found in the windows of the significant SNPs')

    for i in range(len(gene_ids)):
        gene_id = gene_ids[i]
        gene_position = gene_positions_dict[gene_ids[i]]

        process_gene(file_path, uuid, gene_id, gene_position)
        update_status(uuid, f'Colocalization analysis performed for {i + 1}/{len(gene_ids)} genes')            
    
    update_status(uuid, f'Preparing output file')
    create_output_file(uuid)


def process_gene(file_path, uuid, gene_id, gene_position):
    gwas_snps = obtain_gwas_snps(file_path, gene_id, gene_position, uuid)

    if len(gwas_snps) > 0:
        eqtl_snps = obtain_eqtl_snps(gene_id, uuid)

        if len(eqtl_snps) > 0:
            # prepare_coloc_data(gene_id, eqtl_snps, gwas_snps, gene_position, uuid)

            if prepare_coloc_data(gene_id, eqtl_snps, gwas_snps, gene_position, uuid):
                perform_coloc_all_genes(gene_id, uuid)


def merge_windows(arr):
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


def update_status(uuid, message, finished=False):
    conn = sqlite3.connect(settings.DATABASES['default']['NAME'])
    cursor = conn.cursor()
    cursor.execute('''UPDATE coloc_api_colocanalysis SET status_message=?, finished=? WHERE uuid=?;''',(message, finished, uuid,))
    conn.commit()
    conn.close()


def get_windows(file_path, uuid):
    window_dict = {}
    num_significant_snips = 0

    # Loop through GWAS in chunks of 1.000.000
    for df in pd.read_csv(file_path, sep='\t', comment='#', chunksize=1000000):
        
        # Identify significant snips (cutoff of 5 × 10−8)
        df = df.loc[df['p_value'] < 5e-8]

        # Update analysis total number of significant snips found
        num_significant_snips += len(df)

        # For each significant snip, create start and end window columns
        df['window_start'] = df['base_pair_location'] - 10000
        df['window_end'] = df['base_pair_location'] + 10000

        # Group the data by chromosome number
        grouped = df.groupby('chromosome')

        # Loop throug all the chromosomes in the grouped dataset
        for k, v in grouped:

            start_list = v['window_start'].tolist()
            end_list = v['window_end'].tolist()
            windows = [[start_list[i], end_list[i]] for i in range(len(start_list))]

            # Merge/cluster the windows of significant snips so genes don't get colocalized twice
            merged_windows = merge_windows(windows)

            # Add the merged windows to a dictionary where the windows of all chunks are saved
            try:
                window_dict[k] += merged_windows
            except:
                window_dict[k] = merged_windows

        update_status(uuid, f'{num_significant_snips} significant snips found on {len(window_dict)} chromosomes')

    return window_dict
    
    
def get_genes(window_dict, uuid):
     # Open de gencode file (file with position data of snips)
    gene_store = pd.HDFStore(settings.STORAGE['STATIC_FILES']['GENE_POSITIONS'])
    genes = []
    num_genes = 0

    # Retrieve all genes that are within the windows
    for chromosome, windows in window_dict.items():
        for window in windows:
            c = f'chr{chromosome}'
            data = gene_store.select('genes', where='chromosome == c and ((start >= window[0] and start <= window[1]) or (end >= window[0] and end <= window[1]))')
            num_genes += len(data)
            update_status(uuid, f'{num_genes} genes found in windows of significant snips')
            genes.append(data)
            
    gene_store.close()
        
    return genes


def get_gene_info(genes):
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
    gene_chromosomes = [re.sub('\D', '', chromosome) for chromosome in gene_chromosomes]

    # Combine the above created lists
    gene_positions = [(gene_start_positions[i], gene_end_positions[i], gene_chromosomes[i], gene_names[i]) for i in range(len(gene_start_positions))]

    # Put the items of the combined list to a dict with the corresponding gene as key
    gene_positions_dict = {gene_ids[i]: gene_positions[i] for i in range(len(gene_ids))}
    
    return gene_ids, gene_positions_dict


def obtain_gwas_snps(file_path, gene_id, gene_positions_dict, uuid):
    snps = None

    # Loop throug GWAS again, now that we know the location of the genes
    for df in pd.read_csv(file_path, sep='\t', comment='#', chunksize=1000000):

        df = df[['variant_id', 'chromosome', 'base_pair_location', 'p_value', 'beta', 'standard_error']]

        # Do some data processing
        df['varbeta'] = df['standard_error']**2
        # df.loc[df['variant_id'].isna(), "variant_id"] = pd.util.testing.rands_array(10, sum(df['variant_id'].isnull()))
        df.rename(columns={'variant_id': 'snp'}, inplace=True)
        df.loc[df['p_value'] == 0, 'p_value'] = 2.22e-308
        df['logp'] = -np.log10(df['p_value'])
        
        # Loop through all genes
        # for gene in gene_ids:
            
        if gene_positions_dict[2] != '':
            chromosome = int(gene_positions_dict[2])
        else:
            chromosome = gene_positions_dict[2]

        # Select all snips from GWAS data for every gene
        data = df.loc[
                (df['base_pair_location'] >= int(gene_positions_dict[0]) - 100_000) & 
                (df['base_pair_location'] <= int(gene_positions_dict[1]) + 100_000) &
                (df['chromosome'] == chromosome)
            ]
        
        # Add the data to a dict, with gene ids as key
        if len(data) > 0:
            try:
                snps.append(data)
            except:
                snps = data

            # Remove duplicate snips (keep the first for now)
            data.drop_duplicates(subset='snp', inplace=True)


    return snps

def obtain_eqtl_snps(gene_id, uuid):
    print(gene_id)
    data = pd.read_hdf(settings.STORAGE['STATIC_FILES']['EQTLGEN'], key='eqtls', where='Gene == gene_id')    
    return data


def prepare_coloc_data(gene_id, eqtl_snips, gwas_snips, gene_positions, uuid):
    # Store total number of GWAS and eQTL snps in variables
    total_gwas_snips = len(gwas_snips)
    total_eqtl_snips = len(eqtl_snips)

    # Get intersecting snps
    eqtl_snips_list = eqtl_snips['snp'].tolist()
    gwas_snips = gwas_snips.loc[gwas_snips['snp'].isin(eqtl_snips_list)]
    gwas_snips_list = gwas_snips['snp'].tolist()
    gwas_snips = gwas_snips.loc[gwas_snips['snp'].isin(gwas_snips_list)]

    # Store total number of intersecting snips in variable
    total_intersecting_snips = len(gwas_snips)

    # If the length of the data in the GWAS dict and the eQTLGEN dict is not 0, write them to a JSON file, which can later be loaded in the Coloc R script
    if len(gwas_snips) > 0 and len(gwas_snips) > 0:

        meta_data = {
            'total_gwas_snps': total_gwas_snips,
            'total_eqtl_snps': total_eqtl_snips,
            'total_intersecting_snps': total_intersecting_snips,
            'chromosome': gene_positions[2],
            'start_position': gene_positions[0],
            'stop_position': gene_positions[1],
            'gene_name': gene_positions[3],
        }

        d = {
            'meta_data': meta_data,
            'gwas': gwas_snips.to_dict(orient='list'), 
            'eqtls': eqtl_snips.to_dict(orient='list')
                }

        gene_file_name = gene_id + '.json'
        with open(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid / gene_file_name, 'w') as out_file:
            out_file.write(json.dumps(d))

        return True
    
    return False


def perform_coloc_all_genes(gene_id, uuid):
    subprocess.call([settings.R['R_PATH'], settings.R['COLOC_PATH'], gene_id, uuid])


def create_output_file(uuid):
    directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
    res_dict = {
        'genes': {},
    }

    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        gene = f.split('\\').pop().split('.')[0]

        with open(f, 'r') as input:
            data = json.load(input)
            res_dict['genes'][gene] = data['datasets']

    res_dict['total_tested_genes'] = len(os.listdir(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid))
    res_dict['total_colocalizing_genes'] = len(os.listdir(directory))
    res_dict['total_gwas_snps'] = sum([v['meta_data']['total_gwas_snps'] for k, v in res_dict['genes'].items()])
    res_dict['total_eqtl_snps'] = sum([v['meta_data']['total_eqtl_snps'] for k, v in res_dict['genes'].items()])

    with open(settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json', 'w') as output:
        update_status(uuid, f'Done!', True)
        json.dump(res_dict, output)



    