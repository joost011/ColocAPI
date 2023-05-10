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


def process_file(file_name):
    file_path = settings.STORAGE['IN_FILES_PATH'] / file_name
    uuid = file_name.split('.')[0]

    window_dict = {}
    num_significant_snips = 0

    # Loop through GWAS in chunks of 1.000.000
    for df in pd.read_csv(file_path, sep='\t', comment='#', chunksize=1000000):
        
        # Identify significant snips (cutoff of 5 × 10−8)
        df = df.loc[df['p_value'] < 0.00000005]

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
            # break

    # Open de gencode file (file with position data of snips)
    print(settings.STORAGE['STATIC_FILES']['GENE_POSITIONS'])
    gene_store = pd.HDFStore(settings.STORAGE['STATIC_FILES']['GENE_POSITIONS'])
    genes = []
    num_genes = 0

    # Retrieve all genes that are within the windows
    for chromosome, windows in window_dict.items():
        for window in windows:
            c = f'chr{chromosome}'
            data = gene_store.select('genes', where='chromosome == c and start >= window[0] and end <= window[1]')
            num_genes += len(df)
            update_status(uuid, f'{num_genes} genes found in windows of significant snips')
            genes.append(data)


    gene_store.close()

    update_status(uuid, f'A total of {num_genes} was found in the windows of the significant snips')

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

    gwas_gene_dict = {}
    num_gwas_snips = 0

    # Loop throug GWAS again, now that we know the location of the genes
    for df in pd.read_csv(file_path, sep='\t', comment='#', chunksize=1000000):

        # Do some data processing
        df['varbeta'] = df['standard_error']**2
        # df.loc[df['variant_id'].isna(), "variant_id"] = pd.util.testing.rands_array(10, sum(df['variant_id'].isnull()))
        df.rename(columns={'variant_id': 'snp'}, inplace=True)

        # Calculate -log10(p) for snps that have a p-value greater then 0
        df_p_not_zero = df.loc[df['p_value'] != 0]
        df_p_not_zero['logp'] = -np.log10(df_p_not_zero['p_value'])

        # Set -log10(p) to a very low default value for snps that have a p-value of exactly 0
        df_p_zero = df.loc[df['p_value'] == 0]
        df_p_zero['logp'] = 1000 * 10**-8

        # Concatenate dataframes back together
        df = pd.concat([df_p_not_zero, df_p_zero], ignore_index=True)

        print(df)
        
        # Loop through all genes
        for gene in gene_ids:
            
            if gene_positions_dict[gene][2] != '':
                chromosome = int(gene_positions_dict[gene][2])
            else:
                chromosome = gene_positions_dict[gene][2]

            # Select all snips from GWAS data for every gene
            data = df.loc[
                    (df['base_pair_location'] >= int(gene_positions_dict[gene][0]) - 1000000) & 
                    (df['base_pair_location'] <= int(gene_positions_dict[gene][1]) + 1000000) &
                    (df['chromosome'] == chromosome)
                ]
            
            # Add the data to a dict, with gene ids as key
            if len(data) > 0:
                try:
                    gwas_gene_dict[gene].append(data)
                except:
                    gwas_gene_dict[gene] = data

                # Remove duplicate snips (keep the first for now)
                gwas_gene_dict[gene].drop_duplicates(subset='snp', keep='first', inplace=True)

            num_gwas_snips += len(data)
            update_status(uuid, f'A total of {num_gwas_snips} snips selected for {len(gwas_gene_dict)} genes from GWAS data')

            # break

    eqtl_gene_dict = {}

    # For every gene, retrieve the snips from eQTLGEN dataset
    for gene, row in gwas_gene_dict.items():
        print(gene)
        data = pd.read_hdf(settings.STORAGE['STATIC_FILES']['EQTLGEN'], key='eqtls', where='Gene == gene')
        gwas_snips = row['snp'].tolist()
        data = data[data['snp'].isin(gwas_snips)]
        eqtl_gene_dict[gene] = data
        update_status(uuid, f'{len(data)} snips for gene {gene} retrieved from the eQTLGen dataset ({len(eqtl_gene_dict)}/{len(gwas_gene_dict)})')
        # break 

    os.mkdir(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid)

    written_genes = []

    # Finally, loop through all the genes in the GWAS gene dict that is just filled
    for k, v in gwas_gene_dict.items():

        eqtl_snips = eqtl_gene_dict[k]['snp'].tolist()
        v = v.loc[v['snp'].isin(eqtl_snips)]

        # If the length of the data in the GWAS dict and the eQTLGEN dict is not 0, write them to a JSON file, which can later be loaded in the Coloc R script
        if len(v) > 0 and len(eqtl_gene_dict[k]) > 0:
            location_dict = {
                'chromosome': gene_positions_dict[k][2],
                'start_position': gene_positions_dict[k][0],
                'stop_position': gene_positions_dict[k][1],
                'gene_name': gene_positions_dict[k][3],
            }

            d = {'location': location_dict, 'gwas': v.to_dict(orient='list'), 'eqtls': eqtl_gene_dict[k].to_dict(orient='list')}

            gene_file_name = k + '.json'
            with open(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid / gene_file_name, 'w') as out_file:
                out_file.write(json.dumps(d))
                written_genes.append(k)
                update_status(uuid, f'{len(written_genes)} overlapping genes found and prepared for colocalization analysis')


    os.mkdir(settings.STORAGE['OUT_FILES_PATH'] / uuid)

    genes_analysed = 0

    for gene in written_genes:
        subprocess.call([settings.R['R_PATH'], settings.R['COLOC_PATH'], gene, uuid])
        genes_analysed += 1
        update_status(uuid, f'Colocalization analysis performed for {genes_analysed}/{len(written_genes)} genes')

    
    update_status(uuid, f'Preparing output file')

    directory = settings.STORAGE['OUT_FILES_PATH'] / uuid

    res_dict = {}

    for filename in os.listdir(directory):
        f = os.path.join(directory, filename)
        gene = f.split('\\').pop().split('.')[0]

        with open(f, 'r') as input:
            data = json.load(input)
            res_dict[gene] = data['datasets']


    with open(settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json', 'w') as output:
        update_status(uuid, f'Done!', True)
        json.dump(res_dict, output)