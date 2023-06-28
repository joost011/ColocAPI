import pandas as pd
import sys
import os
import json

def read_json(file_path):
    with open(file_path) as file:
        data = json.load(file)
    
    return data['datasets']

def write_file(file_path, content_list):
    with open(file_path, 'w') as outfile:
        for item in content_list:
            outfile.write(item + '\n')

def create_ld_matrix(gene):
    plink = 'C:\\Users\\jobma\\Documents\\School\\DSLS\\semester_2\\integrated_omics\\tests\\susie2\\plink\\plink'
    reference_panel = '..\\susie\\annotated_1000g\\t'
    
    os.system(f'{plink} --bfile {reference_panel} --extract {gene}_stripped.txt --out {gene}_ld --r square  --keep-allele-order --write-snplist')
    
    LD = pd.read_table(f'{gene}_ld.ld', sep='\t', header=None)
    with open(f'{gene}_ld.snplist') as snp_file:
        snp_names = snp_file.read().split('\n')[:-1]
    
    LD.columns = snp_names
    LD.index = snp_names
    
    return LD
    
def get_snp_list(file_path):
    with open(file_path) as file:
        sorter = file.read().splitlines()
        
    return sorter 

def update_df(df, column, snp_list):
    # remove non-intersecting snp_list
    df = df[df[column].isin(snp_list)]
    # sort list according to snp_list
    df = df.sort_values(by=column, key=lambda col: col.map(lambda snp: snp_list.index(snp)))
    
    return df

def has_na_values(matrix):
    has_na = None
    number_na = matrix.isna().sum().sum()
    
    if number_na > 0:
        has_na = True
    else:
        has_na = False
        
    return has_na

def remove_na_values(matrix):
    contains_na = has_na_values(matrix)
    dropped_ids = []
    
    while contains_na:
        id = matrix.isna().sum().sort_values().keys()[-1]  # get id with most na's
        dropped_ids.append(id)
        del matrix[id]
        matrix.drop(index=id, inplace=True)
        contains_na = has_na_values(matrix)
    
    return matrix, dropped_ids 

def drop_ids(df, drop_ids):
    for id in drop_ids:
        df.drop(df[df.snp == id].index, inplace=True)
        
    return df

def write_output(file_path, data):
    with open(file_path, 'w') as out_file:
        out_file.write(data)
    

if __name__ == '__main__':
    # script from Harm-Jan
    json_path = sys.argv[1]
    gene = json_path.rsplit('.', maxsplit=1)[0]
    
    # create df's 
    data = read_json(json_path)
    gwas_df = pd.DataFrame(data['gwas'])
    eqtl_df = pd.DataFrame(data['eqtls'])
    print('got data')
    
    # get snps
    snps = gwas_df['snp'].values.tolist()
    write_file(f'{gene}_stripped.txt', snps)
    print('got SNPs')
    
    # create LD matrix
    print('Creating LD matrix...')
    LD = create_ld_matrix(gene)
    print('got LD matrix')
    
    print('updating data...')
    # get remaining snps after LD --extract 
    sorter = get_snp_list(f'{gene}_ld.snplist')

    # sort the dfs based on the sorter list
    gwas_df = update_df(gwas_df, 'snp', sorter)
    eqtl_df = update_df(eqtl_df, 'snp', sorter)
    
    # drop nans
    LD, dropped_rsids = remove_na_values(LD)
    print(gwas_df.shape)
    gwas_df = drop_ids(gwas_df, dropped_rsids)
    eqtl_df = drop_ids(eqtl_df, dropped_rsids)
    print('updated data')
    print(gwas_df.shape)
    
    # write filtered and sorted list to files
    data['gwas'] = gwas_df.to_dict(orient='list')
    data['eqtls'] = eqtl_df.to_dict(orient='list')
    print('writing files...')
    write_output(f'{gene}_processed.json', json.dumps(data))
    write_output(f'{gene}_processed.snplist', '\n'.join(gwas_df.snp.to_list()))
    write_output(f'{gene}_processed.ld', LD.to_csv(sep='\t', index=False, header=False, lineterminator='\r\n'))
    print('done.')