#!/usr/bin/python3

import pandas as pd
import json
import sys

def extract_data(file_path):
    """Loads JSON data from a file"""
    with open(file_path) as file:
        data = json.load(file)
    
    return data
    
def get_sorter_list(file_path):
    """Loads the sorter list that is obtained from PLINK"""
    with open(file_path) as file:
        sorter = file.read().splitlines()
        
    return sorter   

def sort_df(df, column, sorter):
    """Sorts a dataframe based on the ocurrence of the sorter"""
    df.sort_values(by=column, key=lambda col: col.map(lambda snp: sorter.index(snp)), inplace=True)
    json_sorted = {col: df[col].values.tolist() for col in df}
    
    return json_sorted

def write_output(file_path, data):
    """Writes the output in JSON format to a file"""
    with open(file_path, 'w') as out_file:
        out_file.write(json.dumps(data))


if __name__ == '__main__':
    input_file = sys.argv[1]
    sorter_file = sys.argv[2]
    output_file = '_sorted.'.join(input_file.rsplit('.', maxsplit=1))
    
    # get the data which needs to be sorted
    data = extract_data(input_file)
    gwas_df = pd.DataFrame(data['gwas'])
    eqtl_df = pd.DataFrame(data['eqtls'])
    
    # get the sort data as a list
    sorter = get_sorter_list(sorter_file)

    # remove non-intersecting observations
    gwas_df = gwas_df[gwas_df['snp'].isin(sorter)]
    eqtl_df = eqtl_df[eqtl_df['snp'].isin(sorter)]
    
    # sort the dfs based on the sorter list
    data['gwas'] = sort_df(gwas_df, 'snp', sorter)
    data['eqtls'] = sort_df(eqtl_df, 'snp', sorter)
    
    # write the output in JSON format
    write_output(output_file, data)