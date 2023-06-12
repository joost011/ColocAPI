from django.conf import settings
import json
import pandas as pd

class ExportService:

    @classmethod
    def exportAllToCSV(cls, uuid):
        # Load data from file
        data = cls.__load_data_from_file(uuid)

        # Create output file name and path
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        csv_file_name = uuid + '.csv'
        csv_file_path = file_directory / csv_file_name

        # Create empty dataframe
        df = pd.DataFrame(columns=['gene_name', 'gene_id', 'chromosome', 'tested_snps', 'pph4', 'top_snp'])

        # Format every gene in output file and add to the empty dataframe
        for k, v in data['genes'].items():
            gene_name = v['meta_data']['gene_name']
            gene_id = k
            chromosome = v['meta_data']['chromosome']
            tested_snps = v['meta_data']['total_intersecting_snps']
            pph4 = v['posterior']['PP.H4.abf']
            top_snp = v['meta_data']['top_snp']
            new_row = pd.DataFrame({'gene_name': [gene_name], 'gene_id': [gene_id], 'chromosome': [chromosome], 'tested_snps': [tested_snps], 'pph4': [pph4], 'top_snp': [top_snp]})
            df = pd.concat([df, new_row], ignore_index=True)

        # Sort dataframe by posterior h4 probability
        df = df.sort_values('pph4', ascending=False)

         # Write dataframe dataframe to file
        df.to_csv(csv_file_path, index=False)

        return csv_file_path
    

    @classmethod
    def exportGeneToCSV(cls, uuid, gene):
        # Load data from file
        data = cls.__load_data_from_file(uuid)
        data = data['genes'][gene]

        # Create output file name and path
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        csv_file_name = f'{gene}_{uuid}.csv'
        csv_file_path = file_directory / csv_file_name

        # Process GWAS data
        gwas_data = data['gwas']
        gwas_df = pd.DataFrame.from_dict(gwas_data)
        gwas_df = gwas_df[['snp', 'p_value']]
        gwas_df.rename(columns={'p_value': 'gwas_p_value'}, inplace=True)

        # Process eQTL data
        eqtl_data = data['eqtls']
        eqtl_df = pd.DataFrame.from_dict(eqtl_data)
        eqtl_df = eqtl_df[['snp', 'pvalues', 'Zscore', 'SNPPos']]
        eqtl_df.rename(columns={'pvalues': 'eqtl_p_value', 'Zscore': 'eqtl_z_score', 'SNPPos': 'snp_position'}, inplace=True)

        # Process coloc data
        coloc_data = data['coloc']
        coloc_df = pd.DataFrame.from_dict(coloc_data)
        coloc_df.rename(columns={'SNP.PP.H4': 'pph4'}, inplace=True)

        # Merge dataframes
        merged_df = pd.merge(gwas_df, eqtl_df, on='snp')
        merged_df = pd.merge(merged_df, coloc_df, on='snp')

        # Sort dataframe by posterior h4 probability
        merged_df = merged_df.sort_values('pph4', ascending=False)

        # Write merged dataframe to file
        merged_df.to_csv(csv_file_path, index=False)

        return csv_file_path
    
    
    def __load_data_from_file(uuid):
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        file_path = file_directory / 'output.json'

        with open(file_path, 'r') as in_file:
            data = json.load(in_file)

        return data