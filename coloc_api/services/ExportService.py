from django.conf import settings
import json

class ExportService:

    @classmethod
    def exportAllToCSV(cls, uuid):
        # Load data from file
        data = cls.__load_data_from_file(uuid)

        # Create output file name and path
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        csv_file_name = uuid + '.csv'
        csv_file_path = file_directory / csv_file_name

        rows = [
            ','.join(['Gene name', 'Gene ID', 'Chromosome', 'Tested SNPs', 'PPH4', '\n'])
        ]
        
        for k, v in data['genes'].items():
            gene_name = v['meta_data']['gene_name']
            gene_id = k
            chromosome = v['meta_data']['chromosome']
            tested_snps = v['meta_data']['total_intersecting_snps']
            pph4 = v['posterior']['PP.H4.abf']
            row = [gene_name, gene_id, str(chromosome), str(tested_snps), str(pph4), '\n']
            rows.append(','.join(row))


        with open(csv_file_path, 'w') as out_file:
            out_file.writelines(rows)

        return csv_file_path
    

    @classmethod
    def exportGeneToCSV(cls, uuid, gene):
        # Load data from file
        data = cls.__load_data_from_file(uuid)

        # Create output file name and path
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        csv_file_name = f'{gene}_{uuid}.csv'
        csv_file_path = file_directory / csv_file_name

        # Format data
        columns = ['SNP,Position,PPH4,\n']
        snp_list = data['genes'][gene]['gwas']['snp']
        location_list = data['genes'][gene]['eqtls']['SNPPos']
        pph4_list = data['genes'][gene]['coloc']
        csv_data = [','.join([str(snp_list[i]), str(location_list[i]), str(pph4_list[i]), '\n']) for i in range(len(snp_list))]
        rows = columns + csv_data

        with open(csv_file_path, 'w') as out_file:
            out_file.writelines(rows)

        return csv_file_path
    
    def __load_data_from_file(uuid):
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        file_path = file_directory / 'output.json'

        with open(file_path, 'r') as in_file:
            data = json.load(in_file)

        return data