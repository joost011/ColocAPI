from django.conf import settings
import json
import pandas as pd

class ExportService:

    @staticmethod
    def exportToCSV(uuid):
        file_directory = settings.STORAGE['OUT_FILES_PATH'] / uuid
        file_path = file_directory / 'output.json'
        csv_file_name = uuid + '.csv'
        csv_file_path = file_directory / csv_file_name

        with open(file_path, 'r') as in_file:
            data = json.load(in_file)

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