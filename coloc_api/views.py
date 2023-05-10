from django.core.files.storage import default_storage
from rest_framework.response import Response
from rest_framework import status
from rest_framework.views import APIView
from rest_framework import serializers
import json
import uuid
import sqlite3
from .functions import coloc
from .models import ColocAnalysis
from .serializers import ColocAnalysisSerializer
class FileView(APIView):

    def post(self, request, *args, **kwargs):
        file = request.data['file']
        file_extension = '.'.join(str(file).split('.')[1:])
        file_name = f'{uuid.uuid4()}.{file_extension}'
        default_storage.save(f'storage/in_files/{file_name}', file)

        # Create database colocalization object
        conn = sqlite3.connect('db.sqlite3')
        cursor = conn.cursor()
        cursor.execute(f'''INSERT INTO coloc_api_colocanalysis(
            uuid, extension, status_message, finished) VALUES 
            ('{file_name.split('.')[0]}', '{'.'.join(file_name.split('.')[1:])}', 'Colozalization analysis started', False)''')
        conn.commit()
        conn.close()

        return Response(json.loads(json.dumps({'file_name': file_name})), status=status.HTTP_200_OK)
    

class ColocView(APIView):

    def post(self, request, *args, **kwarts):
        coloc_output = coloc(request.data['file'])
        return Response(coloc_output, status=status.HTTP_200_OK) 
    
    def get(self, request, *args, **kwargs):
        uuid = kwargs['uuid']
        analysis = ColocAnalysis.objects.get(uuid=uuid)
        serializer = ColocAnalysisSerializer(analysis)
        return Response(json.loads(json.dumps(serializer.data)), status=status.HTTP_200_OK) 
    

class ResultView(APIView):
    
    def get(self, request, *args, **kwargs):
        uuid = kwargs['uuid']
        out_files_dir = 'C:\\Users\\joost\\OneDrive\\Documents\\master\\coloc_project\\django\\coloc\\files\\out_files\\'
        output_path = out_files_dir + uuid + '\\output.json'
        with open(output_path, 'r') as file:
            data = json.load(file)
        # data = json.load(output_path)
        # print(type(data))
        return Response(data, status=status.HTTP_200_OK) 
