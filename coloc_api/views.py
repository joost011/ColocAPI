from django.core.files.storage import default_storage
from rest_framework.response import Response
from django.http import FileResponse
from django.db.models import Prefetch
from rest_framework import status
from rest_framework.views import APIView
from rest_framework import serializers
import json
import uuid
import sqlite3
from .functions import coloc
from .models import ColocAnalysis, ColocAnalysisStatus
from .serializers import ColocAnalysisSerializer
from django.conf import settings
from .services.ExportService import ExportService
import django_rq

class FileView(APIView):

    def post(self, request, *args, **kwargs):
        file = request.data['file']
        file_extension = '.'.join(str(file).split('.')[1:])
        file_name = f'{uuid.uuid4()}.{file_extension}'
        default_storage.save(f'storage/in_files/{file_name}', file)

        # Create database colocalization object with initial status object
        file_uuid = file_name.split('.')[0]
        file_extension = '.'.join(file_name.split('.')[1:])
        coloc_analysis_object = ColocAnalysis.objects.create(uuid=file_uuid, extension=file_extension, finished=False)
        ColocAnalysisStatus.objects.create(coloc_analysis=coloc_analysis_object, status_message='Initializing colocalization analysis', status_order=0)

        return Response(json.loads(json.dumps({'file_name': file_name})), status=status.HTTP_200_OK)
    

class ColocView(APIView):

    def post(self, request, *args, **kwargs):
        django_rq.enqueue(coloc, job_timeout=604800, args=(request.data['file'],))
        return Response(json.dumps({'success': True}), status=status.HTTP_200_OK) 
    
    def get(self, request, *args, **kwargs):
        uuid = kwargs['uuid']
        analysis = ColocAnalysis.objects.prefetch_related('status_list').get(uuid=uuid)
        serializer = ColocAnalysisSerializer(analysis)
        return Response(serializer.data, status=status.HTTP_200_OK) 
    

class ResultView(APIView):
    
    def get(self, request, *args, **kwargs):
        uuid = kwargs['uuid']
        output_path = settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json'
        with open(output_path, 'r') as file:
            data = json.load(file)
        return Response(data, status=status.HTTP_200_OK) 


class ExportView(APIView):

    def post(self, request, *args, **kwargs):
        uuid = request.data['uuid']
        export_type = request.data['export_type']

        if export_type == 'csv_all':
            export_file_path = ExportService.exportAllToCSV(uuid)

        if export_type == 'csv_gene':
            gene = request.data['gene']
            export_file_path = ExportService.exportGeneToCSV(uuid, gene)

        file = open(export_file_path, 'rb')
        response = FileResponse(file, content_type='application/json')
        response['Content-Disposition'] = f'attachment; filename="{uuid}.csv"'
        return response