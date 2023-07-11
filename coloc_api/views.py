from rest_framework.response import Response
from django.http import FileResponse
from rest_framework import status
from rest_framework.views import APIView
import json
import uuid
from .models import ColocAnalysis, ColocAnalysisStatus
from .serializers import ColocAnalysisSerializer
from django.conf import settings
from .jobs import initialize_coloc
from .services.ExportService import ExportService
from .services.FileService import FileService
import django_rq


class FileView(APIView):
    '''View that handles the file related requests'''

    def post(self, request, *args, **kwargs):
        '''Function that handles the post request for uploading a GWAS file'''
        file = request.data['file']
        file_extension = '.'.join(str(file).split('.')[1:])
        file_name = f'{uuid.uuid4()}.{file_extension}'
        FileService.save_input_file(file, file_name)

        # Create database colocalization object with initial status object
        file_uuid = file_name.split('.')[0]
        file_extension = '.'.join(file_name.split('.')[1:])
        coloc_analysis_object = ColocAnalysis.objects.create(uuid=file_uuid, extension=file_extension, finished=False)
        ColocAnalysisStatus.objects.create(coloc_analysis=coloc_analysis_object, status_message='Initializing colocalization analysis', status_order=0)

        return Response(json.loads(json.dumps({'file_name': file_name})), status=status.HTTP_200_OK)
    

class ColocView(APIView):
    '''View that handels the coloc related requests'''

    def post(self, request, *args, **kwargs):
        '''Function that handles the post request for starting an analysis'''
        coloc_type = request.data['type']

        if coloc_type == 'cc':
            coloc_args = (
                request.data['file'],
                coloc_type,
                request.data['numCases'],
                request.data['numControls'],
            )
        
        else:
            coloc_args = (
                request.data['file'],
                coloc_type,
                request.data['sampleSize'],
            )

        print(request.data)

        django_rq.enqueue(initialize_coloc, job_timeout=604800, args=coloc_args)
        return Response(json.dumps({'success': True}), status=status.HTTP_200_OK) 
    
    def get(self, request, *args, **kwargs):
        '''Function that handles the get request for getting the status of an analysis process'''
        uuid = kwargs['uuid']
        analysis = ColocAnalysis.objects.prefetch_related('status_list').get(uuid=uuid)
        serializer = ColocAnalysisSerializer(analysis)
        return Response(serializer.data, status=status.HTTP_200_OK)


class ResultView(APIView):
    '''View that handles the result related requests'''
    
    def get(self, request, *args, **kwargs):
        '''Function that handles the get request for getting results of an analysis'''
        uuid = kwargs['uuid']
        output_path = settings.STORAGE['OUT_FILES_PATH'] / uuid / 'output.json'
        with open(output_path, 'r') as file:
            data = json.load(file)
        return Response(data, status=status.HTTP_200_OK) 


class ExportView(APIView):
    '''View that handles the export related requests'''

    def post(self, request, *args, **kwargs):
        '''Function that handles the post request for creating an export of an analysis'''
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