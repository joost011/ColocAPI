from django.core.files.storage import default_storage
from rest_framework.response import Response
from rest_framework import status
from rest_framework.views import APIView
import json
import uuid

class FileView(APIView):

    def post(self, request, *args, **kwargs):
        file = request.data['file']
        file_extension = str(file).split('.')[-1]
        file_name = f'{uuid.uuid4()}.{file_extension}'
        print('file', file_name)
        default_storage.save(f'files/in_files/{file_name}', file)
        return Response(json.loads(json.dumps({'file_name': file_name})), status=status.HTTP_200_OK)
