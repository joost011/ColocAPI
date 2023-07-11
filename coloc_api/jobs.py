from django_rq import job
from .services.ColocService import ColocService

@job('high')
def initialize_coloc(file_name):
    '''Function that starts a coloc process on the Redis server'''
    coloc_service = ColocService()
    coloc_service.start_coloc(file_name)
