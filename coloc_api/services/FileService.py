from django.core.files.storage import default_storage
from django.conf import settings
import os

class FileService:

    @staticmethod
    def save_input_file(file, file_name):
        '''Function that stores an input file in the input file directory'''
        default_storage.save(settings.STORAGE['IN_FILES_PATH'] / file_name, file)

    @staticmethod
    def delete_input_file(file_name):
        '''Function that deletes an input file from the input file directory'''
        default_storage.delete(settings.STORAGE['IN_FILES_PATH'] / file_name)

    @staticmethod
    def delete_processed_files(uuid):
        '''Function that deletes all processed files for a specific uuid'''
        directory = settings.STORAGE['PROCESSED_FILES_PATH'] / uuid
        for file_name in os.listdir(directory):
            file_path = os.path.join(directory, file_name)
            default_storage.delete(file_path)

    @staticmethod
    def delete_processed_file_directory(uuid):
        '''Function that deletes a folder with the specified uuid from the processed files directroy'''
        default_storage.delete(settings.STORAGE['PROCESSED_FILES_PATH'] / uuid)