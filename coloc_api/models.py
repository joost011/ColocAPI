from django.db import models

class ColocAnalysis(models.Model):
    uuid = models.CharField(max_length=200, primary_key=True, unique=True)
    extension = models.CharField(max_length=200)
    finished = models.BooleanField(default=False)

class ColocAnalysisStatus(models.Model):
    coloc_analysis = models.ForeignKey(
        ColocAnalysis, 
        on_delete=models.CASCADE, 
        to_field='uuid',
        related_name='status_list')
    
    status_message = models.CharField(max_length=500)
    status_order = models.SmallIntegerField()