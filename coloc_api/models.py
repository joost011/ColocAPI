from django.db import models

class ColocAnalysis(models.Model):
    uuid = models.CharField(max_length=200)
    extension = models.CharField(max_length=200)
    status_message = models.CharField(max_length=500)
    finished = models.BooleanField(default=False)