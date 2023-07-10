from rest_framework import serializers
from .models import ColocAnalysis, ColocAnalysisStatus


class ColocAnalysisStatusSerializer(serializers.ModelSerializer):
    class Meta:
        model = ColocAnalysisStatus
        fields = '__all__'


class ColocAnalysisSerializer(serializers.ModelSerializer):
    status_list = ColocAnalysisStatusSerializer(many=True, read_only=True)

    class Meta:
        model = ColocAnalysis
        fields = ['uuid', 'extension', 'finished', 'status_list']
