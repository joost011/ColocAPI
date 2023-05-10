from rest_framework import serializers

class ColocAnalysisSerializer(serializers.Serializer):
    uuid = serializers.CharField(max_length=200)
    extension = serializers.CharField(max_length=200)
    status_message = serializers.CharField(max_length=500)
    finished = serializers.BooleanField()