# todo/todo/urls.py : Main urls.py
from django.urls import path, include
from .views import (
    FileView,
)

urlpatterns = [
    path('file', FileView.as_view()),
]