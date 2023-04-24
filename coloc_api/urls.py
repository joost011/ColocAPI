# todo/todo/urls.py : Main urls.py
from django.urls import path, include
from .views import (
    FileView,
    ColocView,
)

urlpatterns = [
    path('file', FileView.as_view()),
    path('coloc', ColocView.as_view()),
]