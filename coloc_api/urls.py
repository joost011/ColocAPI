from django.urls import path
from .views import (
    FileView,
    ColocView,
    ResultView,
)

urlpatterns = [
    path('file', FileView.as_view()),
    path('coloc', ColocView.as_view()),
    path('coloc/<str:uuid>', ColocView.as_view()),
    path('result/<str:uuid>', ResultView.as_view())
]