from django.urls import re_path

from currantDjango.consumers import OperationConsumer

websocket_urlpatterns = [
    re_path(r'ws/operation/(?P<session_id>\w+)/(?P<personal_id>\w+)/$', OperationConsumer.as_asgi()),
]