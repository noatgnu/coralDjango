import django_rq
import numpy as np
import pandas as pd
from alphastats import GenericLoader, DataSet
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from rest_flex_fields import is_expanded
from rest_framework import viewsets, permissions, status
from django.core.files.base import File as djangoFile, ContentFile
from rest_framework.decorators import action
from rest_framework.parsers import MultiPartParser, JSONParser
from rest_framework.response import Response
from rq.job import Job

from currant.models import File, Operation, OperationSequence, Coral
from currant.serializers import FileSerializer, OperationSerializer, OperationSequenceSerializer, CoralSerializer
from json import dumps

from currant.tasks import save_file, alphapeptstats_diff, r_qfeatures_protein, r_qfeatures_normalization, \
    r_qfeatures_imputation, r_correlation_matrix


class FileViewSets(viewsets.ModelViewSet):
    queryset = File.objects.all()
    serializer_class = FileSerializer
    parser_classes = [MultiPartParser, JSONParser]
    permission_classes = (permissions.AllowAny,)
    lookup_field = 'link_id'


    def get_queryset(self):
        return File.objects.all()

    def create(self, request, *args, **kwargs):
        f = File(file_type="raw_file"+";"+request.data['fileType'])
        channel_layer = get_channel_layer()
        async_to_sync(channel_layer.group_send)("1", {
            'type': 'job_message',
            'message': {
                'message': "Start uploaded",
                'senderName': "Server",
                'requestType': "fileUpload"
            }
        })

        extension = self.request.data['file'].name.split(".")[-1]
        if request.data['fileType'] == "table":
            first_line = self.request.data['file'].readline().decode('utf-8').strip()
            columns = []
            if extension == "csv":
                # df = pd.read_csv(f.file.url)
                columns = first_line.split(",")
            #elif extension == "xlsx":
                # df = pd.read_excel(f.file.url)
            elif extension == "txt" or extension == "tsv":
                columns = first_line.split("\t")
                # df = pd.read_csv(f.file.url, sep="\t")
            #df.replace("#VALUE!", np.NAN, inplace=True)
            if len(columns) > 0:
                f.columns = dumps(columns)
            self.request.data['file'].seek(0)
        f.file.save(f"{str(f.link_id)}.{extension}", djangoFile(self.request.data['file']))

        f.save()
        file_json = FileSerializer(f, many=False, context={"request": request})
        async_to_sync(channel_layer.group_send)("1", {
            'type': 'job_message',
            'message': {
                'message': "Finished uploaded",
                'senderName': "Server",
                'requestType': "fileUpload"
            }
        })
        return Response(file_json.data, status=status.HTTP_201_CREATED)

    def destroy(self, request, *args, **kwargs):
        instance = self.get_object()
        instance.file.delete()
        instance.delete()
        return Response(status=status.HTTP_204_NO_CONTENT)

    @action(detail=True, methods=['post'], permission_classes=[permissions.AllowAny])
    def process(self, request, *args, **kwargs):
        instance = self.get_object()
        return Response(status=status.HTTP_204_NO_CONTENT)


class OperationViewSets(viewsets.ModelViewSet):
    queryset = Operation.objects.all()
    serializer_class = OperationSerializer
    parser_classes = [MultiPartParser, JSONParser]
    permission_classes = (permissions.AllowAny,)

    def get_queryset(self):
        return self.queryset

    def create(self, request, *args, **kwargs):
        ot = request.data["operationType"]
        choices = [c[0] for c in Operation._meta.get_field('operation_type').choices]
        if ot in choices:
            o = Operation(value=dumps(request.data["form"]), operation_type=request.data["operationType"])
            o.save()
            o.input_files.set([File.objects.get(link_id=link_id) for link_id in request.data['inputFiles']])
            o.save()
            if ot == "APS-DF":
                channel_layer = get_channel_layer()
                async_to_sync(channel_layer.group_send)(request.data["sessionId"], {
                    'type': 'job_message',
                    'message': {
                        'message': "Created operation",
                        'senderName': "Server",
                        'requestType': "Differential Analysis (AlphaPeptStats)",
                        'operationId': o.id
                    }
                })
                job = alphapeptstats_diff.delay(o.id, request.data["sessionId"])
                o.job_id = job.id
            elif ot == "RQF-PROT" or ot == "RQF-PEP":
                channel_layer = get_channel_layer()
                async_to_sync(channel_layer.group_send)(request.data["sessionId"], {
                    'type': 'job_message',
                    'message': {
                        'message': "Created operation",
                        'senderName': "Server",
                        'requestType': "QFeatures (Protein) (R)",
                        'operationId': o.id
                    }
                })
                job = r_qfeatures_protein.delay(o.id, request.data["sessionId"])
                o.job_id = job.id
            elif ot == "RQF-NORM":
                channel_layer = get_channel_layer()
                async_to_sync(channel_layer.group_send)(request.data["sessionId"], {
                    'type': 'job_message',
                    'message': {
                        'message': "Created operation",
                        'senderName': "Server",
                        'requestType': "QFeatures (Normalization) (R)",
                        'operationId': o.id
                    }
                })
                job = r_qfeatures_normalization.delay(o.id, request.data["sessionId"])
                o.job_id = job.id
            elif ot == "RQF-IMP":
                channel_layer = get_channel_layer()
                async_to_sync(channel_layer.group_send)(request.data["sessionId"], {
                    'type': 'job_message',
                    'message': {
                        'message': "Created operation",
                        'senderName': "Server",
                        'requestType': "QFeatures (Imputation) (R)",
                        'operationId': o.id
                    }
                })
                job = r_qfeatures_imputation.delay(o.id, request.data["sessionId"])
                o.job_id = job.id
            elif ot == "R-CORR":
                channel_layer = get_channel_layer()
                async_to_sync(channel_layer.group_send)(request.data["sessionId"], {
                    'type': 'job_message',
                    'message': {
                        'message': "Created operation",
                        'senderName': "Server",
                        'requestType': "Correlation Matrix (R)",
                        'operationId': o.id
                    }
                })
                job = r_correlation_matrix.delay(o.id, request.data["sessionId"])
                o.job_id = job.id

            o.save()
            o_json = OperationSerializer(o, many=False, context={"request": request})
            return Response(o_json.data, status=status.HTTP_201_CREATED)

        return Response(status=status.HTTP_400_BAD_REQUEST)


class OperationSequenceViewSets(viewsets.ModelViewSet):
    queryset = OperationSequence.objects.all()
    serializer_class = OperationSequenceSerializer
    parser_classes = [MultiPartParser, JSONParser]
    permission_classes = (permissions.AllowAny,)


class CoralViewSets(viewsets.ModelViewSet):
    queryset = Coral.objects.all()
    serializer_class = CoralSerializer
    parser_classes = [MultiPartParser, JSONParser]
    permission_classes = (permissions.AllowAny,)
    lookup_field = 'link_id'

    def create(self, request, *args, **kwargs):
        c = Coral(request.data["name"], request.data["description"])
        c.save()
        if request.user:
            if request.user.is_authenticated:
                c.owners.set([request.user])
        files = [File.objects.get(link_id=f["link_id"]) for f in request.data["files"]]
        if len(files) > 0:
            c.files.set(files)

        operations = [Operation.objects.get(id=o["id"]) for o in request.data["operations"]]
        if len(operations) > 0:
            c.operations.set(operations)

        operation_sequences = [OperationSequence.objects.get(id=os["id"]) for os in request.data["operation_sequences"]]
        if len(operation_sequences) > 0:
            c.operation_sequences.set(operation_sequences)
        c.save()
        c_json = CoralSerializer(c, many=False, context={"request": request})
        return Response(c_json.data, status=status.HTTP_201_CREATED)