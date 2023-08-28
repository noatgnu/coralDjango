from rest_flex_fields import FlexFieldsModelSerializer
from rest_framework import serializers
from json import loads
from currant.models import File, Operation, OperationSequence, Coral


class FileSerializer(serializers.ModelSerializer):
    columns = serializers.SerializerMethodField()
    file_type = serializers.SerializerMethodField()

    def get_columns(self, obj):
        return loads(obj.columns)

    def get_file_type(self, obj):
        return obj.file_type.split(";")

    class Meta:
        model = File
        fields = ["id", "link_id", "file", "columns", "file_type"]
        lookup_field = "link_id"


class OperationSerializer(FlexFieldsModelSerializer):
    input_files = FileSerializer(many=True)
    output_files = FileSerializer(many=True)
    value = serializers.SerializerMethodField()

    def get_value(self, obj):
        return loads(obj.value)

    class Meta:
        model = Operation
        fields = ["id", "name", "value", "created", "input_files", "output_files", "sequence", "job_id", "job_finished", "order", "operation_type" ]


class OperationSequenceSerializer(serializers.ModelSerializer):
    class Meta:
        model = OperationSequence
        fields = ["name", "value", "created", "operations"]


class CoralSerializer(FlexFieldsModelSerializer):
    class Meta:
        modal = Coral
        fields = ["name", "description", "link_id", "created", "owners", "files", "operations", "operation_sequences"]
        lookup_field = "link_id"