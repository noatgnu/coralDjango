import uuid

from django.db import models
from django.utils import timezone


# Create your models here.

class UserExtraData(models.Model):
    userSource = models.TextField()
    userSourceId = models.TextField()
    user = models.OneToOneField('auth.User', related_name='extra', on_delete=models.CASCADE, blank=True, null=True)


class File(models.Model):
    root = models.BooleanField(default=False)
    created = models.DateTimeField(default=timezone.now, editable=False)
    link_id = models.TextField(unique=True, default=uuid.uuid4, editable=False)
    file = models.FileField(upload_to='currant/files/', blank=True, null=True)
    file_type = models.TextField(blank=True, null=True)
    owners = models.ManyToManyField('auth.User', related_name='files', blank=True)
    columns = models.TextField(blank=True, null=True)
    coral = models.ForeignKey('Coral', related_name='files', on_delete=models.RESTRICT, blank=True, null=True)

    def delete(self, using=None, keep_parents=False):
        self.file.delete()
        super().delete(using=using, keep_parents=keep_parents)

class Operation(models.Model):
    name = models.TextField()
    value = models.TextField()
    created = models.DateTimeField(default=timezone.now, editable=False)
    input_files = models.ManyToManyField(File, related_name='input_parameters', blank=True)
    output_files = models.ManyToManyField(File, related_name='output_parameters', blank=True)
    sequence = models.ForeignKey('OperationSequence', related_name='operations', on_delete=models.RESTRICT, blank=True, null=True)
    operation_type_choices = [
        ("APS-DF", "AlphaPeptStats"),
        ("PLCHD", "Placeholder"),
        ("RQF-PROT", "QFeatures (Protein) (R)"),
        ("RQF-PEP", "QFeatures (Peptide) (R)"),
        ("RQF-NORM", "QFeatures (Normalization) (R)"),
        ("RQF-IMP", "QFeatures (Imputation) (R)"),
        ("R-CORR", "Correlation Matrix (R)"),
    ]

    operation_type = models.CharField(
        choices=operation_type_choices,
        default="PLCHD",
    )

    order = models.IntegerField(default=0)
    job_id = models.TextField(blank=True, null=True)
    job_finished = models.BooleanField(default=False)
    coral = models.ForeignKey('Coral', related_name='operations', on_delete=models.RESTRICT, blank=True, null=True)
    job_error = models.TextField(blank=True, null=True)
    job_error_status = models.BooleanField(default=False)
    result = models.TextField(blank=True, null=True)

class OperationSequence(models.Model):
    name = models.TextField()
    value = models.TextField()
    created = models.DateTimeField(default=timezone.now, editable=False)
    coral = models.ForeignKey('Coral', related_name='operation_sequences', on_delete=models.RESTRICT, blank=True, null=True)


class Coral(models.Model):
    name = models.TextField()
    description = models.TextField()
    link_id = models.TextField(unique=True, default=uuid.uuid4, editable=False)
    created = models.DateTimeField(default=timezone.now, editable=False)
    owners = models.ManyToManyField('auth.User', related_name='corals', blank=True)