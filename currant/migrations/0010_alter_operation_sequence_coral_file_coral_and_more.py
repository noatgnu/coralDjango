# Generated by Django 4.2.4 on 2023-08-21 10:37

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion
import django.utils.timezone
import uuid


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('currant', '0009_operation_job_finished_operation_job_id'),
    ]

    operations = [
        migrations.AlterField(
            model_name='operation',
            name='sequence',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='operations', to='currant.operationsequence'),
        ),
        migrations.CreateModel(
            name='Coral',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.TextField()),
                ('description', models.TextField()),
                ('link_id', models.TextField(default=uuid.uuid4, editable=False, unique=True)),
                ('created', models.DateTimeField(default=django.utils.timezone.now, editable=False)),
                ('owners', models.ManyToManyField(blank=True, related_name='corals', to=settings.AUTH_USER_MODEL)),
            ],
        ),
        migrations.AddField(
            model_name='file',
            name='coral',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='files', to='currant.coral'),
        ),
        migrations.AddField(
            model_name='operation',
            name='coral',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='operations', to='currant.coral'),
        ),
        migrations.AddField(
            model_name='operationsequence',
            name='coral',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.RESTRICT, related_name='operation_sequences', to='currant.coral'),
        ),
    ]
