# Generated by Django 4.2.4 on 2023-08-24 12:23

from django.conf import settings
from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        migrations.swappable_dependency(settings.AUTH_USER_MODEL),
        ('currant', '0012_userextradata_alter_operation_operation_type'),
    ]

    operations = [
        migrations.AddField(
            model_name='userextradata',
            name='user',
            field=models.OneToOneField(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='extra', to=settings.AUTH_USER_MODEL),
        ),
    ]
