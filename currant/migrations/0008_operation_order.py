# Generated by Django 4.2.4 on 2023-08-18 14:07

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('currant', '0007_operation_operation_type'),
    ]

    operations = [
        migrations.AddField(
            model_name='operation',
            name='order',
            field=models.IntegerField(default=0),
        ),
    ]
