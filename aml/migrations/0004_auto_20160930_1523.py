# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-09-30 15:23
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('aml', '0003_auto_20160628_1208'),
    ]

    operations = [
        migrations.CreateModel(
            name='FragmentAnalysis',
            fields=[
                ('result_id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('sample', models.CharField(default='Unknown', max_length=100)),
                ('run', models.CharField(default='Unknown', max_length=100)),
                ('itd', models.IntegerField(default=0)),
                ('ab', models.FloatField(default='0.0')),
            ],
            options={
                'db_table': 'Frag',
            },
        ),
        migrations.CreateModel(
            name='Variants',
            fields=[
                ('result_id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('chrom', models.CharField(default='00', max_length=6)),
                ('pos', models.IntegerField(default=0)),
                ('ref', models.CharField(default='.', max_length=300)),
                ('alt', models.CharField(default='.', max_length=300)),
            ],
            options={
                'db_table': 'Variants',
            },
        ),
        migrations.AddField(
            model_name='results',
            name='precision',
            field=models.CharField(default='', max_length=10),
        ),
    ]