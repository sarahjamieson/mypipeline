# -*- coding: utf-8 -*-
# Generated by Django 1.9.2 on 2016-06-27 11:22
from __future__ import unicode_literals

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Results',
            fields=[
                ('result_id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('run', models.CharField(default='Unknown', max_length=100)),
                ('sample', models.CharField(default='Unknown', max_length=100)),
                ('caller', models.CharField(default='Unknown', max_length=20)),
                ('chrom', models.CharField(default='00', max_length=6)),
                ('pos', models.IntegerField(default=0)),
                ('ref', models.CharField(default='.', max_length=300)),
                ('alt', models.CharField(default='.', max_length=300)),
                ('chr2', models.CharField(default='00', max_length=6)),
                ('end', models.IntegerField(default=0)),
                ('region', models.CharField(default='.', max_length=20)),
                ('sv_type', models.CharField(default='Unknown', max_length=15)),
                ('size', models.IntegerField(default=0)),
                ('gt', models.CharField(default='.', max_length=5)),
                ('depth', models.IntegerField(default=0)),
                ('alleles', models.CharField(default='0,0', max_length=10)),
                ('ab', models.FloatField(default='0.0')),
                ('gene', models.CharField(default='Unknown', max_length=50)),
                ('func', models.CharField(default='Unknown', max_length=200)),
                ('exonic_func', models.CharField(default='Unknown', max_length=200)),
            ],
            options={
                'db_table': 'Results',
            },
        ),
        migrations.CreateModel(
            name='Runs',
            fields=[
                ('run_id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('run', models.CharField(default='Unknown', max_length=100)),
            ],
        ),
        migrations.CreateModel(
            name='Samples',
            fields=[
                ('sample_id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('sample', models.CharField(default='Unknown', max_length=100)),
                ('run', models.CharField(default='Unknown', max_length=100)),
            ],
        ),
    ]
