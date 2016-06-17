from __future__ import unicode_literals

from django.db import models
import django_tables2 as table


class Results(models.Model):
    sample_id = models.AutoField(primary_key=True, unique=True)
    run = models.CharField(max_length=100, default='Unknown')
    sample = models.CharField(max_length=100, default='Unknown')
    caller = models.CharField(max_length=20, default='Unknown')
    chrom = models.CharField(max_length=6, default='00')
    pos = models.IntegerField(default=0)
    ref = models.CharField(max_length=300, default='.')
    alt = models.CharField(max_length=300, default='.')
    chr2 = models.CharField(max_length=6, default='00')
    end = models.IntegerField(default=0)
    region = models.CharField(max_length=20, default='.')
    sv_type = models.CharField(max_length=15, default='Unknown')
    size = models.IntegerField(default=0)
    gt = models.CharField(max_length=5, default='.')
    depth = models.IntegerField(default=0)
    alleles = models.CharField(max_length=10, default='0,0')
    ab = models.FloatField(default='0.0')
    gene = models.CharField(max_length=50, default='Unknown')
    func = models.CharField(max_length=200, default='Unknown')
    exonic_func = models.CharField(max_length=200, default='Unknown')

    class Meta:
        app_label = 'aml'
        db_table = 'Results'


class ResultTable(table.Table):
    class Meta:
        model = Results
