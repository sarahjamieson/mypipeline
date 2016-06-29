from __future__ import unicode_literals

from django.db import models
import django_tables2 as table


class Results(models.Model):
    result_id = models.AutoField(primary_key=True, unique=True)
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
    freq = models.CharField(max_length=30, default='Unknown')
    af = models.FloatField(default='0.0')
    precision = models.CharField(max_length=10, default='')

    class Meta:
        app_label = 'aml'
        db_table = 'Results'


class PindelTable(table.Table):
    result_id = table.Column(visible=False)
    run = table.Column(visible=False)
    sample = table.Column(verbose_name="SAMPLE", default='', visible=False)
    caller = table.Column(verbose_name="CALLER", default='', visible=False)
    chrom = table.Column(verbose_name="CHROM", default='')
    pos = table.Column(verbose_name="POS", default='')
    ref = table.Column(verbose_name="REF", default='')
    alt = table.Column(verbose_name="ALT", default='')
    chr2 = table.Column(verbose_name="CHR2", default='', visible=False)
    end = table.Column(verbose_name="END", default='')
    region = table.Column(verbose_name="REGION", default='')
    sv_type = table.Column(verbose_name="TYPE", default='')
    size = table.Column(verbose_name="SIZE", default='')
    gt = table.Column(verbose_name="GT", default='')
    depth = table.Column(verbose_name="DEPTH", default='')
    alleles = table.Column(verbose_name="ALLELES", default='')
    ab = table.Column(verbose_name="AB", default='')
    gene = table.Column(verbose_name="GENE", default='')
    func = table.Column(verbose_name="FUNC", default='')
    exonic_func = table.Column(verbose_name="EXONIC_FUNC", default='')
    freq = table.Column(verbose_name="FREQ", default='')
    af = table.Column(verbose_name="AF", default='')
    precision = table.Column(verbose_name="PRECISION", default='', visible=False)

    class Meta:
        model = Results


class DellyTable(table.Table):
    result_id = table.Column(visible=False)
    run = table.Column(visible=False)
    sample = table.Column(verbose_name="SAMPLE", default='', visible=False)
    caller = table.Column(verbose_name="CALLER", default='', visible=False)
    precision = table.Column(verbose_name="PRECISION", default='')
    chrom = table.Column(verbose_name="CHROM", default='')
    pos = table.Column(verbose_name="POS", default='')
    ref = table.Column(verbose_name="REF", default='', visible=False)
    alt = table.Column(verbose_name="ALT", default='', visible=False)
    chr2 = table.Column(verbose_name="CHR2", default='')
    end = table.Column(verbose_name="END", default='')
    region = table.Column(verbose_name="REGION", default='')
    sv_type = table.Column(verbose_name="TYPE", default='')
    size = table.Column(verbose_name="SIZE", default='')
    gt = table.Column(verbose_name="GT", default='')
    depth = table.Column(verbose_name="DEPTH", default='')
    alleles = table.Column(verbose_name="ALLELES", default='')
    ab = table.Column(verbose_name="AB", default='')
    gene = table.Column(verbose_name="GENE", default='')
    func = table.Column(verbose_name="FUNC", default='')
    exonic_func = table.Column(verbose_name="EXONIC_FUNC", default='')
    freq = table.Column(verbose_name="FREQ", default='')
    af = table.Column(verbose_name="AF", default='')

    class Meta:
        model = Results


class Runs(models.Model):
    run_id = models.AutoField(primary_key=True, unique=True)
    run = models.CharField(max_length=100, default='Unknown')

    class Meta:
        app_label = 'aml'
        db_table = 'Runs'


class Samples(models.Model):
    sample_id = models.AutoField(primary_key=True, unique=True)
    sample = models.CharField(max_length=100, default='Unknown')
    run = models.CharField(max_length=100, default='Unknown')

    class Meta:
        app_label = 'aml'
        db_table = 'Samples'
