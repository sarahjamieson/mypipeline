from django.shortcuts import render
from aml.models import Runs, Results, Samples, PindelTable, DellyTable
from django.http import HttpResponse
from django_tables2 import RequestConfig
import os


def index(request):
    run = Runs.objects.all()
    return render(request, 'aml/index.html', {'runs': run})


def get_samples_for_run(request, run):
    samples = Samples.objects.filter(run__icontains=run).order_by('sample')
    return render(request, 'aml/samples.html', {'samples': samples, 'run': run})


def get_results_for_sample(request, sample, run):
    pindel = PindelTable(Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Pindel'))
    RequestConfig(request).configure(pindel)
    return render(request, 'aml/results.html', {'pindel': pindel, 'sample': sample, 'run': run})


def get_delly_for_sample(request, sample, run):
    delly = DellyTable(Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Delly'))
    RequestConfig(request).configure(delly)
    var = False
    if os.path.isfile("/home/cuser/PycharmProjects/django_apps/mypipeline/aml/static/rcircos/%s/%s.png" % (run, sample)):
        var = True
    return render(request, 'aml/delly.html', {'sample': sample, 'delly': delly, 'run': run, 'var': var})


def get_interop_for_run(request, run):
    with open('/home/cuser/PycharmProjects/django_apps/mypipeline/aml/static/aml/%s/%s_InterOp_Results.pdf' % (
            run, run), 'r') as pdf:
        response = HttpResponse(pdf.read(), content_type='application/pdf')
        response['Content-disposition'] = 'filename=%s_InterOp_Results.pdf' % run
        return response


def get_sample_quality(request, sample, run):
    with open('/home/cuser/PycharmProjects/django_apps/mypipeline/aml/static/aml/%s/%s_sample_quality.pdf' % (
            run, sample), 'r') as pdf:
        response = HttpResponse(pdf.read(), content_type='application/pdf')
        response['Content-disposition'] = 'filename=%s_sample_quality.pdf' % sample
        return response


def get_circos(request, sample, run):
    return render(request, 'aml/circos.html', {'sample': sample, 'run': run})
