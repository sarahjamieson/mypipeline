from django.shortcuts import render
from aml.models import Runs, Results, Samples, PindelTable, DellyTable
from django.http import HttpResponse
from django_tables2 import RequestConfig


def index(request):
    run = Runs.objects.all()
    return render(request, 'aml/index.html', {'runs': run})


def get_samples_for_run(request, run):
    samples = Samples.objects.filter(run__icontains=run).order_by('sample')
    return render(request, 'aml/samples.html', {'samples': samples, 'run': run})


def get_results_for_sample(request, sample):
    pindel = PindelTable(Results.objects.filter(sample__icontains=sample, caller='Pindel'))
    RequestConfig(request).configure(pindel)
    return render(request, 'aml/results.html', {'pindel': pindel, 'sample': sample})


def get_delly_for_sample(request, sample):
    delly = DellyTable(Results.objects.filter(sample__icontains=sample, caller='Delly'))
    RequestConfig(request).configure(delly)
    return render(request, 'aml/delly.html', {'sample': sample, 'delly': delly})


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
