from django.shortcuts import render
from aml.models import Runs, Results, ResultTable, Samples


def index(request):
    run = Runs.objects.all()
    return render(request, 'aml/index.html', {'run': run})


def get_samples_for_run(request, run):
    samples = Samples.objects.filter(run__icontains=run)
    return render(request, 'aml/samples.html', {'samples': samples, 'run': run})


def get_results_for_sample(request, sample):
    results = ResultTable(Results.objects.filter(sample__icontains=sample))
    return render(request, 'aml/results.html', {'results': results})
