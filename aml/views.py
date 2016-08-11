from django.shortcuts import render
from aml.models import Runs, Results, Samples, PindelTable, DellyTable
from django.http import HttpResponse
from django_tables2 import RequestConfig
import pandas as pd
from collections import defaultdict


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
    colors = ['Blue', 'Red', 'Lime', 'Orange', 'DarkOrchid', 'Black', 'DarkCyan', 'DeepPink', 'SeaGreen', 'Yellow']
    color_no = 0
    color_dict = {}
    variant = None
    delly = DellyTable(Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Delly'))
    row = Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Delly')
    if request.method == 'POST':
        variant = request.POST.get('var')
    for item in row:
        color_dict[item.result_id] = colors[color_no]
        color_no += 1
    RequestConfig(request).configure(delly)
    return render(request, 'aml/delly.html', {'sample': sample, 'delly': delly, 'run': run, 'variant': variant,
                                              'color_dict': color_dict})


def get_interop_for_run(request, run):
    with open('/home/shjn/PycharmProjects/mypipeline/aml/static/aml/%s/%s_InterOp_Results.pdf' % (
            run, run), 'r') as pdf:
        response = HttpResponse(pdf.read(), content_type='application/pdf')
        response['Content-disposition'] = 'filename=%s_InterOp_Results.pdf' % run
        return response


def get_sample_quality(request, sample, run):
    with open('/home/shjn/PycharmProjects/mypipeline/aml/static/aml/%s/%s/%s_sample_quality.pdf' % (
            run, sample, sample), 'r') as pdf:
        response = HttpResponse(pdf.read(), content_type='application/pdf')
        response['Content-disposition'] = 'filename=%s_sample_quality.pdf' % sample
        return response


def get_flt3_only(request, sample, run):
    row = Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains='FLT3')
    d = defaultdict(list)
    flt3_sizes = []
    flt3_abs = []
    flt3_counts = []
    flt3_singles = []
    for item in row:
        if item.size in flt3_sizes:
            list_index = flt3_sizes.index(item.size)
            flt3_abs[list_index] = flt3_abs[list_index] + item.ab
            flt3_counts[list_index] += 1
            flt3_singles.append(item.size)
        else:
            flt3_sizes.append(item.size)
            flt3_abs.append(item.ab)
            flt3_counts.append(1)
    flt3_set = set(flt3_singles)
    for size in flt3_set:
        row = Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Pindel',
                                     gene__icontains='FLT3', size__icontains=size)
        for item in row:
            d[item.size].append(item.pos)
    flt3_combined = zip(flt3_sizes, flt3_abs, flt3_counts)
    var = bool(d)
    return render(request, 'aml/flt3.html', {'list': flt3_combined, 'sample': sample, 'run': run, 'd': dict(d),
                                             'var': var})


def get_fastqc(request, sample, run):
    fastqc_list = []
    stats_df = pd.read_table('aml/static/aml/%s/%s/fastqc_R1/fastqc_data.txt' % (run, sample), header=None,
                             names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
    properties = stats_df['Property'].tolist()
    values = stats_df['Value'].tolist()
    stats_dict = dict(zip(properties, values))
    stats_trim_df = pd.read_table('aml/static/aml/%s/%s/fastqc_R1_trim/fastqc_data.txt' % (run, sample), header=None,
                                  names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
    trim_properties = stats_trim_df['Property'].tolist()
    trim_values = stats_trim_df['Value'].tolist()
    stats_trim_dict = dict(zip(trim_properties, trim_values))

    fastqc_list.append(stats_dict["Filename"])
    fastqc_list.append(stats_dict["File type"])
    fastqc_list.append(stats_dict["Encoding"])
    fastqc_list.append(stats_dict["Total Sequences"])
    fastqc_list.append(stats_trim_dict["Total Sequences"])
    fastqc_list.append(stats_dict["Sequences flagged as poor quality"])
    fastqc_list.append(stats_trim_dict["Sequences flagged as poor quality"])
    fastqc_list.append(stats_dict["Sequence length"])
    fastqc_list.append(stats_trim_dict["Sequence length"])
    fastqc_list.append(stats_dict["%GC"])
    fastqc_list.append(stats_trim_dict["%GC"])
    fastqc_list.append(format((float(stats_trim_dict["Total Sequences"])/float(stats_dict["Total Sequences"])) * 100, '.1f'))

    r1 = get_summary('fastqc_R1', run, sample)
    r1_trim = get_summary('fastqc_R1_trim', run, sample)
    r2 = get_summary('fastqc_R2', run, sample)
    r2_trim = get_summary('fastqc_R2_trim', run, sample)

    return render(request, 'aml/fastqc.html', {'stats': fastqc_list, 'sample': sample, 'run': run, 'r1': r1, 'r2': r2,
                                               'r1_trim': r1_trim, 'r2_trim': r2_trim})


def get_summary(folder, run, sample):
    scores_df = pd.read_table('aml/static/aml/%s/%s/%s/summary.txt' % (run, sample, folder), header=None,
                              names=['Score', 'Parameter'], usecols=[0, 1])
    scores = scores_df['Score'].tolist()
    params = scores_df['Parameter'].tolist()
    scores_dict = dict(zip(params, scores))

    return scores_dict


def get_bamstats(request, run, sample):
    stats_file = '/home/shjn/PycharmProjects/mypipeline/aml/static/aml/%s/%s/bamstats/%s.bam.stats' \
                 % (run, sample, sample)
    bamstats_list = get_bamstat_scores(stats_file)
    return render(request, 'aml/bamstats.html', {'bamstats_list': bamstats_list, 'run': run, 'sample': sample})


def get_bamstat_scores(infile):
    bamstats_list = []
    scores_df = pd.read_table(infile, header=None, skiprows=7, nrows=30, usecols=[1, 2], names=['Parameter', 'Score'])
    scores = scores_df['Score'].tolist()
    params = scores_df['Parameter'].tolist()
    scores_dict = dict(zip(params, scores))

    filtered_reads = scores_dict["filtered sequences:"]
    total_reads = scores_dict["raw total sequences:"]
    remaining_reads = total_reads - filtered_reads
    mapped_reads = scores_dict["reads mapped:"]
    dup_reads = scores_dict["reads duplicated:"]
    mq0_reads = scores_dict["reads MQ0:"]
    mapped_bases = scores_dict["bases mapped (cigar):"]
    total_bases = scores_dict["total length:"]
    percent_filtered = format((filtered_reads * 100) / total_reads, '.1f')
    percent_mapped = format((mapped_reads * 100) / remaining_reads, '.1f')
    percent_dup = format((dup_reads * 100) / remaining_reads, '.1f')
    percent_mq0 = format((mq0_reads * 100) / mapped_reads, '.1f')
    bpercent_mapped = format(mapped_bases * 100 / total_bases, '.1f')
    error_rate = format(scores_dict["error rate:"] * 100, '.2f')

    bamstats_list.append("{:,}".format(int(scores_dict["reads paired:"])))
    bamstats_list.append("{:,}".format(int(filtered_reads)))
    bamstats_list.append(percent_filtered)
    bamstats_list.append("{:,}".format(int(scores_dict["non-primary alignments:"])))
    bamstats_list.append("{:,}".format(int(dup_reads)))
    bamstats_list.append(percent_dup)
    bamstats_list.append("{:,}".format(int(mapped_reads)))
    bamstats_list.append(percent_mapped)
    bamstats_list.append("{:,}".format(int(mq0_reads)))
    bamstats_list.append(percent_mq0)
    bamstats_list.append("{:,}".format(int(scores_dict["average length:"])))
    bamstats_list.append("{:,}".format(int(total_bases)))
    bamstats_list.append(bpercent_mapped)
    bamstats_list.append("{:,}".format(int(mapped_bases)))
    bamstats_list.append(error_rate)

    return bamstats_list
