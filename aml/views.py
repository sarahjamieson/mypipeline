from django.shortcuts import render
from aml.models import Runs, Results, Samples, PindelTable, DellyTable, Variants, VariantTable
from django.http import HttpResponse
from django_tables2 import RequestConfig
import pandas as pd
from django.template.defaulttags import register
import difflib
from collections import defaultdict
import os


def index(request):
    run = Runs.objects.all()
    return render(request, 'aml/index.html', {'runs': run})


def get_samples_for_run(request, run):
    samples = Samples.objects.filter(run__icontains=run).order_by('sample')
    return render(request, 'aml/samples.html', {'samples': samples, 'run': run})


def get_results_for_sample(request, sample, run):
    gene = None
    ab_min = None
    if 'gene_filter' in request.GET:
        gene = request.GET['gene_filter'].upper()
    else:
        pass
    if 'ab_min' in request.GET:
        ab_min = request.GET['ab_min']
    else:
        pass
    if gene and ab_min:
        pindel = PindelTable(Results.objects.filter(
            sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains=gene, ab__gte=ab_min
        ))
    elif gene and not ab_min:
        pindel = PindelTable(Results.objects.filter(
            sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains=gene
        ))
    elif not gene and ab_min:
        pindel = PindelTable(Results.objects.filter(
            sample__icontains=sample, run__icontains=run, caller='Pindel', ab__gte=ab_min
        ))
    else:
        pindel = PindelTable(Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Pindel'))

    RequestConfig(request).configure(pindel)
    return render(request, 'aml/results.html', {'pindel': pindel, 'sample': sample, 'run': run, 'gene': gene})


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
        if color_no > 9:
            color_no = 0
        else:
            color_no = color_no
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
    pindel = PindelTable(
        Results.objects.filter(
            sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains='FLT3'
        ).order_by('size', 'pos')
    )
    runs_dict = dict()
    ab_dict = dict()
    samples_dict = dict()
    results_duplicate_sizes = []
    counts = []
    d = defaultdict(list)
    seq_length = 0
    for item in row:
        total_ab = 0
        size_count = Results.objects.filter(
            sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains='FLT3', size=item.size
        )
        if len(size_count) > 1:
            d[item.size].append(item.pos)
            results_duplicate_sizes.append(item.result_id)
            counts.append(len(size_count))
            filename = "%s_%s.txt" % (sample, item.size)
            with open("%s" % filename, "w+") as seq_file:
                for result in size_count:
                    total_ab += result.ab
                    seq_file.write(">chr13:%s\n" % result.pos)
                    seq_file.write("%s\n" % result.alt)
                    seq_length = len(result.alt)
            seq_file.close()
            os.system("/home/shjn/meme/bin/meme %s -dna -mod oops -w %s" % (filename, seq_length))
            os.system("cp /home/shjn/PycharmProjects/mypipeline/aml/meme_out/logo1.png "
                      "/media/sf_S_DRIVE/MiSeq_data/Nextera_Rapid_Capture/Sarah_STP_Project_AML/%s/%s/%s_%s.meme.png"
                      % (run, sample, sample, item.size))
            os.system("rm %s" % filename)
        ab_dict[item.size] = total_ab
        wide_search = Results.objects.filter(chrom=item.chrom)
        run_list = []
        sample_list = []
        for result in wide_search:
            seq = difflib.SequenceMatcher(a=item.alt.lower(), b=result.alt.lower())
            if seq.ratio() > 0.90 and (item.pos - 60) < result.pos < (item.pos + 60) and item.size == result.size:
                run_list.append(result.run)
                sample_list.append(result.sample)
            else:
                pass
        no_of_unique_runs = len(set(run_list))
        no_of_unique_samples = len(set(sample_list))
        runs_dict[item.result_id] = no_of_unique_runs
        samples_dict[item.result_id] = no_of_unique_samples

    counts_list = [results_duplicate_sizes[i] for i in range(0, len(results_duplicate_sizes), 2)]
    no_counts_list = [results_duplicate_sizes[i] for i in range(1, len(results_duplicate_sizes), 2)]
    counts_altered = [counts[i] for i in range(0, len(counts), 2)]
    counts_dict = dict()
    counter = 0
    for item in no_counts_list:
        counts_dict[item] = counts_altered[counter]
        counter += 1
    return render(request, 'aml/flt3.html', {'sample': sample, 'run': run, 'pindel': pindel, 'runs': runs_dict,
                                             'samples': samples_dict, 'ab': ab_dict, 'counts': counts_list,
                                             'no_counts': no_counts_list, 'counts_d': counts_dict, 'd': dict(d)})


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


def view_variants(request):
    variants = Variants.objects.all()
    return render(request, 'aml/var_table.html', {'var_table': variants})


def view_polyphen(request, result_id):
    variants = Variants.objects.filter(result_id=result_id)



@register.filter
def lookup(dictionary, key):
    return dictionary.get(key)
