from django.shortcuts import render
from aml.models import Runs, Results, Samples, PindelTable, DellyTable, Variants, FragmentAnalysis
from django.http import HttpResponse
from django_tables2 import RequestConfig
from django.template.defaulttags import register
import difflib
import sqlite3 as sql
import os
import polyphen
import time
from chartit import DataPool, Chart
from get_quality_results import get_fastqc_results, get_bamstat_scores
import pandas as pd
from generate_R_plots import generate_ngs_vs_frag, generate_ngs_vs_frag_exclude_outliers

script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))
config_df = pd.read_table('%s/config.txt' % parent_dir, header=None, names=['Setting', 'Value'], sep='=')
settings = config_df['Setting'].tolist()
values = config_df['Value'].tolist()
config_dict = dict(zip(settings, values))


def index(request):
    """Returns all entries in Runs model to index.html page.
    Note: could remove this model and pull out only unique runs from Results model.
    """
    runs = Runs.objects.all()
    return render(request, 'aml/index.html', {'runs': runs})


def get_samples_for_run(request, run):
    """Takes a run and returns all collated entires from the Samples model to samples.html page.
    Note: as above, could remove this model and pull out only unique samples from Results model.
    """
    samples = Samples.objects.filter(run__icontains=run).order_by('sample')
    return render(request, 'aml/samples.html', {'samples': samples, 'run': run})


def get_results_for_sample(request, sample, run):
    """Takes a run and sample and returns matching Pindel entries from Results model to results.html page. Users can
    filter data by gene name or minimum allelic ratio (ab) which is handled by request.GET.

    """
    gene = None
    ab_min = None

    # Need to treat parameters separately in case the user only wants to filter by one of them.
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
    RequestConfig(request).configure(pindel)  # allows sorting of django_tables2 columns
    return render(request, 'aml/results.html', {'pindel': pindel, 'sample': sample, 'run': run, 'gene': gene})


def get_delly_for_sample(request, sample, run):
    """Takes a run and sample and returns matching Delly entries from Results model to delly.html page. Each variant is
    matched to a colour corresponding to the colour in the sample's circos plot.

    Note: need to update to colour match only PRECISE results -- remove IMPRECISE results??
    """
    color_no = 0
    color_dict = {}
    colors = ['Blue', 'Red', 'Lime', 'Orange', 'DarkOrchid', 'Black', 'DarkCyan', 'DeepPink', 'SeaGreen', 'Yellow']
    delly_results = Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Delly')
    delly_table = DellyTable(Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Delly'))
    for result in delly_results:
        if color_no > 9:
            color_no = 0
        else:
            color_no = color_no
        color_dict[result.result_id] = colors[color_no]
        color_no += 1
    RequestConfig(request).configure(delly_table)  # allows sorting of django_tables2 columns
    return render(request, 'aml/delly.html', {'sample': sample, 'delly': delly_table, 'run': run,
                                              'color_dict': color_dict})


def get_interop_for_run(request, run):
    """Takes a run and downloads the InterOp PDF file saved in the static folder for that run.

    Note: make path to static folder general.

    """
    with open('%s/static/aml/%s/%s_InterOp_Results.pdf' % (
            script_dir, run, run), 'r') as pdf:
        response = HttpResponse(pdf.read(), content_type='application/pdf')
        response['Content-disposition'] = 'filename=%s_InterOp_Results.pdf' % run
        return response


def get_sample_quality(request, sample, run):
    """Takes a sample and run and downloads the Quality PDF file saved in the static folder for that sample.

    Note: "request" is required as argument for browser connection even though not used.

    """
    with open('%s/static/aml/%s/%s/%s_sample_quality.pdf' % (script_dir, run, sample, sample), 'r') as pdf:
        response = HttpResponse(pdf.read(), content_type='application/pdf')
        # browser will treat file as attachment
        response['Content-disposition'] = 'filename=%s_sample_quality.pdf' % sample
        return response


def get_flt3_only(request, sample, run):
    """Takes a sample and run and provides a more detailed analysis of the FLT3 gene results for that sample.

    (1) Pulls FLT3 entries from Results model for sample and searches for each result across all samples and runs to
        calculate the number of occurrences. A results is considered the same if it has a sequence match of >90%
        (using difflib.SequenceMatcher), the ITD size is identical, and the ITDs are within 60bp of each other.

    (2) Part of the analysis combines ITDs of the same size in the sample and generates a combined allelic ratio (ab).

    (3) All FLT3 data for sample is added to a dictionary and then this dictionary is added to a list to be displayed
        as a table in the template.

    (3) For results where there are other results of the same size, these are written to a .txt file and used as input
        to "meme" which identifies sequence motifs. This image is then displayed in the html page.

        Meme:
            "meme <input.txt> -dna -mod oops -w <length>"

            -dna : knows to recognise DNA sequences
            -mod : oops = expect one occurrence of the motif in the aligned sequences
            -w : search for motif with length <length>

        Duplicate results are also added to a dictionary so the size and positions can be displayed in the template.

    """
    flt3_results = []
    duplicate_dict = {}

    # Find all FLT3 results for sample in run.
    flt3_results_for_sample = Results.objects.filter(
        sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains='FLT3'
    ).order_by('size', 'pos')

    # For each result...
    for f in flt3_results_for_sample:
        run_list = []
        sample_list = []
        ab_total = 0
        seq_length = 0

        # ...search for similar results across all samples and runs.
        all_results_for_chrom = Results.objects.filter(chrom=f.chrom)
        for c in all_results_for_chrom:
            seq = difflib.SequenceMatcher(a=f.alt.lower(), b=c.alt.lower())
            if seq.ratio() > 0.90 and (f.pos - 60) < c.pos < (f.pos + 60) and f.size == c.size:
                run_list.append(c.run)
                sample_list.append(c.sample)
            else:
                pass
        no_of_unique_runs = len(set(run_list))
        no_of_unique_samples = len(set(sample_list))

        # ...find other results with the same size in this sample/run and combine to calculate a total AB.
        flt3_results_with_same_size = Results.objects.filter(
            sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains='FLT3', size=f.size
        )
        for s in flt3_results_with_same_size:
            ab_total += s.ab

        # ...if there are multiple results for the size, add to dictionary and run meme.
        if len(flt3_results_with_same_size) > 1:
            if f.size in duplicate_dict:
                duplicate_dict[f.size].append(f.pos)
            else:
                duplicate_dict[f.size] = [f.pos]

            filename = "%s_%s.txt" % (sample, f.size)
            with open("%s" % filename, "w+") as seq_file:
                for s in flt3_results_with_same_size:
                    seq_file.write(">chr13:%s\n%s\n" % (s.pos, s.alt))
                    seq_length = len(s.alt)

            seq_file.close()
            os.system("%s %s -dna -mod oops -w %s" % (config_dict['meme'], filename, seq_length))
            os.system("cp %s/meme_out/logo1.png %s/%s/%s/%s_%s.meme.png"
                      % (parent_dir, config_dict['static_link_folder'], run, sample, sample, f.size))
            os.system("rm %s" % filename)

        # ...add all data to a dictionary which is then added to a list of dictionaries.
        results = {'id': f.result_id,
                   'pos': f.pos,
                   'ref': f.ref,
                   'alt': f.alt,
                   'size': f.size,
                   'ab': ab_total,
                   'runs': no_of_unique_runs,
                   'samples': no_of_unique_samples}
        flt3_results.append(results)

    return render(request, 'aml/flt3.html', {'run': run, 'sample': sample, 'flt3': flt3_results, 'd': duplicate_dict})


def get_fastqc(request, sample, run):
    """Takes a sample and run, obtains quality data from fastqc_data.txt and summary.txt files generated by FastQC and
    returns to fastqc.html page.

        fastqc_list: list of results for each parameter e.g. file type, total sequences.
        r1/r2: dictionary of scores for each measurement before trimming e.g. per base sequence quality.
        r1_trim/r2_trim: dictionary of scores for each measurement after trimming.

    """
    r1, r2, r1_trim, r2_trim, stats_dict, stats_trim_dict = get_fastqc_results(run, sample)

    return render(request, 'aml/fastqc.html', {'sample': sample, 'run': run, 'r1': r1, 'r2': r2, 'r1_trim': r1_trim,
                                               'r2_trim': r2_trim, 'stats_dict': stats_dict,
                                               'trim_dict': stats_trim_dict})


def get_bamstats(request, run, sample):
    """Takes a run and sample and returns BAMStats from .bam.stats file to bamstats.html page."""
    stats_file = '%s/static/aml/%s/%s/bamstats/%s.bam.stats' % (script_dir, run, sample, sample)
    bamstats_list = get_bamstat_scores(stats_file)
    return render(request, 'aml/bamstats.html', {'bamstats_list': bamstats_list, 'run': run, 'sample': sample})


def view_variants(request):
    """View for dummy page to represent how PolyPhen-2 could be integrated. Variants model contains two BRCA2 SNVs as
    example.

    This returns the variants to the var_table.html page.
    """
    variants = Variants.objects.all()
    return render(request, 'aml/var_table.html', {'var_table': variants})


def view_polyphen(request, result_id):
    """View for dummy page to represent how PolyPhen-2 could be integrated. For the selected variant, a connection to
    PolyPhen API is made using curl and the session ID is stored. Then the results are obtained using the session ID
    and urllib2 and parsed into two dictionaries (HumDiv results and HumVar results).

    Note: delay of 30 seconds required for session ID to be recognised through urllib2 connection.

    """
    variants = Variants.objects.filter(result_id=result_id)
    div_dict = {}
    var_dict = {}
    for item in variants:
        filename = polyphen.create_polyphen_input(item.chrom, item.pos, item.ref, item.alt)
        humdiv, humvar = polyphen.run_polyphen(filename)
        time.sleep(30)
        div_dict = polyphen.get_polyphen_results(humdiv)
        var_dict = polyphen.get_polyphen_results(humvar)
    polyphen.get_polyphen_colormap(div_dict.get('prob_score'), 'humdiv')
    polyphen.get_polyphen_colormap(var_dict.get('prob_score'), 'humvar')

    return render(request, 'aml/polyphen.html', {'div_dict': div_dict, 'var_dict': var_dict})


def view_flt3_combined(request):
    """Returns a table and graphs with combined FLT3 results from all samples.

    Three parts to the template:
        (1) Form to input Fragment Analysis data (users selects run, sample, ITD size and allelic ratio).
        (2) Table of results comparing Fragment Analysis to NGS data.
        (3) Graphs showing the results.

    Note: currently multiple runs with same samples in this version, so results are filtered for the two main runs to
    avoid lots of duplicates (16053 and 160805).
    """
    show_plot = False
    show_matches = False
    samples = []
    runs = []
    unique_samples = []
    unique_runs = []
    matching_samples = []

    # Generate unique sample list and unique run list for drop down menus in template.
    results = Results.objects.all().order_by('sample')
    frag_results = FragmentAnalysis.objects.all()
    for result in results:
        samples.append(result.sample)
        runs.append(result.run)
    for sample in samples:
        if sample not in unique_samples:
            unique_samples.append(sample)
    for run in runs:
        if run not in unique_runs:
            unique_runs.append(run)

    # If the form is filled in, results are added to the database and FragmentAnalysis django model.
    if 'option' in request.GET and 'itd' in request.GET and 'ab' in request.GET and 'run' in request.GET:
        sample = request.GET['option']
        itd = request.GET['itd']
        ab = request.GET['ab']
        run = request.GET['run']
        con = sql.connect('%s/db.sqlite3' % parent_dir)
        curs = con.cursor()
        curs.execute("CREATE TABLE IF NOT EXISTS Frag(result_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "
                     "sample TEXT, run TEXT, itd INTEGER, ab REAL)")
        curs.execute("INSERT INTO Frag (sample, run, itd, ab) VALUES (?, ?, ?, ?)", (sample, run, itd, ab))
        con.commit()
    else:
        pass

    # If user selects button in template named "graph", get all fragment analysis results and if a fragment analysis
    # itd was also detected by NGS then results are added to a list. Then an R script is written to produce a ggplot of
    # fragment analysis ratios versus NGS ratios with size as an additional variable. The resulting graph is saved in
    # the static folder.
    if 'graph' in request.GET:
        show_plot = True
        final = []
        frag = []
        ngs = []
        itd = []
        frags = FragmentAnalysis.objects.all()
        for item in frags:
            ngs_temp = []
            sample_results = Results.objects.filter(sample=item.sample, run=item.run, size=item.itd)
            if sample_results:
                for result in sample_results:
                    ar = float(format(float(result.alleles.split(',')[1]) / float(result.alleles.split(',')[0]), '.3f'))
                    ngs_temp.append(ar)
                frag.append(item.ab)
                ngs.append(format(sum(ngs_temp), '.3f'))
                itd.append(item.itd)
                final.append([item.itd, item.ab, format(sum(ngs_temp), '.3f')])
            else:
                pass
        generate_ngs_vs_frag(frag, ngs, itd)  # ggplot function
        generate_ngs_vs_frag_exclude_outliers(frag, ngs, itd)  # ggplot function
    else:
        pass

    # Comparison table to look like this:
    #           |      |Fragment Analysis|       NGS       |
    #           |      |-----------------|-----------------|
    #           |      |  Size  |   AR   |  Size  |   AR   |
    #           |------|--------|--------|--------|--------|
    #           |Sample|        |        |        |        |
    # For every sample in database, get corresponding fragment analysis results.
    #  - If results are found, combine all allelic ratios in a list and all sizes in a list.
    #  - Then add these two lists to a nested dictionary with the sample name as the key.
    #  - Now search for FLT3 Pindel NGS results for sample.
    #  - If results are found, combine all allelic ratios in a list and all sizes in a list. For results with the same
    #    size, add together the allelic ratios.
    #  - Add the two lists to a nested dictionary with the sample name as the key.
    # These two dictionaries will then be called in the template to generate the table.
    if 'toggle' in request.GET:
        show_plot = True
        show_matches = True
        ngs_dict, frag_dict, matching_samples = get_matching_flt3(unique_samples)
    else:
        ngs_dict, frag_dict = get_all_flt3(unique_samples)

    # Get list of unique samples with FLT3 results to use as keys for dictionaries to generate table.
    flt3_samples = []
    flt3_results = Results.objects.filter(gene__icontains='FLT3', caller='Pindel',
                                          run='160628_merged') | Results.objects.filter(gene__icontains='FLT3',
                                                                                        caller='Pindel', run='160805')
    for item in flt3_results:
        flt3_samples.append(item.sample)
    for item in frag_results:
        if item.sample not in flt3_samples:
            flt3_samples.append(item.sample)
    flt3_samples = list(set(flt3_samples))
    flt3_samples.sort()
    if 'toggle' in request.GET:
        flt3_samples = matching_samples
    else:
        flt3_samples = flt3_samples

    # Django chartit example
    # 1) Create DataPool with data for chart
    frag_data = DataPool(
        series=
        [{'options': {'source': FragmentAnalysis.objects.all()}, 'terms': ['itd', 'ab']}])
    # 2) Create chart object
    cht = Chart(
        datasource=frag_data,
        series_options=[
            {'options': {
                'type': 'scatter'},
                'terms': {
                    'itd': [
                        'ab'
                    ]
                }}
        ],
        chart_options={
            'title': {
                'text': 'Fragment Analysis size vs allelic ratio'
            },
            'xAxis': {
                'title': {
                    'text': 'ITD size'
                }
            }
        }
    )

    return render(request, 'aml/combined.html', {'samples': unique_samples, 'runs': unique_runs, 'show': show_plot,
                                                 'frag_d': frag_dict, 'ngs_d': ngs_dict, 'flt3_samples': flt3_samples,
                                                 'fragchart': cht, 'show_matches': show_matches})


@register.filter
def lookup(dictionary, key):
    """Django custom filter. In template, do:
            {{ dictionary|lookup:key }}
    """
    return dictionary.get(key)


@register.simple_tag
def nested_lookup(dictionary, key1, key2):
    """Django custom filter. In template, do:
            {{ nested_lookup dictionary key1 key2 }}
    """
    return dictionary[key1][key2]


def get_all_flt3(unique_samples):
    frag_dict = {}
    ngs_dict = {}
    for sample in unique_samples:
        abs = []
        sizes = []
        abs_ngs = []
        sizes_ngs = []
        frags = FragmentAnalysis.objects.filter(sample=sample).order_by('sample', 'itd')
        if frags:
            for item in frags:
                abs.append(item.ab)
                sizes.append(item.itd)
            frag_dict[sample] = {
                'size': sizes,
                'ab': abs
            }
        else:
            pass
        results = Results.objects.filter(
            sample=sample, gene__icontains='FLT3', caller='Pindel').order_by('sample', 'size')
        if results:
            for item in results:
                if item.run == '160628_merged' or item.run == '160805':  # ideally remove this if statement
                    ar = float(item.alleles.split(',')[1]) / float(item.alleles.split(',')[0])
                    ar = float(format(ar, '.3f'))
                    if item.size in sizes_ngs:
                        size_index = sizes_ngs.index(item.size)
                        abs_ngs[size_index] += ar
                    else:
                        sizes_ngs.append(item.size)
                        abs_ngs.append(ar)
                else:
                    pass
            if sizes_ngs:
                ngs_dict[sample] = {
                    'size': sizes_ngs,
                    'ab': abs_ngs
                }

    return ngs_dict, frag_dict


def get_matching_flt3(unique_samples):
    ngs_sample_sizes = []
    frag_sample_sizes = []
    frag_dict = {}
    ngs_dict = {}
    samples_with_matches = []

    for result in Results.objects.filter(caller='Pindel', gene__icontains='FLT3',
                                         run='160628_merged') | Results.objects.filter(caller='Pindel',
                                                                                       gene__icontains='FLT3',
                                                                                       run='160805'):
        ngs_sample_sizes.append("%s_%s" % (result.sample, result.size))
    for frag in FragmentAnalysis.objects.all():
        frag_sample_sizes.append("%s_%s" % (frag.sample, frag.itd))

    for sample in unique_samples:
        abs = []
        sizes = []
        abs_ngs = []
        sizes_ngs = []
        frags = FragmentAnalysis.objects.filter(sample=sample).order_by('sample', 'itd')
        if frags:
            for item in frags:
                check = "%s_%s" % (sample, item.itd)
                if check in ngs_sample_sizes and check in frag_sample_sizes:
                    samples_with_matches.append(sample)
                    abs.append(item.ab)
                    sizes.append(item.itd)
            frag_dict[sample] = {
                'size': sizes,
                'ab': abs
            }
        else:
            pass
        results = Results.objects.filter(
            sample=sample, gene__icontains='FLT3', caller='Pindel').order_by('sample', 'size')
        if results:
            for item in results:
                check = "%s_%s" % (sample, item.size)
                if check in ngs_sample_sizes and check in frag_sample_sizes:
                    if item.run == '160628_merged' or item.run == '160805':  # ideally remove this if statement
                        ar = float(item.alleles.split(',')[1]) / float(item.alleles.split(',')[0])
                        ar = float(format(ar, '.3f'))
                        if item.size in sizes_ngs:
                            size_index = sizes_ngs.index(item.size)
                            abs_ngs[size_index] += ar
                        else:
                            sizes_ngs.append(item.size)
                            abs_ngs.append(ar)
                if sizes_ngs:
                    ngs_dict[sample] = {
                        'size': sizes_ngs,
                        'ab': abs_ngs
                    }

    return ngs_dict, frag_dict, list(set(samples_with_matches))
