import pandas as pd
import os
script_dir = os.path.dirname(os.path.abspath(__file__))


def get_fastqc_results(run, sample):
    """Reads results before and after trimming from fastqc_data.txt file and extracts values into lists. Also gets
    results from summary.txt file for each read before and after trimming (r1, r2, r1_trim, r2_trim).
    """

    # (1) Data into dataframe, pull out property and value data and create dictionary.
    stats_df = pd.read_table('%s/static/aml/%s/%s/fastqc_R1/fastqc_data.txt' % (script_dir, run, sample), header=None,
                             names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
    properties = stats_df['Property'].tolist()
    values = stats_df['Value'].tolist()
    stats_dict = dict(zip(properties, values))
    stats_trim_df = pd.read_table('%s/static/aml/%s/%s/fastqc_R1_trim/fastqc_data.txt' % (script_dir, run, sample),
                                  header=None, names=['Property', 'Value'], usecols=[0, 1], skiprows=3, nrows=7)
    trim_properties = stats_trim_df['Property'].tolist()
    trim_values = stats_trim_df['Value'].tolist()
    stats_trim_dict = dict(zip(trim_properties, trim_values))

    # (2) Get summary data into dictionaries.
    r1 = get_summary('fastqc_R1', run, sample)
    r1_trim = get_summary('fastqc_R1_trim', run, sample)
    r2 = get_summary('fastqc_R2', run, sample)
    r2_trim = get_summary('fastqc_R2_trim', run, sample)

    return r1, r2, r1_trim, r2_trim, stats_dict, stats_trim_dict


def get_summary(folder, run, sample):
    """Get data from .txt file and put into dictionary."""
    scores_df = pd.read_table('%s/static/aml/%s/%s/%s/summary.txt' % (script_dir, run, sample, folder), header=None,
                              names=['Score', 'Parameter'], usecols=[0, 1])
    scores = scores_df['Score'].tolist()
    params = scores_df['Parameter'].tolist()
    scores_dict = dict(zip(params, scores))

    return scores_dict


def get_bamstat_scores(infile):
    """Scores obtained from .bam.stats file and put into dictionary.

    Scores are then extracted and put into a list as many parameters require modifications which are easier to do in
    Python that in the Django template.

    """
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
