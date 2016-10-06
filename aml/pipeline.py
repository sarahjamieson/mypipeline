import sqlite3 as sql
from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import os
import datetime
from ruffus import *
import glob
from create_interop_pdf import CreatePDF
import vcf
import pandas as pd
import re
from pandas import ExcelWriter
import parse_vcfs
import argparse
from argparse import RawTextHelpFormatter
from parse_sample_sheet import ParseSampleSheet
import shelve
from collections import defaultdict
import warnings
from sorting import sort_pipeline_data

warnings.simplefilter(action="ignore", category=FutureWarning)
curr_datetime = datetime.datetime.now().isoformat()
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.abspath(os.path.join(script_dir, os.pardir))

# ----------------------------------------------------------------------------------------------------------------------
# Command line parameters for pipeline script:
#   -s = Sample Sheet
#   -d = Directory with results (InterOp and Data folders)
#   -o = Name of output directory
# ----------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Runs pipeline for Illumina MiSeq Nextera AML data.",
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-s', '--sheet', action="store", dest='sample_sheet', help='An illumina sample sheet for this run.',
                    default='/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/SampleSheet.csv')
parser.add_argument('-d', '--result_dir', action='store', dest='result_dir',
                    help='Directory containing Illumina MiSeq InterOp and Data folders.',
                    default='/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/')
parser.add_argument('-o', '--output_dir', action='store', dest='output_dir',
                    help='Directory for output to be stored in',
                    default='/media/sf_sarah_share/MiSeq_Nextera_Results/')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------------------------------------
# Get values from config.txt file.
# ----------------------------------------------------------------------------------------------------------------------
InterOp = '%sInterOp/' % args.result_dir
config_df = pd.read_table('%s/config.txt' % parent_dir, header=None, names=['Setting', 'Value'], sep='=')
settings = config_df['Setting'].tolist()
values = config_df['Value'].tolist()
config_dict = dict(zip(settings, values))

# ----------------------------------------------------------------------------------------------------------------------
# 1) Copy fastqc files to current directory.
# 2) Parse sample sheet to get run and sample info.
# 3) Make new directory for results in user-inputted output directory.
# ----------------------------------------------------------------------------------------------------------------------
os.system("cp %s/Data/Intensities/BaseCalls/*.fastq.gz %s/" % (args.result_dir, script_dir))
parse_sheet = ParseSampleSheet(args.sample_sheet)
run_dict, sample_dict = parse_sheet.parse_sample_sheet()
worksheet = run_dict.get('worksheet')
os.system("mkdir -p %s%s" % (args.output_dir, worksheet))


# ----------------------------------------------------------------------------------------------------------------------
# Assess run quality based on InterOp data and output PDF file.
# ----------------------------------------------------------------------------------------------------------------------
def assess_quality():
    """Obtains quality statistics from MiSeq InterOp files and produces PDF summary report.
    Notes:
        Uses Python module Illuminate 0.6.2. to read InterOp file data (https://pypi.python.org/pypi/illuminate/).
        See http://rpackages.ianhowson.com/bioc/savR/ for details on column meanings in InterOp files.
    """
    tilemetrics = InteropTileMetrics('%sTileMetricsOut.bin' % InterOp)
    controlmetrics = InteropControlMetrics('%sControlMetricsOut.bin' % InterOp)
    errormetrics = InteropErrorMetrics('%sErrorMetricsOut.bin' % InterOp)
    extractionmetrics = InteropExtractionMetrics('%sExtractionMetricsOut.bin' % InterOp)
    indexmetrics = InteropIndexMetrics('%sIndexMetricsOut.bin' % InterOp)
    qualitymetrics = InteropQualityMetrics('%sQMetricsOut.bin' % InterOp)
    corintmetrics = InteropCorrectedIntensityMetrics('%sCorrectedIntMetricsOut.bin' % InterOp)

    pdf = CreatePDF(tilemetrics, controlmetrics, errormetrics, extractionmetrics, indexmetrics, qualitymetrics,
                    corintmetrics, worksheet, config_dict['pdflatex'])
    pdf.create_pdf()

    os.system('cp %s_InterOp_Results.pdf %s%s/' % (worksheet, args.output_dir, worksheet))
    filelist = glob.glob("*tex")
    for f in filelist:
        os.remove(f)

assess_quality()


# ----------------------------------------------------------------------------------------------------------------------
# Run pipeline:
#   1) Run FastQC on raw fastqc files.
#   2) Run Trimmomatic [.fastq.gz --> .qfilter.fastq.gz]
#   3) Run FastQC on trimmed fastqc files.
#   4) Run BWA to align reads to hg19 [.qfilter.fastq.gz --> .bwa.drm.bam].
#   5) Run SAMtools "Sort" on aligned data [.bwa.drm.bam --> .bwa.drm.sorted.bam].
#   6) Run SAMtools "Index" on sorted data to generate index file [.bwa.drm.sorted.bam --> .bwa.drm.sorted.bam.bai].
#   7) Run SAMtools "Stats" and "Plot-bamstats" & combine with FastQC results.
#   8) Create breakdancer config and run breakdancer.
#   9) Create pindel config, run pindel and convert pindel output to VCF file.
#  10) Run delly and convert BCF output to VCF file.
#  11) Merge all VCFs.
#  12) Annotate merged VCF using ANNOVAR.
#  13) Calculate frequency of each variant in data so far by compressing and combining all variants into database, and
#      add to annotated VCF file.
#  14) Convert annotated VCF file to an excel spreadsheet and add to db.sqlite3 database for use in web app.
# ----------------------------------------------------------------------------------------------------------------------

@collate("*.fastq.gz", formatter("([^/]+)R[12]_001.fastq.gz$"), "{path[0]}/{1[0]}.fastq.gz")
def run_fastqc(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (config_dict['fastqc'], fastq1))
    os.system("%s --extract %s" % (config_dict['fastqc'], fastq2))


@follows(run_fastqc)
@collate("*.fastq.gz", formatter("([^/]+)R[12]_001.fastq.gz$"), "{path[0]}/{1[0]}.qfilter.fastq.gz")
def quality_trim(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]

    os.system("java -Xmx4g -jar %s "
              "PE "  # Paired-end mode
              "-threads 2 "
              "-phred33 "  # Phred33 is used in Illumina versions >= 1.8.
              "%s %s "  # Input files (R1 & R2)
              "%s.qfilter.fastq.gz %s.unpaired.fastq.gz "  # Paired and unpaired output for R1
              "%s.qfilter.fastq.gz %s.unpaired.fastq.gz "  # Paired and unpaired output for R2
              "ILLUMINACLIP:%s:2:30:10 "  # Adaptor trimming [adapters : allow 2 mismatches : PE >= Q30 : SE >= Q10]
              "CROP:150 "  # Crop reads to max 150 from end
              "SLIDINGWINDOW:4:25 "  # Perform trim on every 4bp for minimum quality of 25.
              "MINLEN:50"  # All reads should be minimum length of 50bp.
              % (config_dict['Trimmomatic'], fastq1, fastq2, fastq1[:-9], fastq1, fastq2[:-9], fastq2,
                 config_dict['trim_adapters']))

    # glob module: finds path names matching specified pattern (https://docs.python.org/2/library/glob.html).
    # Delete unpaired data no longer required.
    filelist = glob.glob("*unpaired*")
    for f in filelist:
        os.remove(f)


@follows(quality_trim)
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12]_001.qfilter.fastq.gz$"), "{path[0]}/{1[0]}.fastq.gz")
def run_fastqc_trimmed(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (config_dict['fastqc'], fastq1))
    os.system("%s --extract %s" % (config_dict['fastqc'], fastq2))


@follows(run_fastqc_trimmed)
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12]_001.qfilter.fastq.gz$"), "{path[0]}/{1[0]}.bwa.drm.bam")
def align_bwa(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    sample_name = fastq1[:-20]
    rg_header = '@RG\tID:%s\tSM:%s\tCN:WMRGL\tDS:TruSight_Myeloid_Nextera_v1.1\tDT:%s' \
                % (sample_name, sample_name, curr_datetime)

    os.system("%s mem -t 2 "  # number of threads
              "-k 18 "  # minimum seed length
              "-aM "  # output alignments for SE and unpaired PE & mark shorter split hits as secondary
              "-R \'%s\' "  # read group header
              "%s "  # genome (bwa indexed hg19)
              "%s %s"  # input R1, input R2
              "| sed \'s/@PG\tID:bwa/@PG\tID:%s/\' - "  # "-" requests standard output, useful when combining tools
              "| %s --removeDups 2> "
              "%s.samblaster.log --splitterFile %s.splitreads.sam --discordantFile %s.discordant.sam "
              "| %s view -Sb - > %s"  # -Sb puts output through samtools sort
              % (config_dict['bwa'], rg_header, config_dict['genome'], fastq1, fastq2, sample_name,
                 config_dict['samblaster'], sample_name, sample_name, sample_name, config_dict['samtools'], outfile))


@follows(align_bwa)
@transform(["*.bwa.drm.bam"], suffix(".bwa.drm.bam"), ".bwa.drm.sorted.bam")
def sort_bam(infile, outfile):
    os.system("%s sort "
              "-@ 2 "  # number of threads = 2
              "-m 2G "  # memory per thread = 2G
              "-O bam "  # output file format = .bam
              "-T %s "  # temporary file name
              "-o %s "  # output file
              "%s"  # input file
              % (config_dict['samtools'], infile, outfile, infile))


@follows(sort_bam)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.bai")
def index_bam(infile, outfile):
    os.system("%s index %s" % (config_dict['samtools'], infile))


@follows(index_bam)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.stats")
def run_samtools_stats(infile, outfile):
    os.system("%s stats %s > %s" % (config_dict['samtools'], infile, outfile))
    os.system("%s -p %s %s" % (config_dict['plot_bamstats'], outfile, outfile))

    joint_sample_name = infile[:-19]
    from create_sample_quality_pdf import CreateFastQCPDF
    pdf = CreateFastQCPDF(joint_sample_name, config_dict['pdflatex'])
    pdf.create_pdf()


@follows(run_samtools_stats)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".breakdancer_config.txt")
def create_breakdancer_config(infile, outfile):
    sample_name = infile[:-19]
    config_file = open("%s" % outfile, "w")  # write into output file
    config_file.write("map:%s\tmean:160\tstd:50\treadlen:150\tsample:%s\texe:bwa-0.7.12\n" % (infile, sample_name))
    config_file.close()


@follows(create_breakdancer_config)
@transform(create_breakdancer_config, suffix(".breakdancer_config.txt"), ".breakdancer_output.sv")
def run_breakdancer(infile, outfile):
    os.system("%s -q 10 %s > %s" % (config_dict['breakdancer'], infile, outfile))


@follows(run_breakdancer)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".pindel_config")
def create_pindel_config(infile, outfile):
    sample_name = infile[:-19]
    config_file = open("%s" % outfile, "w+")  # write into output file
    config_file.write("%s\t240\t%s\n" % (infile, sample_name))  # name of BAM; insert size; sample_label
    config_file.close()


@follows(create_pindel_config)
@transform("*.pindel_config", formatter(), ["{path[0]}/{basename[0]}._BP",
                                            "{path[0]}/{basename[0]}._D",
                                            "{path[0]}/{basename[0]}._INV",
                                            "{path[0]}/{basename[0]}._LI",
                                            "{path[0]}/{basename[0]}._RP",
                                            "{path[0]}/{basename[0]}._SI",
                                            "{path[0]}/{basename[0]}._TD"], "{path[0]}/{basename[0]}.")
def run_pindel(infile, outfile, pindel_prefix):
    bd_file = infile[:-14]
    os.system("%s "  # pindel program
              "-f %s "  # reference genome
              "-i %s "
              "-c ALL "  # ??output_prefix
              "-o %s "  # look at all chromosomes
              "--minimum_support_for_event 5 "
              "-T 2 "  # number of threads
              "-w 50 "  # window size
              "-x 4 "  # maximum size of SVs to detect, 4 = 8092    why 4???
              "--sensitivity 0.95 "
              "-v 6 "  # minimum inversion size in bases     why 6???
              "-e 0.01 "  # expected sequencing error rate
              "--report_breakpoints TRUE "
              "--report_long_insertions TRUE "
              "--name_of_logfile %spindel.log "
              "-b %s.breakdancer_output.sv"  # file name with BreakDancer results     or -Q???
              % (config_dict['pindel'], config_dict['genome'], infile, pindel_prefix, pindel_prefix, bd_file))


# https://github.com/ELIXIR-ITA-training/VarCall2015/blob/master/Tutorials/T5.2_variantcalling_stucturalvariants_tutorial.md
# for pindel2vcf parameters
@follows(run_pindel)
@transform(run_pindel, formatter(), "{path[0]}/{basename[0]}.pindel.merged.vcf", "{path[0]}/{basename[0]}.")
def pindel_to_vcf(infile, outfile, pindel_prefix):
    os.system("%s "  # pindel2vcf program
              "-P %s "  # input file
              "-r %s "  # reference genome on computer
              "-R hg19 "  # reference genome
              "-d 2009 "
              "--min_size 5 "  # minimum size of events (what measurement??)
              "--min_coverage 10 "  # minimum number of reads
              "--min_supporting_reads 5 "
              "--max_size 8092 "
              "-v %s "  # output file
              "--het_cutoff 0.1 "
              "--hom_cutoff 0.85 "
              "-G"
              % (config_dict['pindel2vcf'], pindel_prefix, config_dict['genome'], outfile))


@follows(pindel_to_vcf)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), r"\1.bcf")
def call_delly(infile, outfile):
    sample_name = infile[:-19]
    blacklisted_regions = '%s/excluded_regions.excl.bed' % script_dir
    os.system("%s call -t DEL -x %s -o %s.del.bcf -g %s %s" % (
        config_dict['delly'], blacklisted_regions, sample_name, config_dict['genome'], infile))
    os.system("%s call -t TRA -x %s -o %s.tra.bcf -g %s %s" % (
        config_dict['delly'], blacklisted_regions, sample_name, config_dict['genome'], infile))
    os.system("%s call -t INV -x %s -o %s.inv.bcf -g %s %s" % (
        config_dict['delly'], blacklisted_regions, sample_name, config_dict['genome'], infile))
    os.system("%s call -t INS -x %s -o %s.ins.bcf -g %s %s" % (
        config_dict['delly'], blacklisted_regions, sample_name, config_dict['genome'], infile))
    os.system("%s call -t DUP -x %s -o %s.dup.bcf -g %s %s" % (
        config_dict['delly'], blacklisted_regions, sample_name, config_dict['genome'], infile))


@follows(call_delly)
@transform(["*.bcf"], suffix(".bcf"), r"\1.delly.vcf")
def delly_to_vcf(infile, outfile):
    os.system("%s view %s > %s" % (config_dict['bcftools'], infile, outfile))


@follows(delly_to_vcf)
@collate(["*.vcf"], regex(r"(.+)\.(.+)\.(.+)\.vcf"), r"\1.unified.vcf")
def combine_vcfs(infile, outfile):
    sample_name = outfile[:-12]
    delly_del_data = {}
    delly_dup_data = {}
    delly_inv_data = {}
    delly_ins_data = {}
    delly_tra_data = {}
    pindel_data = {}
    for f in infile:
        if f.endswith(".del.delly.vcf"):
            delly_del_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".dup.delly.vcf"):
            delly_dup_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".inv.delly.vcf"):
            delly_inv_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".ins.delly.vcf"):
            delly_ins_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".tra.delly.vcf"):
            delly_tra_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".pindel.merged.vcf"):
            pindel_data = parse_vcfs.get_pindel_output(f)
    header = ['##fileformat=VCFv4.0\n',
              '##fileDate=%s\n' % curr_datetime,
              '##source=Delly,Pindel\n',
              '##FILTER=<ID=PASS,Description="All filters passed">\n',
              '##INFO=<ID=Precision,Number=1,Type=String,Description="Precise structural variation">\n',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
              '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n',
              '##INFO=<ID=Caller,Number=1,Type=String,Description="Caller used">\n',
              '##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">\n',
              '##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n',
              '##INFO=<ID=Region,Number=1,Type=String,Description="Chromosomal region of variant">\n',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
              '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample_name]
    with open("%s" % outfile, 'w+') as output_file:
        for line in header:
            output_file.write(line)
        for var in delly_del_data.keys():
            if delly_del_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_del_data, var))
        for var in delly_dup_data.keys():
            if delly_dup_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_dup_data, var))
        for var in delly_ins_data.keys():
            if delly_ins_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_ins_data, var))
        for var in delly_inv_data.keys():
            if delly_inv_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_inv_data, var))
        for var in delly_tra_data.keys():
            if delly_tra_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_tra_data, var))
        for var in pindel_data.keys():
            output_file.write(parse_vcfs.print_vcf(pindel_data, var))


@follows(combine_vcfs)
@transform(["*.unified.vcf"], suffix(".unified.vcf"), ".annovar.vcf")
def annotate_vcf(infile, outfile):
    os.system("%s "  # table_annovar.pl
              "%s "  # infile
              "%s "  # directory database files are stored in
              "-buildver hg19 "  # genome build
              "-out %s "  # outfile
              "-remove "  # removes all temporary files
              "-protocol refGene,knownGene,ensgene,esp6500siv2_all,snp138, "  # databases
              "-arg '--exonicsplicing --hgvs --splicing 30',,,, "  # optional arguments, splicing threshold 30bp
              "-operation g,g,g,f,f "  # gene-based or filter-based for protocols
              "-nastring . "  # if variant doesn't have a score from tools (e.g. intronic variant & SIFT), position="."
              "-vcfinput"  # required if input is vcf file
              % (config_dict['annovar'], infile, config_dict['humandb'], outfile))

    os.rename("%s.hg19_multianno.vcf" % outfile, "%s" % outfile)


@follows(annotate_vcf)
@transform(["*.annovar.vcf"], suffix(".annovar.vcf"), ".variant_id.db")
def compress_variants(infile, outfile):
    # shelve object containing variant dictionary for each sample, contained within a database file.
    shelf = shelve.open(outfile, writeback=False)
    variant_dict_for_sample = defaultdict(list)
    sample_name = infile[:-12]

    with open(infile, 'r') as inputfile:
        for line in inputfile:
            line_split = line.strip().split('\t')
            if line_split[0].startswith("#"):
                pass
            else:
                variant_id = "%s_%s_%s_%s" % (line_split[0], line_split[1], line_split[3], line_split[4])
                variant_dict_for_sample[variant_id] = [sample_name]
        shelf.update(variant_dict_for_sample)
        shelf.close()


@follows(compress_variants)
@merge(compress_variants, "master_variants.txt")
def combine_variants(infile, outfile):
    # produces dictionary where each variant is the key and the value is a list containing the worksheet and sample id.
    variant_dict_for_run = defaultdict(list)
    for inputfile in infile:
        sample_variants = shelve.open(inputfile)
        for k, v in sample_variants.items():
            variant_dict_for_run[k].append((worksheet, ','.join(v)))
        sample_variants.close()  # need to close to save

    master_database = shelve.open("AML_variant_frequency.db")
    # only populate if empty to have something to check against, will add on later if not.
    if len(master_database) == 0:
        master_database.update(variant_dict_for_run)
    else:
        pass

    # if already in master database, add sample to existing sample(s) so it doesn't overwrite; otherwise just add.
    combined_dict = defaultdict(list)
    for k, v in variant_dict_for_run.items():
        if k in master_database:
            combined_dict[k] = master_database[k] + v
        else:
            combined_dict[k] = v

    master_database.update(combined_dict)
    master_database.close()

    # Create new dictionary with unique sample numbers to remove duplicate sample entries which may be caused by
    # running the script multiple times.
    analysis_dict = defaultdict(list)
    for k, v in combined_dict.items():
        value_set = set()
        for item in v:
            value_set.add(item)
        value_list = list(value_set)
        analysis_dict[k] = value_list

    shelf = shelve.open("AML_variant_frequency.db", writeback=False)
    shelf.update(analysis_dict)
    shelf.close()


@follows(combine_variants)
@transform(["*.annovar.vcf"], suffix(".annovar.vcf"), ".annovar.final.vcf")
def add_freq_to_vcf(infile, outfile):
    master_database = shelve.open("AML_variant_frequency.db")
    total_samples_list = []
    for k, v in master_database.items():
        total_samples_list.extend([item[1] for item in v])
        total_samples = len(set(total_samples_list))

    with open("%s" % infile, 'r') as input_vcf:
        with open("%s" % outfile, 'w+') as output_vcf:
            for row in input_vcf:
                if row.startswith("#"):
                    output_vcf.write(row)
                else:
                    row_split = row.strip().split('\t')
                    variant = "%s_%s_%s_%s" % (row_split[0], row_split[1], row_split[3], row_split[4])
                    chrom = row_split[0]
                    pos = row_split[1]
                    id_col = row_split[2]
                    ref = row_split[3]
                    alt = row_split[4]
                    qual = row_split[5]
                    filter_col = row_split[6]
                    info = row_split[7]
                    format_col = row_split[8]
                    format_val = row_split[9]
                    if master_database[variant]:
                        total_runs = len(set(item[0] for item in master_database[variant]))
                        total_vars = len(set(item[1] for item in master_database[variant]))
                        pct = round((float(total_vars) / float(total_samples)) * 100, 1)
                        freq = "%s|%s|%s|%s" % (total_vars, total_samples, pct, total_runs)
                    else:
                        freq = "."
                    output_vcf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s;FREQ=%s\t%s\t%s\n"
                                     % (chrom, pos, id_col, ref, alt, qual, filter_col, info, freq, format_col,
                                        format_val))


@follows(add_freq_to_vcf)
@transform(["*.annovar.final.vcf"], suffix(".annovar.final.vcf"), ".annovar.xlsx")
def vcf_to_excel(infile, outfile):
    colors = ['blue', 'red', 'green', 'orange', 'purple', 'black', 'cyan4', 'deeppink', 'seagreen', 'yellow']
    color_no = 0
    con = sql.connect('%s/db.sqlite3' % parent_dir)
    curs = con.cursor()
    curs.execute("CREATE TABLE IF NOT EXISTS Results(result_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, run TEXT, "
                 "sample TEXT, caller TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, chr2 TEXT, end INTEGER, "
                 "region TEXT, sv_type TEXT, size INTEGER, gt TEXT, depth INTEGER, alleles TEXT, ab REAL, gene TEXT, "
                 "func TEXT, exonic_func TEXT, freq TEXT, af REAL, precision TEXT)")
    curs.execute("CREATE TABLE IF NOT EXISTS Runs(run_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, run TEXT)")
    curs.execute("CREATE TABLE IF NOT EXISTS Samples(sample_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, "
                 "sample TEXT, run TEXT)")

    vcf_reader = vcf.Reader(open(infile, 'r'))
    col_list = ['run', 'sample', 'caller', 'chrom', 'pos', 'ref', 'alt', 'chr2', 'end', 'region', 'sv_type', 'size',
                'gt', 'depth', 'alleles', 'ab', 'gene', 'func', 'exonic_func', 'freq', 'af', 'precision']
    annovar_df = pd.DataFrame(columns=col_list)
    sample_name = infile[:-18]
    rcircos_file = open("%s.txt" % sample_name[3:12:], "w")  # write into output file
    for record in vcf_reader:
        for sample in record:
            # PyVCF reader
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = ",".join(str(a) for a in record.ALT)
            gt = sample['GT']
            ad = sample['AD']
            info_dict = record.INFO
            chr2 = info_dict.get("CHR2")
            end_pos = info_dict.get("END")
            sv_type = info_dict.get("SVTYPE")
            gene = ",".join(str(g) for g in info_dict.get("Gene.refGene"))
            if re.match("None", gene.upper()):
                gene_mod = "."
            else:
                gene_mod = gene
            func = ",".join(str(f) for f in info_dict.get("Func.refGene"))
            exonic_func = ",".join(str(f) for f in info_dict.get("ExonicFunc.refGene"))
            ref_reads = ad[0]
            alt_reads = ad[1]
            total_reads = ref_reads + alt_reads
            freq = ",".join(str(a) for a in info_dict.get("FREQ"))
            freq_split = freq.split('|')
            allele_freq = float(freq_split[2])
            if alt_reads != 0 and ref_reads != 0:
                ab = format((float(alt_reads) / float(total_reads) * 100), '.2f')
            else:
                ab = format(int(0), '.2f')
            if sv_type == "DUP:TANDEM":
                size = len(alt) - len(ref)
            else:
                size = info_dict.get("SVLEN")
            ad_str = '%s,%s' % (str(ref_reads), str(alt_reads))
            caller = info_dict.get("Caller")
            prec = ",".join(str(a) for a in info_dict.get("PRECISION"))
            if re.match("None", exonic_func):
                exonic_func_mod = '.'
            else:
                exonic_func_mod = exonic_func
            if re.match("None", func):
                func_mod = '.'
            else:
                func_mod = func
            region = info_dict.get("Region")
            if caller == 'Delly':
                if color_no > 9:
                    color_no = 0
                if prec == 'PRECISE':
                    rcircos_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, pos, chr2, end_pos, end_pos, colors[color_no]))
                    color_no += 1
                output_df = pd.DataFrame([[worksheet, sample_name[3:12:], caller, chrom, pos, ref, alt, chr2,
                                           end_pos, region, sv_type, size, gt, total_reads, ad_str, ab, gene_mod,
                                           func_mod, exonic_func_mod, freq, allele_freq, prec]], columns=col_list)
                annovar_df = annovar_df.append(output_df)
            else:
                if re.match("(.*)exonic(.*)", func) or re.match("(.*)splicing(.*)", func):
                    output_df = pd.DataFrame([[worksheet, sample_name[3:12:], caller, chrom, pos, ref, alt, chr2,
                                               end_pos, region, sv_type, size, gt, total_reads, ad_str, ab, gene_mod,
                                               func_mod, exonic_func_mod, freq, allele_freq, prec]], columns=col_list)
                    annovar_df = annovar_df.append(output_df)
                else:
                    pass

    rcircos_file.close()
    os.system("Rscript %s/rcircos_link.R %s.txt %s.png"
              % (script_dir, sample_name[3:12:], sample_name[3:12:]))

    writer = ExcelWriter('%s' % outfile)
    annovar_df.to_excel(writer, sheet_name="Variants-all", index=False)
    writer.save()

    curs.execute("SELECT * FROM Runs WHERE run LIKE '%s'" % worksheet)
    run_result = curs.fetchone()
    if run_result is None:
        curs.execute("INSERT INTO Runs (run) VALUES (?)", (worksheet,))
    curs.execute("SELECT * FROM Samples WHERE sample LIKE '%s' and run LIKE '%s'" % (sample_name[3:12:], worksheet))
    sample_result = curs.fetchone()
    if sample_result is None:
        curs.execute("INSERT INTO Samples (sample, run) VALUES (?,?)", (sample_name[3:12:], worksheet))
    annovar_df.to_sql("Results", con=con, if_exists='append', index=False)
    con.commit()

    sort_pipeline_data(script_dir, config_dict['static_link_folder'], worksheet, sample_name, args.output_dir)

pipeline_run()
