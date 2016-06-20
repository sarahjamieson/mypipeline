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
curr_datetime = datetime.datetime.now().isoformat()

# ----------------------------------------------------------------------------------------------------------------------
# Command line parameters for pipeline script:
#   -s = Sample Sheet
#   -d = Directory with results (InterOp and Data folders)
#   -o = Name of output directory
# ----------------------------------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Runs pipeline for Illumina MiSeq Nextera AML data.",
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-s', '--sheet', action="store", dest='sample_sheet', help='An illumina sample sheet for this run.',
                    required=True, default='/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/SampleSheet.csv')
parser.add_argument('-d', '--result_dir', action='store', dest='result_dir', required=True,
                    help='Directory containing Illumina MiSeq InterOp and Data folders.',
                    default='/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/')
parser.add_argument('-o', '--output_dir', action='store', dest='output_dir', required=True,
                    help='Directory for output to be stored in',
                    default='/media/sf_sarah_share/MiSeq_Nextera_Results/')
args = parser.parse_args()

# ----------------------------------------------------------------------------------------------------------------------
# Paths to all tools used in pipeline. Put in config file?
# ----------------------------------------------------------------------------------------------------------------------
Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'
trim_adapters = '/home/cuser/programs/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'
InterOp = '/%s/InterOp/' % args.result_dir
bwa = '/home/cuser/programs/bwa/bwa'
samblaster = '/home/cuser/programs/samblaster/samblaster'
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
plot_bamstats = '/home/cuser/programs/samtools/samtools/bin/plot-bamstats'
bam2cfg = '/home/cuser/programs/breakdancer/breakdancer-max1.4.5/bam2cfg.pl'
breakdancer = '/home/cuser/programs/breakdancer/breakdancer-max'
pindel = '/home/cuser/programs/pindel/pindel'
genome = '/media/genomicdata/ucsc_hg19/ucsc.hg19.fasta'
pindel2vcf = '/home/cuser/programs/pindel/pindel2vcf'
annovar = '/media/sf_sarah_share/ANNOVAR/annovar/table_annovar.pl'
varscan = '/home/cuser/programs/VarScan2/VarScan.v2.4.0.jar'
fastqc = '/home/cuser/programs/FastQC/fastqc'
delly = '/home/cuser/programs/delly_v0.7.3/delly'
bcftools = '/home/cuser/programs/samtools/bcftools-1.3.1/bcftools'

# ----------------------------------------------------------------------------------------------------------------------
# 1) Copy fastqc files to current directory.
# 2) Parse sample sheet to get run and sample info.
# 3) Make new directory for results in user-inputted output directory.
# ----------------------------------------------------------------------------------------------------------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
os.system("cp %s/Data/Intensities/BaseCalls/*.fastq.gz %s/" % (args.result_dir, script_dir))
parse_sheet = ParseSampleSheet(args.sample_sheet)
run_dict, sample_dict = parse_sheet.parse_sample_sheet()
worksheet = run_dict.get('worksheet')
os.system("mkdir %s%s" % (args.output_dir, worksheet))


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
                    corintmetrics, worksheet)
    pdf.create_pdf()

    os.system('cp %s_InterOp_Results.pdf %s%s/' % (worksheet, args.output_dir, worksheet))


# ----------------------------------------------------------------------------------------------------------------------
# Run pipeline:
#   1) Run FastQC on raw fastqc files.
#   2) Run Trimmomatic [.fastq.gz --> .qfilter.fastq.gz]
#   3) Run FastQC on trimmed fastqc files.
#   4) Run BWA to align reads to hg19 [.qfilter.fastq.gz --> .bwa.drm.bam].
#   5) Run SAMtools "Sort" on aligned data [.bwa.drm.bam --> .bwa.drm.sorted.bam].
#   6) Run SAMtools "Index" on sorted data to generate index file [.bwa.drm.sorted.bam --> .bwa.drm.sorted.bam.bai].
#   7) Run SAMtools "Stats" and "Plot-bamstats" & combine with FastQC results.
#   8)
# ----------------------------------------------------------------------------------------------------------------------
@collate("*.fastq.gz", formatter("([^/]+)R[12]_001.fastq.gz$"), "{path[0]}/{1[0]}.fastq.gz")
def run_fastqc(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (fastqc, fastq1))
    os.system("%s --extract %s" % (fastqc, fastq2))


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
              % (Trimmomatic, fastq1, fastq2, fastq1[:-9], fastq1, fastq2[:-9], fastq2, trim_adapters))

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
    os.system("%s --extract %s" % (fastqc, fastq1))
    os.system("%s --extract %s" % (fastqc, fastq2))


@follows(run_fastqc_trimmed)
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12]_001.qfilter.fastq.gz$"), "{path[0]}/{1[0]}.bwa.drm.bam")
def align_bwa(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    name = fastq1[:-20]
    rg_header = '@RG\tID:%s\tSM:%s\tCN:WMRGL\tDS:TruSight_Myeloid_Nextera_v1.1\tDT:%s' % (name, name, curr_datetime)

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
              % (bwa, rg_header, genome, fastq1, fastq2, name, samblaster, name, name, name, samtools, outfile))


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
              % (samtools, infile, outfile, infile))


@follows(sort_bam)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.bai")
def index_bam(infile, outfile):
    os.system("%s index %s" % (samtools, infile))


@follows(index_bam)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.stats")
def run_samtools_stats(infile, outfile):
    os.system("%s stats %s > %s" % (samtools, infile, outfile))
    os.system("%s -p %s %s" % (plot_bamstats, outfile, outfile))

    joint_sample_name = infile[:-19]
    from create_sample_quality_pdf import CreateFastQCPDF
    pdf = CreateFastQCPDF(joint_sample_name)
    pdf.create_pdf()


con = sql.connect('/home/cuser/PycharmProjects/django_apps/mypipeline/db.sqlite3')
curs = con.cursor()

curs.execute("CREATE TABLE IF NOT EXISTS Results(result_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, run TEXT, "
             "sample TEXT, caller TEXT, chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, chr2 TEXT, end INTEGER, "
             "region TEXT, sv_type TEXT, size INTEGER, gt TEXT, depth INTEGER, alleles TEXT, ab REAL, gene TEXT, "
             "func TEXT, exonic_func TEXT)")
curs.execute("CREATE TABLE IF NOT EXISTS Runs(run_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, run TEXT)")
curs.execute("CREATE TABLE IF NOT EXISTS Samples(sample_id INTEGER PRIMARY KEY AUTOINCREMENT NOT NULL, sample TEXT, "
             "run TEXT")


# along with adding to excel spreadsheet
run = ''
if curs.execute("SELECT * FROM Runs WHERE run LIKE '%s'" % run) is None:
    curs.execute("INSERT INTO Runs (run) VALUES (?)", (run,))
else:
    pass
sample = ''
if curs.execute("SELECT * FROM Samples WHERE sample LIKE '%s'" % sample) is None:
    curs.execute("INSERT INTO Samples (sample) VALUES (?)", (sample,))
else:
    pass
'''
curs.execute("INSERT INTO Results (run, sample, caller, chrom, pos, ref, alt, chr2, end, region, sv_type, size, gt, "
             "depth, alleles, ab, gene, func, exonic_func) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
             (run, sample, caller, chrom, pos, ref, alt, chr2, end, region, sv_type, size, gt, depth, alleles, ab,
              gene, func, exonic_func))
con.commit()
'''