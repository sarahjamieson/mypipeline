from scipy import stats
import sqlite3 as sql
from ruffus import *
import vcf
import os
import django
from scipy.optimize import curve_fit
import numpy as np
from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
from create_interop_pdf import CreatePDF

# os.environ['DJANGO_SETTINGS_MODULE'] = 'mypipeline.settings'
# django.setup()
# from aml.models import Results
import datetime
import matplotlib.pyplot as plt
import glob

'''
worksheet = '160628_merged'

samples = ['D15-18331', 'D15-21584', 'D15-22373', 'D15-20343', 'D15-25430', 'D03-21521', 'D15-08791', 'D15-08798',
           'D15-04183', 'D15-00899', 'D15-50424', 'D15-45066', 'D14-45300', 'D15-35262', 'D14-33938', 'D13-42537',
           'D14-16565', 'D15-31492', 'D15-41762', 'D14-30832', 'D14-27112', 'D15-02217', 'D15-26810'
           ]

for sample in samples:
    os.system("cp /media/sf_sarah_share/MiSeq_Nextera_Results/16053/*%s*/Data/*%s*.annovar.final.vcf /home/cuser/PycharmProjects/django_apps/mypipeline/aml/"
              % (sample, sample))


script_dir = os.path.dirname(os.path.abspath(__file__))
worksheet = '16053'


@transform(["*.annovar.final.vcf"], suffix(".annovar.final.vcf"), ".annovar.xlsx")
def vcf_to_excel(infile, outfile):
    vcf_reader = vcf.Reader(open(infile, 'r'))
    sample_name = infile[:-18]
    rcircos_file = open("%s.txt" % sample_name[3:12:], "w")
    for record in vcf_reader:
        for sample in record:
            # PyVCF reader
            chrom = record.CHROM
            pos = record.POS
            info_dict = record.INFO
            chr2 = info_dict.get("CHR2")
            end_pos = info_dict.get("END")
            caller = info_dict.get("Caller")
            prec = ",".join(str(a) for a in info_dict.get("PRECISION"))

            if caller == 'Delly':
                if prec == 'PRECISE':
                    rcircos_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, pos, chr2, end_pos, end_pos))

    rcircos_file.close()
    if os.stat("%s.txt" % sample_name[3:12:]).st_size != 0:
        os.system("Rscript /home/cuser/PycharmProjects/django_apps/mypipeline/aml/rcircos_link.R %s.txt %s.png"
                  % (sample_name[3:12:], sample_name[3:12:]))
    os.system("mkdir -p %s/static/rcircos/%s/" % (script_dir, worksheet))
    os.system("mv %s.png %s/static/rcircos/%s/" % (sample_name[3:12:], script_dir, worksheet))


pipeline_run()


from aml.models import Results
colors = ['blue', 'red', 'green', 'orange', 'purple', 'black', 'cyan4', 'deeppink', 'seagreen', 'yellow']
sample = 'D14-16565'
run = '16053'
color_no = 0
row = Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Delly')
rcircos_file = open("%s.txt" % sample, "w+")
for item in row:
    rcircos_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n"
                       % (item.chrom, item.pos, item.pos, item.chr2, item.end, item.end, colors[color_no]))
    color_no += 1

rcircos_file.close()
os.system("Rscript /home/shjn/PycharmProjects/mypipeline/rcircos_link.R %s.txt %s.png" % (sample, sample))

colors = ['blue', 'red', 'green', 'orange', 'purple', 'black', 'cyan4', 'deeppink', 'seagreen', 'yellow', 'black', 'black', 'black', 'black', 'black']
color_no = 0

samples = ['02-D15-18331-AR-Nextera-Myeloid-Val1-Repeat_S2_L001_',
           '03-D15-21584-RD-Nextera-Myeloid-Val1-Repeat_S3_L001_',
           '04-D15-22373-HT-Nextera-Myeloid-Val1-Repeat_S4_L001_',
           '05-D15-20343-JR-Nextera-Myeloid-Val1-Repeat_S5_L001_',
           '06-D15-25430-BH-Nextera-Myeloid-Val1-Repeat_S6_L001_',
           '07-D03-21521-KT-Nextera-Myeloid-Val1-Repeat_S7_L001_',
           '08-D15-08791-AD-Nextera-Myeloid-Val1-Repeat_S8_L001_',
           '09-D15-08798-MK-Nextera-Myeloid-Val1-Repeat_S9_L001_',
           '10-D15-04183-CS-Nextera-Myeloid-Val1-Repeat_S10_L001_',
           '11-D15-00899-DR-Nextera-Myeloid-Val1-Repeat_S11_L001_',
           '12-D15-50424-LB-Nextera-Myeloid-Val1-Repeat_S12_L001_',
           '13-D15-45066-AG-Nextera-Myeloid-Val1-Repeat_S13_L001_',
           '14-D14-45300-CB-Nextera-Myeloid-Val1-Repeat_S14_L001_',
           '15-D15-35262-GP-Nextera-Myeloid-Val1-Repeat_S15_L001_',
           '16-D14-33938-ES-Nextera-Myeloid-Val1-Repeat_S16_L001_',
           '17-D13-42537-RB-Nextera-Myeloid-Val1-Repeat_S17_L001_',
           '19-D15-31492-AC-Nextera-Myeloid-Val1-Repeat_S19_L001_',
           '20-D15-41762-MC-Nextera-Myeloid-Val1-Repeat_S20_L001_',
           '21-D14-30832-BY-Nextera-Myeloid-Val1-Repeat_S21_L001_',
           '22-D14-27112-BD-Nextera-Myeloid-Val1-Repeat_S22_L001_',
           '23-D15-02217-LT-Nextera-Myeloid-Val1-Repeat_S23_L001_',
           '24-D15-26810-FM-Nextera-Myeloid-Val1-Repeat_S24_L001_'
           ]
for sample in samples:
    os.system("cp /media/sf_sarah_share/MiSeq_Nextera_Results/160725/%s/Data/*.bcf /home/shjn/PycharmProjects/mypipeline/aml/" % sample)
    os.system("cp /media/sf_sarah_share/MiSeq_Nextera_Results/160725/%s/Data/*.delly.vcf /home/shjn/PycharmProjects/mypipeline/aml/" % sample)



xvals = [0.24, 0.17, 0.04, 0.21, 0.77, 0.5, 0.34, 0.42, 0.3, 0.1, 0.36, 0.4, 0.04]
yvals = [0.0643, 0.0583, 0.0193, 0.0461, 0.2825, 0.2625, 0.2095, 0.1669, 0.0272, 0.0347, 0.1658, 0.2261, 0.0134]


# http://stackoverflow.com/questions/19612348/break-x-axis-in-r


def sigmoid(x, x0, k):
    y = 1 / (1 + np.exp(-k*(x-x0)))
    return y

popt, pcov = curve_fit(sigmoid, xvals, yvals)
print popt

x = np.linspace(-1, 15, 50)
y = sigmoid(x, *popt)

plt.plot(xvals, yvals, 'o', label='data')
plt.plot(x, y, label='fit')
plt.ylim(0, 0.3)
plt.xlim(0, 1.0)
plt.legend(loc='best')
plt.show()

from collections import defaultdict

new_dict = {"1": ["val1", "val2"], "2": ["abc", "def"]}
for key, value in new_dict.items():
    for v in new_dict[key]:
        print v

sample = "D15-08798"
run = "16053"
row = Results.objects.filter(sample__icontains=sample, run__icontains=run, caller='Pindel', gene__icontains='FLT3',
                             size__icontains="57")
d = defaultdict(list)
for item in row:
    d[item.size].append(item.pos)
for ke, va in d.items():
    print va

samples = [
    '02-D06-20828-SG-Val2-Enrich3-Run4_S2_L001_',
    '03-D13-25437-PG-Val2-Enrich3-Run4_S3_L001_',
    '04-D13-28660-BH-Val2-Enrich3-Run4_S4_L001_',
    '05-D15-41762-MC-Val2-Enrich3-Run4_S5_L001_',
    '06-D16-17437-MY-Val2-Enrich3-Run4_S6_L001_',
    '07-D15-37942-GH-Val2-Enrich3-Run4_S7_L001_',
    '08-D14-31377-LM-Val2-Enrich3-Run4_S8_L001_'
]
for sample in samples:
    name = sample[3:12:]
    os.system("mv /home/shjn/PycharmProjects/mypipeline/aml/static/aml/160805/%s/%ssample_quality.pdf "
              "/home/shjn/PycharmProjects/mypipeline/aml/static/aml/160805/%s/%s_sample_quality.pdf"
              % (name, name, name, name))
'''
from bs4 import BeautifulSoup
import re
import requests
import pandas as pd
import urllib2

# write chrom location and new base into file
# run curl polyphen-2 for humdiv and humvar and save to two different output files

html_doc = open("output.txt", "r")
soup = BeautifulSoup(html_doc, 'html.parser')
session_id = None

f = file("output.txt").read()
for word in f.split():
    if re.match("polyphenweb2=.{40,};", word):
        session_id = word[13:53]
print session_id

url = "http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/pph2-short.txt" % session_id

opener = urllib2.build_opener()
opener.addheaders = [('User-agent', 'Mozilla/5.0')]
response = opener.open(url)
new_df = pd.read_table(response.read())
print new_df
