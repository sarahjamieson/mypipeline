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
import pandas as pd
os.environ['DJANGO_SETTINGS_MODULE'] = 'mypipeline.settings'
django.setup()
from aml.models import Results
import datetime
import matplotlib.pyplot as plt
import glob
from aml.models import FragmentAnalysis, Results

'''
samples = ['D15-18331', 'D15-21584', 'D15-22373', 'D15-20343', 'D15-25430', 'D03-21521', 'D15-08791', 'D15-08798',
           'D15-04183', 'D15-00899', 'D15-50424', 'D15-45066', 'D14-45300', 'D15-35262', 'D14-33938', 'D13-42537',
           'D14-16565', 'D15-31492', 'D15-41762', 'D14-30832', 'D14-27112', 'D15-02217', 'D15-26810'
           ]

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

samples = [
    '02-D06-20828-SG-Val2-Enrich3-Run4_S2_L001_',
    '03-D13-25437-PG-Val2-Enrich3-Run4_S3_L001_',
    '04-D13-28660-BH-Val2-Enrich3-Run4_S4_L001_',
    '05-D15-41762-MC-Val2-Enrich3-Run4_S5_L001_',
    '06-D16-17437-MY-Val2-Enrich3-Run4_S6_L001_',
    '07-D15-37942-GH-Val2-Enrich3-Run4_S7_L001_',
    '08-D14-31377-LM-Val2-Enrich3-Run4_S8_L001_'
]
'''
samples = ['D03-21521', 'D06-20828']
cols = ['SAMPLE', 'FRAG ITD SIZE', 'FRAG AR', 'NGS ITD SIZE', 'NGS AR']
combined_df = pd.DataFrame(columns=cols)
for sample in samples:
    frag_dict = {}
    ngs_dict = {}
    for item in FragmentAnalysis.objects.filter(sample=sample):
        frag_dict[item.itd] = item.ab
    for item in Results.objects.filter(sample=sample, run='160628_merged', caller='Pindel', gene='FLT3') | Results.objects.filter(sample=sample, run='160805', caller='Pindel', gene='FLT3'):
        ar = float(item.alleles.split(',')[1]) / float(item.alleles.split(',')[0])
        ar = float(format(ar, '.3f'))
        if item.size in ngs_dict:
            ngs_dict[item.size] += ar
        else:
            ngs_dict[item.size] = ar

    for key, value in frag_dict.items():
        if key in ngs_dict:
            ngs_size = key
            ngs_ar = ngs_dict.get(key)
            del ngs_dict[key]
        else:
            ngs_size = "-"
            ngs_ar = "-"
        combined_df_temp = pd.DataFrame([[sample, key, value, ngs_size, ngs_ar]], columns=cols)
        combined_df = combined_df.append(combined_df_temp)
    if ngs_dict:
        for key, value in ngs_dict.items():
            combined_df_temp = pd.DataFrame([[sample, "-", "-", key, value]], columns=cols)
            combined_df = combined_df.append(combined_df_temp)
combined_df = combined_df.reset_index(drop=True)

# for sample in samples, get indexes where SAMPLE == sample, xlsx writer -> merge cells
sample_index_dict = {}
for sample in samples:
    index_list = combined_df[combined_df['SAMPLE'] == sample].index.tolist()
    sample_index_dict[sample] = index_list

writer = pd.ExcelWriter("combined_flt3_results.xlsx", engine="xlsxwriter")
combined_df.to_excel(writer, sheet_name="Sheet1", index=False)
workbook = writer.book
worksheet = writer.sheets["Sheet1"]
merge_format = workbook.add_format({
    'valign': 'vcenter'
})
for key, value in sample_index_dict.items():
    if len(value) > 1:
        start_index = min(value) + 2
        end_index = max(value) + 2
        worksheet.merge_range('A%s:A%s' % (start_index, end_index), key, merge_format)
writer.save()
# works, now just adjust columns and add to app





