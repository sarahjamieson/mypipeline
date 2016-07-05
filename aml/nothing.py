import sqlite3 as sql
from ruffus import *
import vcf
import os

'''
worksheet = '160628_merged'

samples = ['D15-18331', 'D15-21584', 'D15-22373', 'D15-20343', 'D15-25430', 'D03-21521', 'D15-08791', 'D15-08798',
           'D15-04183', 'D15-00899', 'D15-50424', 'D15-45066', 'D14-45300', 'D15-35262', 'D14-33938', 'D13-42537',
           'D14-16565', 'D15-31492', 'D15-41762', 'D14-30832', 'D14-27112', 'D15-02217', 'D15-26810'
           ]

for sample in samples:
    con = sql.connect('/home/cuser/PycharmProjects/django_apps/mypipeline/db.sqlite3')
    curs = con.cursor()
    curs.execute("INSERT INTO Samples (sample, run) VALUES (?,?)", (sample, worksheet))
    con.commit()
'''
script_dir = os.path.dirname(os.path.abspath(__file__))
worksheet = ''


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
    os.system("rcircos_link.R %s %s.png" % (rcircos_file, sample_name[3:12:]))
    os.system("mkdir %s/static/rcircos/%s/" % (script_dir, worksheet))
    os.system("mv %s.png %s/static/rcircos/%s/" % (sample_name[3:12:], script_dir, worksheet))
