import vcf


def get_delly_output(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    delly_dict = {}
    for record in vcf_reader:
        for sample in record:
            info_dict = record.INFO
            if info_dict.get("IMPRECISE") is True:
                precision = 'IMPRECISE'
                ad = "%s,%s" % (sample['DR'], sample['DV'])
            else:
                precision = 'PRECISE'
                ad = "%s,%s" % (sample['RR'], sample['RV'])
            filter_str = ";".join(record.FILTER)
            if filter_str == '':
                filter = 'PASS'
            else:
                filter = 'LowQual'
            alt = ",".join(str(a) for a in record.ALT)
            svlen = info_dict.get("END") - record.POS
            variant_key = "%s_%s_%s_%s" % (record.CHROM, record.POS, record.REF, record.ALT[0])
            if variant_key in delly_dict:
                pass
            else:
                delly_dict[variant_key] = {"Caller": "Delly",
                                           "Chrom": record.CHROM,
                                           "Position": record.POS,
                                           "Ref": record.REF,
                                           "Alt": alt,
                                           "Qual": sample['GQ'],
                                           "Filter": filter,
                                           "Precision": precision,
                                           "SVTYPE": info_dict.get("SVTYPE"),
                                           "SVLEN": svlen,
                                           "CHR2": info_dict.get("CHR2"),
                                           "END": info_dict.get("END"),
                                           "MAPQ": info_dict.get("MAPQ"),
                                           "PE": info_dict.get("PE"),
                                           "GT": sample['GT'],
                                           "AD": ad}
    return delly_dict


def get_pindel_output(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    pindel_dict = {}
    for record in vcf_reader:
        for sample in record:
            info_dict = record.INFO
            alt_temp = ",".join(str(a) for a in record.ALT)
            if alt_temp > 50:
                alt = alt_temp[:50]
            else:
                alt = alt_temp
            if record.QUAL is None:
                qual = "."
            else:
                qual = record.QUAL
            filter_str = ";".join(record.FILTER)
            if filter_str == '':
                filter = "."
            else:
                filter = filter_str
            variant_key = "%s_%s_%s_%s" % (record.CHROM, record.POS, record.REF, record.ALT[0])
            ad = ",".join(str(a) for a in sample['AD'])
            if variant_key in pindel_dict:
                pass
            else:
                pindel_dict[variant_key] = {"Caller": "Pindel",
                                            "Chrom": record.CHROM,
                                            "Position": record.POS,
                                            "Ref": record.REF,
                                            "Alt": alt,
                                            "Qual": qual,
                                            "Filter": filter,
                                            "Precision": ".",
                                            "SVTYPE": info_dict.get("SVTYPE"),
                                            "SVLEN": info_dict.get("SVLEN"),
                                            "CHR2": record.CHROM,
                                            "END": info_dict.get("END"),
                                            "MAPQ": 0,
                                            "PE": 0,
                                            "GT": sample['GT'],
                                            "AD": ad}

    return pindel_dict


def print_vcf(caller_dict, variant):
    variant_in_vcf_format = "%s\t%s\t%s\t%s\t%s\t%s\t%s\tPRECISION=%s;SVTYPE=%s;SVLEN=%s;CHR2=%s;END=%s;MAPQ=%s;" \
                            "Caller=%s;PE=%s\tGT:AD\t%s:%s\n" \
                            % (caller_dict[variant]['Chrom'],
                               caller_dict[variant]['Position'],
                               ".",
                               caller_dict[variant]['Ref'],
                               caller_dict[variant]['Alt'],
                               caller_dict[variant]['Qual'],
                               caller_dict[variant]['Filter'],
                               caller_dict[variant]['Precision'],
                               caller_dict[variant]['SVTYPE'],
                               caller_dict[variant]['SVLEN'],
                               caller_dict[variant]['CHR2'],
                               caller_dict[variant]['END'],
                               caller_dict[variant]['MAPQ'],
                               caller_dict[variant]['Caller'],
                               caller_dict[variant]['PE'],
                               caller_dict[variant]['GT'],
                               caller_dict[variant]['AD'])

    return variant_in_vcf_format
