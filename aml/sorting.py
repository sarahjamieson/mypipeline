import os


def sort_pipeline_data(script_dir, worksheet, sample, output_dir):
    name = sample[3:12:]

    os.system("mkdir -p %s/static/aml/%s/%s" % (script_dir, worksheet, name))
    os.system('mv %s_InterOp_Results.pdf %s/static/aml/%s/' % (worksheet, script_dir, worksheet))
    os.system('mv %s.png %s/static/aml/%s/%s/' % (name, script_dir, worksheet, name))
    os.system("mv %s_sample_quality.pdf %s/static/aml/%s/%s/" % (name, script_dir, worksheet, name))

    os.system("mkdir -p %s/static/aml/%s/%s/bamstats" % (script_dir, worksheet, name))
    os.system("mkdir -p %s/static/aml/%s/%s/fastqc_R1" % (script_dir, worksheet, name))
    os.system("mkdir -p %s/static/aml/%s/%s/fastqc_R1_trim" % (script_dir, worksheet, name))
    os.system("mkdir -p %s/static/aml/%s/%s/fastqc_R2" % (script_dir, worksheet, name))
    os.system("mkdir -p %s/static/aml/%s/%s/fastqc_R2_trim" % (script_dir, worksheet, name))

    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-acgt-cycles.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-acgt-cycles.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-coverage.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-coverage.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-gc-content.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-gc-content.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-gc-depth.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-gc-depth.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-indel-cycles.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-indel-cycles.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-indel-dist.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-indel-dist.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-insert-size.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-insert-size.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-quals.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals2.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-quals2.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals3.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-quals3.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals-hm.png %s/static/aml/%s/%s/bamstats/%s.bam.stats-quals-hm.png"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats.html %s/static/aml/%s/%s/bamstats/%s.bam.stats.html"
        % (sample, script_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats %s/static/aml/%s/%s/bamstats/%s.bam.stats"
        % (sample, script_dir, worksheet, name, name)
    )

    os.system("mv %sR1_001_fastqc/fastqc_data.txt %s/static/aml/%s/%s/fastqc_R1/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR1_001_fastqc/summary.txt %s/static/aml/%s/%s/fastqc_R1/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR1_001_fastqc/Images/adapter_content.png %s/static/aml/%s/%s/fastqc_R1/%s_adapter_content.png"
              % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR1_001_fastqc/Images/per_base_quality.png %s/static/aml/%s/%s/fastqc_R1/%s_per_base_quality.png"
              % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR1_001_fastqc/Images/per_sequence_gc_content.png %s/static/aml/%s/%s/fastqc_R1/"
              "%s_per_sequence_gc_content.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR1_001_fastqc/Images/sequence_length_distribution.png %s/static/aml/%s/%s/fastqc_R1/"
              "%s_sequence_length_distribution.png" % (sample, script_dir, worksheet, name, name))

    os.system("mv %sR1_001.qfilter_fastqc/fastqc_data.txt %s/static/aml/%s/%s/fastqc_R1_trim/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR1_001.qfilter_fastqc/summary.txt %s/static/aml/%s/%s/fastqc_R1_trim/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/adapter_content.png %s/static/aml/%s/%s/fastqc_R1_trim/"
              "%s_adapter_content.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/per_base_quality.png %s/static/aml/%s/%s/fastqc_R1_trim/"
              "%s_per_base_quality.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/per_sequence_gc_content.png %s/static/aml/%s/%s/fastqc_R1_trim/"
              "%s_per_sequence_gc_content.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/sequence_length_distribution.png %s/static/aml/%s/%s/fastqc_R1_trim/"
              "%s_sequence_length_distribution.png" % (sample, script_dir, worksheet, name, name))

    os.system("mv %sR2_001_fastqc/fastqc_data.txt %s/static/aml/%s/%s/fastqc_R2/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR2_001_fastqc/summary.txt %s/static/aml/%s/%s/fastqc_R2/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR2_001_fastqc/Images/adapter_content.png %s/static/aml/%s/%s/fastqc_R2/%s_adapter_content.png"
              % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR2_001_fastqc/Images/per_base_quality.png %s/static/aml/%s/%s/fastqc_R2/%s_per_base_quality.png"
              % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR2_001_fastqc/Images/per_sequence_gc_content.png %s/static/aml/%s/%s/fastqc_R2/"
              "%s_per_sequence_gc_content.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR2_001_fastqc/Images/sequence_length_distribution.png %s/static/aml/%s/%s/fastqc_R2/"
              "%s_sequence_length_distribution.png" % (sample, script_dir, worksheet, name, name))

    os.system("mv %sR2_001.qfilter_fastqc/fastqc_data.txt %s/static/aml/%s/%s/fastqc_R2_trim/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR2_001.qfilter_fastqc/summary.txt %s/static/aml/%s/%s/fastqc_R2_trim/"
              % (sample, script_dir, worksheet, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/adapter_content.png %s/static/aml/%s/%s/fastqc_R2_trim/"
              "%s_adapter_content.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/per_base_quality.png %s/static/aml/%s/%s/fastqc_R2_trim/"
              "%s_per_base_quality.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/per_sequence_gc_content.png %s/static/aml/%s/%s/fastqc_R2_trim/"
              "%s_per_sequence_gc_content.png" % (sample, script_dir, worksheet, name, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/sequence_length_distribution.png %s/static/aml/%s/%s/fastqc_R2_trim/"
              "%s_sequence_length_distribution.png" % (sample, script_dir, worksheet, name, name))

    os.system("mkdir -p %s%s/%s/Data/" % (output_dir, worksheet, sample))
    os.system("mkdir -p %s%s/%s/Results/" % (output_dir, worksheet, sample))
    os.system("mv %s.annovar.xlsx %s%s/%s/Results/" % (sample, output_dir, worksheet, sample))
    os.system("mv %s* %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR1_001.qfilter_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR1_001_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR2_001.qfilter_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR2_001_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("rm %s.txt" % name)

