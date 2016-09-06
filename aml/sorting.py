import os


def sort_pipeline_data(script_dir, static_dir, worksheet, sample, output_dir):
    name = sample[3:12:]

    # 1) Create directory in S: drive where static files are stored
    os.system("mkdir -p %s/%s/%s" % (static_dir, worksheet, name))

    # 2) Move summary quality results to static folders
    os.system("mv %s_InterOp_Results.pdf %s/%s/" % (worksheet, static_dir, worksheet))
    os.system("mv %s.png %s/%s/%s/" % (name, static_dir, worksheet, name))
    os.system("mv %s_sample_quality.pdf %s/%s/%s/" % (name, static_dir, worksheet, name))

    # 3) Move bamstat and fastqc folders to static folders
    os.system("mkdir -p %s/%s/%s/bamstats" % (static_dir, worksheet, name))
    os.system("mkdir -p %s/%s/%s/fastqc_R1" % (static_dir, worksheet, name))
    os.system("mkdir -p %s/%s/%s/fastqc_R1_trim" % (static_dir, worksheet, name))
    os.system("mkdir -p %s/%s/%s/fastqc_R2" % (static_dir, worksheet, name))
    os.system("mkdir -p %s/%s/%s/fastqc_R2_trim" % (static_dir, worksheet, name))

    # 3) Move bamstat images to static folders
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-acgt-cycles.png %s/%s/%s/bamstats/%s.bam.stats-acgt-cycles.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-coverage.png %s/%s/%s/bamstats/%s.bam.stats-coverage.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-gc-content.png %s/%s/%s/bamstats/%s.bam.stats-gc-content.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-gc-depth.png %s/%s/%s/bamstats/%s.bam.stats-gc-depth.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-indel-cycles.png %s/%s/%s/bamstats/%s.bam.stats-indel-cycles.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-indel-dist.png %s/%s/%s/bamstats/%s.bam.stats-indel-dist.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-insert-size.png %s/%s/%s/bamstats/%s.bam.stats-insert-size.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals.png %s/%s/%s/bamstats/%s.bam.stats-quals.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals2.png %s/%s/%s/bamstats/%s.bam.stats-quals2.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals3.png %s/%s/%s/bamstats/%s.bam.stats-quals3.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats-quals-hm.png %s/%s/%s/bamstats/%s.bam.stats-quals-hm.png"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats.html %s/%s/%s/bamstats/%s.bam.stats.html"
        % (sample, static_dir, worksheet, name, name)
    )
    os.system(
        "mv %s.bwa.drm.sorted.bam.stats %s/%s/%s/bamstats/%s.bam.stats"
        % (sample, static_dir, worksheet, name, name)
    )

    # 4) Move fastqc images to static folders
    os.system("mv %sR1_001_fastqc/fastqc_data.txt %s/%s/%s/fastqc_R1/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR1_001_fastqc/summary.txt %s/%s/%s/fastqc_R1/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR1_001_fastqc/Images/adapter_content.png %s/%s/%s/fastqc_R1/%s_adapter_content.png"
              % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR1_001_fastqc/Images/per_base_quality.png %s/%s/%s/fastqc_R1/%s_per_base_quality.png"
              % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR1_001_fastqc/Images/per_sequence_gc_content.png %s/%s/%s/fastqc_R1/"
              "%s_per_sequence_gc_content.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR1_001_fastqc/Images/sequence_length_distribution.png %s/%s/%s/fastqc_R1/"
              "%s_sequence_length_distribution.png" % (sample, static_dir, worksheet, name, name))

    os.system("mv %sR1_001.qfilter_fastqc/fastqc_data.txt %s/%s/%s/fastqc_R1_trim/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR1_001.qfilter_fastqc/summary.txt %s/%s/%s/fastqc_R1_trim/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/adapter_content.png %s/%s/%s/fastqc_R1_trim/"
              "%s_adapter_content.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/per_base_quality.png %s/%s/%s/fastqc_R1_trim/"
              "%s_per_base_quality.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/per_sequence_gc_content.png %s/%s/%s/fastqc_R1_trim/"
              "%s_per_sequence_gc_content.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR1_001.qfilter_fastqc/Images/sequence_length_distribution.png %s/%s/%s/fastqc_R1_trim/"
              "%s_sequence_length_distribution.png" % (sample, static_dir, worksheet, name, name))

    os.system("mv %sR2_001_fastqc/fastqc_data.txt %s/%s/%s/fastqc_R2/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR2_001_fastqc/summary.txt %s/%s/%s/fastqc_R2/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR2_001_fastqc/Images/adapter_content.png %s/%s/%s/fastqc_R2/%s_adapter_content.png"
              % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR2_001_fastqc/Images/per_base_quality.png %s/%s/%s/fastqc_R2/%s_per_base_quality.png"
              % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR2_001_fastqc/Images/per_sequence_gc_content.png %s/%s/%s/fastqc_R2/"
              "%s_per_sequence_gc_content.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR2_001_fastqc/Images/sequence_length_distribution.png %s/%s/%s/fastqc_R2/"
              "%s_sequence_length_distribution.png" % (sample, static_dir, worksheet, name, name))

    os.system("mv %sR2_001.qfilter_fastqc/fastqc_data.txt %s/%s/%s/fastqc_R2_trim/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR2_001.qfilter_fastqc/summary.txt %s/%s/%s/fastqc_R2_trim/"
              % (sample, static_dir, worksheet, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/adapter_content.png %s/%s/%s/fastqc_R2_trim/"
              "%s_adapter_content.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/per_base_quality.png %s/%s/%s/fastqc_R2_trim/"
              "%s_per_base_quality.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/per_sequence_gc_content.png %s/%s/%s/fastqc_R2_trim/"
              "%s_per_sequence_gc_content.png" % (sample, static_dir, worksheet, name, name))
    os.system("mv %sR2_001.qfilter_fastqc/Images/sequence_length_distribution.png %s/%s/%s/fastqc_R2_trim/"
              "%s_sequence_length_distribution.png" % (sample, static_dir, worksheet, name, name))

    # 6) Create symbolic link to worksheet folder into static project directory
    os.system("ln -s %s/%s/ %s/static/aml/%s" % (static_dir, worksheet, script_dir, worksheet))

    # 7) Move additional files to shared folder sf_sarah_share
    os.system("mkdir -p %s%s/%s/Data/" % (output_dir, worksheet, sample))
    os.system("mkdir -p %s%s/%s/Results/" % (output_dir, worksheet, sample))
    os.system("mv %s.annovar.xlsx %s%s/%s/Results/" % (sample, output_dir, worksheet, sample))
    os.system("mv %s* %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR1_001.qfilter_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR1_001_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR2_001.qfilter_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("mv %sR2_001_fastqc/ %s%s/%s/Data/" % (sample, output_dir, worksheet, sample))
    os.system("rm %s.txt" % name)

