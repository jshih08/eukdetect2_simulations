import pandas as pd
from pathlib import Path
import numpy as np
import glob

# get input files
data_path = Path("/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads")

# SAMPLES = [f.name.replace("_1.fastq.gz", "") for f in data_path.glob("*_1.fastq.gz") if f.is_file()]
SAMPLES = ["ERR4097207","ERR4097206","ERR4097205","ERR4097204","ERR4097171","ERR4097175"]

print(SAMPLES)
rule target:
    input:
        expand("zymo_results/eukdetect/{sample}/{sample}_filtered_hits_eukfrac.txt", sample = SAMPLES),
        expand("zymo_results/kraken2_pluspf/{sample}.report.txt", sample = SAMPLES),
        # expand("zymo_results/eukdetect/{sample}/{sample}_normalized.tsv", sample = SAMPLES),
        expand("zymo_results/metaphlan/{sample}_profile.txt", sample = SAMPLES),
        expand("zymo_results/metaphlan_relaxed/{sample}_profile.txt", sample = SAMPLES)

rule eukdetect2:
    input:
        r1 = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_1.fastq.gz",
        r2 = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_2.fastq.gz"
    output:
        ed_out = "zymo_results/eukdetect/{sample}/{sample}_filtered_hits_eukfrac.txt",
    params:
        outdir = "zymo_results/eukdetect/{sample}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/eukdb_v5"
    resources:
        cpus_per_task = 20,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        eukdetect single \
            -1 {input.r1} \
            -2 {input.r2} \
            -n {wildcards.sample} \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --cores {resources.cpus_per_task} \
            --force || true
        """

rule eukdetect2_norm:
    input:
        ed_out = "zymo_results/eukdetect/{sample}/{sample}_filtered_hits_eukfrac.txt",
        r1 = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_1.fastq.gz",
        r2 = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_2.fastq.gz"
    output:
        norm_out = "zymo_results/eukdetect/{sample}/{sample}_normalized.tsv"
    resources:
        cpus_per_task = 20,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    params:
        outdir = "zymo_results/eukdetect/{sample}/",
        eukdb = "/storage/project/r-alind6-0/jshih41/tools/eukdetect2_db"
    shell:
        """
        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "{wildcards.sample}\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {input.ed_out} \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {output.norm_out} || true

        """

rule kraken2_classify_pluspf:
    input:
        fwd = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_1.fastq.gz",
        rev = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_2.fastq.gz",
        db = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/k2_pluspf_db"
    output:
        report = "zymo_results/kraken2_pluspf/{sample}.report.txt"
    params:
        output = "zymo_results/kraken2_pluspf"
    resources:
        cpus_per_task = 20,
        mem_mb = 128000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        k2 classify --db {input.db} \
            --threads {resources.cpus_per_task} \
            --paired \
            --classified-out {params.output}/{wildcards.sample}_classified#.fastq \
            --output {params.output}/{wildcards.sample}.kraken     \
            --report {output.report} \
            --use-names \
            {input.fwd} {input.rev}
        """
        # rm {params.output}


rule run_metaphlan:
    input:
        fwd = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_1.fastq.gz",
        rev = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_2.fastq.gz"
    output:
        profile = "zymo_results/metaphlan/{sample}_profile.txt",
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        metaphlan {input.fwd},{input.rev} \
            --mapout zymo_results/metaphlan/{wildcards.sample}_mapout.txt.bz2 \
            --input_type fastq \
            --db_dir /storage/home/hcoda1/1/jshih41/r-alind6-0/tools/metaphlan_db \
            --output_file {output.profile} \
            --nproc {resources.cpus_per_task}
        """

rule run_metaphlan_relaxed:
    input:
        fwd = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_1.fastq.gz",
        rev = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/zymo_reads/{sample}_2.fastq.gz"
    output:
        profile = "zymo_results/metaphlan_relaxed/{sample}_profile.txt",
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        metaphlan {input.fwd},{input.rev} \
            --mapout zymo_results/metaphlan_relaxed/{wildcards.sample}_mapout.txt.bz2 \
            --input_type fastq \
            --db_dir /storage/home/hcoda1/1/jshih41/r-alind6-0/tools/metaphlan_db \
            --output_file {output.profile} \
            --stat_q 0 \
            --nproc {resources.cpus_per_task}
        """
