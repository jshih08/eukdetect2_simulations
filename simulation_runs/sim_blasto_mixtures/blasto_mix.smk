import pandas as pd
from pathlib import Path
import numpy as np
import glob

# get input files

genome_map = pd.read_csv("blasto_mix_mapping.csv")
genome_map["mix_id"] = genome_map["mix_id"].astype(str)

GENOMES = genome_map["full_acc"].astype(str).unique().tolist()
MIX_ID = genome_map["mix_id"].unique().tolist()
fcov_by_mix_genome = genome_map.set_index(["mix_id", "full_acc"])["sim_fcov"].to_dict()
GENOME_DIR = Path("/storage/home/hcoda1/1/jshih41/scratch/euk_genomes")

def samples_in_batch(mix):
    genomes = genome_map.loc[
        genome_map["mix_id"].astype(str) == str(mix), "full_acc"
    ].astype(str).tolist()
    return genomes

def genome_fasta(wildcards):
    genome = wildcards.genome
    filename = genome if genome.endswith(".fna") else f"{genome}.fna"
    return GENOME_DIR / filename

# MIX_ID = ["1"]

rule target:
    input:
        expand("results/combined_reads/fcov_mix/{mix}_simulated_community_R1.fastq", mix = MIX_ID),
        expand("results/combined_reads/fcov_mix/{mix}_simulated_community_R2.fastq", mix = MIX_ID),
        # expand("results/kraken2_pluspf/{mix}.report.txt", mix = MIX_ID),
        # expand("results/metaphlan/{mix}_profile.txt", mix = MIX_ID),
        expand("results/eukdetect/fcov_mix/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt",  mix = MIX_ID)



# # combining the eukaryotic reads witht he crd mixture of reads
# # there should be one sample per eukaryote?


rule sim_all_reads:
    input:
        genome = genome_fasta
    output:
        r1 = "results/simulated_reads/{mix}/{genome}_R1.fastq",
        r2 = "results/simulated_reads/{mix}/{genome}_R2.fastq"
    params:
        fcov = lambda wildcards: fcov_by_mix_genome[(wildcards.mix, wildcards.genome)],
        outprefix = "results/simulated_reads/{mix}/{genome}"
    resources:
        cpus_per_task = 4,
        mem_mb = 16000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        mkdir -p results/simulated_reads/{wildcards.mix}
        art_modern \
            --mode wgs \
            --lc pe \
            --i-seed 123 \
            --parallel 2 \
            --o-fasta-num_threads 2 \
            --builtin_qual_file HiSeq2500_150bp \
            --pe_frag_dist_mean 300 \
            --pe_frag_dist_std_dev 2 \
            --read_len_1 150 --read_len_2 150 \
            --i-file {input.genome} \
            --o-fastq {params.outprefix}.fastq \
            --i-fcov {params.fcov}

        seqtk seq -1 {params.outprefix}.fastq > {output.r1}
        seqtk seq -2 {params.outprefix}.fastq > {output.r2}
        
        rm {params.outprefix}.fastq
        """

# # combining the eukaryotic reads witht he crd mixture of reads
# # there should be one sample per eukaryote?

rule combine_reads:
    input:
        r1_in = lambda wildcards: expand(
            "results/simulated_reads/{mix}/{genome}_R1.fastq",
            mix=wildcards.mix, genome=samples_in_batch(wildcards.mix)),
        r2_in = lambda wildcards: expand(
            "results/simulated_reads/{mix}/{genome}_R2.fastq",
            mix=wildcards.mix, genome=samples_in_batch(wildcards.mix))
    output:
        r1 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R2.fastq"
    resources:
        cpus_per_task = 20,
        mem_mb = 64000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        mkdir -p results/combined_reads/fcov_mix
        cat {input.r1_in} > {output.r1}
        cat {input.r2_in} > {output.r2}
        """


rule eukdetect2:
    input:
        r1 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R2.fastq"
    output:
        ed_out = "results/eukdetect/fcov_mix/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/eukdetect/fcov_mix/{mix}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/eukdb_v5"
    resources:
        cpus_per_task = 12,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        eukdetect single \
            -1 {input.r1} \
            -2 {input.r2} \
            -n simulated_community \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --cores {resources.cpus_per_task} \
            --force

        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv
        """


rule kraken2_classify_pluspf:
    input:
        r1 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R2.fastq",
        db = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/k2_pluspf_db"  
    output:
        report = "results/kraken2_pluspf/{mix}.report.txt"
    params:
        output = "results/kraken2_pluspf"
    resources:
        cpus_per_task = 20,
        mem_mb = 128000,
        runtime = "16h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        k2 classify --db {input.db} \
            --threads {resources.cpus_per_task} \
            --paired \
            --output {params.output}/{wildcards.mix}.kraken     \
            --report {output.report} \
            --use-names \
            {input.r1} {input.r2}
        """


rule run_metaphlan:
    input:
        r1 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/fcov_mix/{mix}_simulated_community_R2.fastq"
    output:
        profile = "results/metaphlan/{mix}_profile.txt",
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        metaphlan {input.r1},{input.r2} \
            --mapout results/metaphlan/{wildcards.mix}_mapout.txt.bz2 \
            --input_type fastq \
            --db_dir /storage/home/hcoda1/1/jshih41/r-alind6-0/tools/metaphlan_db \
            --output_file {output.profile} \
            --nproc {resources.cpus_per_task}
        """
