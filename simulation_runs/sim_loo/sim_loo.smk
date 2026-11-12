import pandas as pd
from pathlib import Path
import numpy as np
import glob


# simulate samples of  candida_albicans  candida_dubliniensis  entamoeba_dispar  entamoeba_histolytica at quarterly coverages and run agains each of the loo eukdb databases

# get input files

genome_map = pd.read_csv("sim_loo_mapping.csv")
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


rule target:
    input:
        expand("results/combined_reads/base/{mix}_simulated_community_R1.fastq", mix = MIX_ID),
        expand("results/combined_reads/base/{mix}_simulated_community_R2.fastq", mix = MIX_ID),
        expand("results/cdub_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt",  mix = MIX_ID),
        expand("results/calb_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt",  mix = MIX_ID),
        expand("results/edis_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt",  mix = MIX_ID),
        expand("results/ehist_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt",  mix = MIX_ID)



# # combining the eukaryotic reads with the crd mixture of reads
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
        r1 = "results/combined_reads/base/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/base/{mix}_simulated_community_R2.fastq"
    resources:
        cpus_per_task = 20,
        mem_mb = 64000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        mkdir -p results/combined_reads/base
        cat {input.r1_in} > {output.r1}
        cat {input.r2_in} > {output.r2}
        """


rule eukdetect2_cdub_db:
    input:
        r1 = "results/combined_reads/base/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/base/{mix}_simulated_community_R2.fastq"
    output:
        ed_out = "results/cdub_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/cdub_db/base/{mix}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/leave_one_out_dbs/candida_dubliniensis_with_index_copy"
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
            -n simulated_community \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --database-prefix eukdb_dubliniensis \
            --cores {resources.cpus_per_task} \
            --force || true

        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv || true
        """

rule eukdetect2_calb_db:
    input:
        r1 = "results/combined_reads/base/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/base/{mix}_simulated_community_R2.fastq"
    output:
        ed_out = "results/calb_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/calb_db/base/{mix}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/leave_one_out_dbs/candida_albicans"
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
            -n simulated_community \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --database-prefix  eukdb_albicans \
            --cores {resources.cpus_per_task} \
            --force || true

        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv || true
        """

rule eukdetect2_edis_db:
    input:
        r1 = "results/combined_reads/base/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/base/{mix}_simulated_community_R2.fastq"
    output:
        ed_out = "results/edis_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/edis_db/base/{mix}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/leave_one_out_dbs/entamoeba_dispar"
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
            -n simulated_community \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --database-prefix eukdb_dispar \
            --cores {resources.cpus_per_task} \
            --force || true

        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv || true
        """

rule eukdetect2_ehist_db:
    input:
        r1 = "results/combined_reads/base/{mix}_simulated_community_R1.fastq",
        r2 = "results/combined_reads/base/{mix}_simulated_community_R2.fastq"
    output:
        ed_out = "results/ehist_db/base/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/ehist_db/base/{mix}/",
        eukdb = "/storage/home/hcoda1/1/jshih41/shared/data_swap/johnathan/leave_one_out_dbs/entamoeba_histolytica"
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
            -n simulated_community \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --database-prefix eukdb_histolytica \
            --cores {resources.cpus_per_task} \
            --force || true

        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv || true
        """
