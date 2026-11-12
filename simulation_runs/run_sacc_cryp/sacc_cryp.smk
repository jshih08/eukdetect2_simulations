import pandas as pd
from pathlib import Path
import numpy as np
import glob


# simulate 50/50 of s. cerevisiae and c. neoformans

GENOMES = ["GCF_000146045.2_R64_genomic","GCF_000149245.1_CNA3_genomic"]

SIM_ID = list(map(str, (np.arange(0,3) + 1).tolist()))


rule target:
    input:
        # "results/kraken2_pluspf/sacc_cryp.report.txt",
        expand("results/kraken2_pluspf/{sim}/sacc_cryp.report.txt", sim = SIM_ID),
        expand("results/eukdetect2/{sim}/sacc_cryp/simulated_community_filtered_hits_eukfrac.txt", sim = SIM_ID),
        expand("results/metaphlan/{sim}_profile.txt", sim = SIM_ID),
        expand("results/metaphlan_relaxed/{sim}_profile.txt", sim = SIM_ID)

rule sim_all_reads:
    input:
        genome = "/storage/home/hcoda1/1/jshih41/scratch/euk_genomes/{genome}.fna"
    output:
        r1 = "results/simulated_reads/{sim}/{genome}_R1.fastq",
        r2 = "results/simulated_reads/{sim}/{genome}_R2.fastq"
    params:
        # fcov = lambda wildcards: get_params_df(wildcards.sim_config).loc[lambda df: df["Assembly"] == wildcards.genome, "fcov"].values[0],
    resources:
        cpus_per_task = 8,
        mem_mb = 16000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6-postpaid"
    shell:
        """
        art_modern \
            --mode wgs \
            --lc pe \
            --i-seed {wildcards.sim} \
            --parallel 8 \
            --o-fasta-num_threads 8 \
            --builtin_qual_file HiSeq2500_150bp \
            --pe_frag_dist_mean 300 \
            --pe_frag_dist_std_dev 2 \
            --read_len_1 150 --read_len_2 150 \
            --i-file {input.genome} \
            --o-fastq results/simulated_reads/{wildcards.sim}/{wildcards.genome}.fastq \
            --i-fcov 1

        seqtk seq -1 results/simulated_reads/{wildcards.sim}/{wildcards.genome}.fastq > results/simulated_reads/{wildcards.sim}/{wildcards.genome}_R1.fastq
        seqtk seq -2 results/simulated_reads/{wildcards.sim}/{wildcards.genome}.fastq > results/simulated_reads/{wildcards.sim}/{wildcards.genome}_R2.fastq
        
        rm results/simulated_reads/{wildcards.sim}/{wildcards.genome}.fastq
        """

# combining the eukaryotic reads witht he crd mixture of reads
# there should be one sample per eukaryote?


rule combine_reads:
    input:
        r1_in = expand("results/simulated_reads/{{sim}}/{genome}_R1.fastq", genome = GENOMES),
        r2_in = expand("results/simulated_reads/{{sim}}/{genome}_R2.fastq", genome = GENOMES)
    output:
        r1 = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R1.fastq",
        r2 = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R2.fastq"
    resources:
        cpus_per_task = 4,
        mem_mb = 32000,
        runtime = "4h",
        partition = "short",
        slurm_account = "gts-alind6-postpaid"
    shell:
        """
        cat {input.r1_in} > {output.r1}
        cat {input.r2_in} > {output.r2}
        """

rule eukdetect2:
    input:
        r1 = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R1.fastq",
        r2 = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R2.fastq"
    output:
        ed_out = "results/eukdetect2/{sim}/sacc_cryp/simulated_community_filtered_hits_eukfrac.txt"
    params:
        outdir = "results/eukdetect2/{sim}/sacc_cryp/",
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
            -n simulated_community \
            --outdir {params.outdir} \
            --database {params.eukdb} \
            --cores {resources.cpus_per_task} \
            --force || true

        echo -e "sample_name\ttotal_bases" > {params.outdir}library_sizes.tsv
        
        bases=$(cat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv || true
        """


rule kraken2_classify_pluspf:
    input:
        fwd = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R1.fastq",
        rev = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R2.fastq",
        db = "/storage/home/hcoda1/1/jshih41/r-alind6-0/data/k2_pluspf_db"
    output:
        report = "results/kraken2_pluspf/{sim}/sacc_cryp.report.txt"
    params:
        output = "results/kraken2_pluspf"
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
            --classified-out {params.output}/sacc_cryp_classified#.fastq \
            --output {params.output}/sacc_cryp.kraken     \
            --report {output.report} \
            --use-names \
            {input.fwd} {input.rev}

        """

rule run_metaphlan:
    input:
        fwd = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R1.fastq",
        rev = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R2.fastq",
    output:
        profile = "results/metaphlan/{sim}_profile.txt",
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        metaphlan {input.fwd},{input.rev} \
            --mapout results/metaphlan/{wildcards.sim}_mapout.txt.bz2 \
            --input_type fastq \
            --db_dir /storage/home/hcoda1/1/jshih41/r-alind6-0/tools/metaphlan_db \
            --output_file {output.profile} \
            --nproc {resources.cpus_per_task}
        """

rule run_metaphlan_relaxed:
    input:
        fwd = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R1.fastq",
        rev = "results/combined_reads/{sim}/sacc_cryp_simulated_community_R2.fastq",
    output:
        profile = "results/metaphlan_relaxed/{sim}_profile.txt",
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6-postpaid"
    shell:
        """
        metaphlan {input.fwd},{input.rev} \
            --mapout results/metaphlan_relaxed/{wildcards.sim}_mapout.txt.bz2 \
            --input_type fastq \
            --db_dir /storage/home/hcoda1/1/jshih41/r-alind6-0/tools/metaphlan_db \
            --output_file {output.profile} \
            --stat_q 0 \
            --nproc {resources.cpus_per_task}
        """
