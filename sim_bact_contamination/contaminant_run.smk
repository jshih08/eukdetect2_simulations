import pandas as pd

bact_genomes = pd.read_csv("gtdb_genome_map_all_full.csv")
bact_genomes["mix_id"] = bact_genomes["mix_id"].astype(str)
MIX_ID = bact_genomes["mix_id"].unique().tolist()

MIX_ID = [str(i) for i in range(1, 39)]

def samples_in_batch(mix):
    return bact_genomes.loc[bact_genomes["mix_id"] == str(mix), "full_acc"].tolist()


rule target:
    input:
        expand("results/base/eukdetect/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt",mix=MIX_ID)

rule sim_all_reads:
    input:
        genome = "/storage/home/hcoda1/1/jshih41/scratch/gtdb_genomes/all_genomes/{genome}.fna"
    output:
        r1 = "results/base/simulated_reads/{genome}_R1.fastq",
        r2 = "results/base/simulated_reads/{genome}_R2.fastq"
    params:
        # fcov = lambda wildcards: get_params_df(wildcards.sim_config).loc[lambda df: df["Assembly"] == wildcards.genome, "fcov"].values[0],
        outprefix = "results/base/simulated_reads/{genome}"
    resources:
        cpus_per_task = 20,
        mem_mb = 32000,
        runtime = "12h",
        partition = "short",
        slurm_account = "gts-alind6-postpaid"
    shell:
        """
        art_modern \
            --mode wgs \
            --lc pe \
            --i-seed 123 \
            --parallel 20 \
            --o-fasta-num_threads 20 \
            --builtin_qual_file HiSeq2500_150bp \
            --pe_frag_dist_mean 300 \
            --pe_frag_dist_std_dev 2 \
            --read_len_1 150 --read_len_2 150 \
            --i-file {input.genome} \
            --o-fastq results/base/simulated_reads/{wildcards.genome}.fastq \
            --i-fcov 1

        seqtk seq -1 results/base/simulated_reads/{wildcards.genome}.fastq > results/base/simulated_reads/{wildcards.genome}_R1.fastq
        seqtk seq -2 results/base/simulated_reads/{wildcards.genome}.fastq > results/base/simulated_reads/{wildcards.genome}_R2.fastq

        rm results/base/simulated_reads/{wildcards.genome}.fastq
        """

rule combine_reads:
    input:
        r1_all = lambda wildcards: expand(
            "/storage/scratch1/1/jshih41/bact_contam_run/results/base/simulated_reads/{genome}_R1.fastq",
            genome=samples_in_batch(wildcards.mix)
        ),
        r2_all = lambda wildcards: expand(
            "/storage/scratch1/1/jshih41/bact_contam_run/results/base/simulated_reads/{genome}_R2.fastq",
            genome=samples_in_batch(wildcards.mix)
        )
    output:
        r1 = "results/base/combined_reads/{mix}_simulated_community_R1.fastq",
        r2 = "results/base/combined_reads/{mix}_simulated_community_R2.fastq"
    resources:
        cpus_per_task = 4,
        mem_mb = 64000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        cat {input.r1_all} > {output.r1}
        cat {input.r2_all} > {output.r2}
        """


rule eukdetect2:
    input:
        r1 = "results/base/combined_reads/{mix}_simulated_community_R1.fastq",
        r2 = "results/base/combined_reads/{mix}_simulated_community_R2.fastq"
    output:
        ed_out = "results/base/eukdetect/{mix}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/base/eukdetect/{mix}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/eukdb_v5"
    resources:
        cpus_per_task = 20,
        mem_mb = 96000,
        runtime = "8h",
        partition = "short",
        slurm_account = "gts-alind6"
    shell:
        """
        mkdir -p {params.outdir}
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