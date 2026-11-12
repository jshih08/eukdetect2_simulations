import pandas as pd
from pathlib import Path
import numpy as np
import glob

input_samples = pd.read_csv("/storage/home/hcoda1/1/jshih41/scratch/final_eukdetect2_sims/simulation_runs/run_baltic_plant/samples.tsv", sep="\t")


SAMPLES = input_samples["sample_name"].tolist()

print(SAMPLES)

rule target:
    input:
        expand("results/eukdetect/{sample}/filtering/simulated_community_read_counts_and_mismatches.txt", sample = SAMPLES)



rule eukdetect2:
    input:
        r1 = "/storage/home/hcoda1/1/jshih41/scratch/globus/{sample}_1.fastq.gz",
        r2 = "/storage/home/hcoda1/1/jshih41/scratch/globus/{sample}_2.fastq.gz"
    output:
        ed_out = "results/eukdetect/{sample}/filtering/simulated_community_read_counts_and_mismatches.txt"
    params:
        outdir = "results/eukdetect/{sample}/",
        eukdb = "/storage/project/r-alind6-0/shared/data_swap/johnathan/eukdb_v5"
    resources:
        cpus_per_task = 20,
        mem_mb = 64000,
        runtime = "6h",
        partition = "short",
        slurm_account = "gts-alind6-postpaid"
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
        
        bases=$(zcat {input.r1} {input.r2} | awk 'NR%4==2{{sum+=length($0)}}END{{print sum}}'); echo -e "simulated_community\t${{bases}}" >> {params.outdir}library_sizes.tsv

        eukdetect-normalize \
            --eukfrac {params.outdir}simulated_community_filtered_hits_eukfrac.txt \
            --library-sizes {params.outdir}library_sizes.tsv \
            --output {params.outdir}simulated_community_normalized.tsv || true
        """