#!/bin/bash
#
#SBATCH --job-name=MK_snakemake
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err
#SBATCH --mail-type=ALL
#SBATCH --time=30-00:00:00

#module load miniconda
module unload python3
module load snakemake

#snakemake -n
snakemake --snakefile Snakefile --executor cluster-generic\
                                --cluster-generic-submit-cmd "sbatch --mem {resources.mem_mb}\
                                                  --cpus-per-task {threads}\
                                                  --time={resources.time}\
                                                  --output=log/slurm_{rule}-%A.out\
                                                  --error=log/slurm_{rule}-%A.err"\
                                --jobs 40\
                                --cores 5000\
                                --keep-going #-R diamond_search
