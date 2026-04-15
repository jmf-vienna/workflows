#!/bin/bash
#
#SBATCH --job-name=run_snakemake
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err
#SBATCH --time=5-00:00:00

#module load miniconda
module load Conda
module load snakemake/9.19.0-foss-2025b
#conda activate /lisc/data/scratch/jmf/conda_envs/snakemake_8_28/

snakemake --snakefile $1 --executor slurm \
                                --jobs 20 \
				--use-conda \
				--cores 200 \
				--conda-frontend conda \
				--conda-prefix /lisc/data/scratch/jmf/internal/Analyses_Jay/snakemake_envs \
				--sdm env-modules
                                #--keep-going
