#!/bin/bash
#
#SBATCH --job-name=run_snakemake
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err
#SBATCH --time=30-00:00:00

#module load miniconda
module load Conda
conda activate /lisc/scratch/jmf/conda_envs/snakemake

snakemake --snakefile $1 --executor slurm \
                                --jobs 20 \
				--use-conda \
				--cores 200 \
				--conda-frontend conda \
				--conda-prefix /lisc/scratch/jmf/internal/Analyses_Jay/snakemake_envs
                                #--keep-going
