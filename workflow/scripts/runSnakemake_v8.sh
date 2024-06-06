#!/bin/bash
#
#SBATCH --job-name=run_snakemake
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --output=snakemake.out
#SBATCH --error=snakemake.err
#SBATCH --time=30-00:00:00

#module load miniconda
module load snakemake/8.12.0-3.12.3
module load conda

snakemake --snakefile $1 --executor cluster-generic --cluster-generic-submit-cmd "sbatch --mem {resources.mem_mb}\
                                                  --cpus-per-task {threads}\
						  --partition={resources.partition}\
                                                  --time={resources.time}\
                                                  --output=log/slurm_{rule}-%A.out\
                                                  --error=log/slurm_{rule}-%A.err"\
                                --jobs 20\
				--use-conda
                                #--keep-going #-R diamond_search


