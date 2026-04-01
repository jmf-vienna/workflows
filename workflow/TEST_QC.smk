#This script will take interleave paired-end DNA libraries from .bam files to interleaved .fastq.gz files

#grab names for samples from data directory
FILES = glob_wildcards('data/DNA/{name}.bam')
NAMES = FILES.name


#Request all necessary outputs. For this workflow these are a interleave .fastq.gz and corresponding fastqc files
rule all:
    input:
        "results/coassembly/final_bins/data_tables/genomeInformation.csv"
        
        #expand("results/{sample}_spades/final_bins/data_tables/genomeInformation.csv", sample=NAMES)

#DNA BAM to interleaved reads
include: "rules/TEST_link.smk"
include: "rules/TEST_bamtofq.smk"
include: "rules/TEST_fastp.smk"
include: "rules/TEST_multiqc.smk"
include: "rules/TEST_mvdna.smk"

#assembly
include: "rules/TEST_spades.smk"
include: "rules/TEST_metabatnocov.smk"
include: "rules/TEST_metabat_cov.smk"
include: "rules/TEST_drep.smk"
#include: "rules/TEST_gtdbtk.smk"

#coassembly
include: "rules/TEST_megahit.smk"
include: "rules/TEST_co_metabat_cov.smk"
include: "rules/TEST_co_metabatnocov.smk"
include: "rules/TEST_co_drep.smk"
