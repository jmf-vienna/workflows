#This script will take interleave paired-end DNA libraries from .bam files to interleaved .fastq.gz files
configfile: "config/general_configs.yaml"

#grab names for samples from data directory
FILES = glob_wildcards('data/DNA/{name}.1.fastq.gz')
NAMES = FILES.name


#Request all necessary outputs. For this workflow these are a interleave .fastq.gz and corresponding fastqc files
rule all:
    input:
        "intermediates/samples.txt",
        "intermediates/sylph_abundance.tsv",
        "intermediates/singleM_microbialfrac.tsv"


#globDB sylph
rule sylph:
    input:
        "intermediates/{sample}.1.phiXclean.fastq.gz",
        "intermediates/{sample}.2.phiXclean.fastq.gz"
    output:
        "intermediates/{sample}.1.phiXclean.fastq.gz.sylphmpa"
    conda:
        "envs/sylph.yaml"
    threads: 8
    resources: mem_mb=100000, runtime="1d", partition="basic"
    log: "log/sylph_{sample}.log"
    shell:
        """
        sylph profile {config[SYLPHDB]} -o intermediates/{wildcards.sample}.sylph.tsv -u --read-seq-id 99 -t {threads} -1 {input[0]} -2 {input[1]} 2> {log}
        sylph-tax taxprof intermediates/{wildcards.sample}.sylph.tsv -t GlobDB_r226 -o intermediates/ 2> {log}
        """

rule sylph_tax:
    input:
        expand("intermediates/{sample}.1.phiXclean.fastq.gz.sylphmpa", sample=NAMES)
    output:
        "intermediates/sylph_abundance.tsv"
    conda:
        "envs/sylph.yaml"
    threads: 8
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        sylph-tax merge intermediates/*.sylphmpa --column sequence_abundance -o intermediates/sylph_abundance.tsv
        """


#globDB singleM
rule singleM:
    input:
        "intermediates/{sample}.1.phiXclean.fastq.gz",
        "intermediates/{sample}.2.phiXclean.fastq.gz"
    output:
        "intermediates/{sample}.spf.tsv"
    conda:
        "envs/singleM.yaml"
    threads: 16
    resources: mem_mb=100000, runtime="1d", partition="basic"
    log: "log/singlem_{sample}.log"
    shell:
        """
        singlem pipe --threads 16 -1 {input[0]} -2 {input[1]} --metapackage {config[SINGLEMDB]} --taxonomic-profile intermediates/{wildcards.sample}.profile 2> {log}

        touch {output}

        singlem microbial_fraction --metapackage {config[SINGLEMDB]} --forward {input[0]} --reverse {input[1]} -p intermediates/{wildcards.sample}.profile > {output} 2> {log}
        """

#singleM merge
rule singleM_merge:
    input:
        expand("intermediates/{sample}.spf.tsv", sample=NAMES)
    output:
        "intermediates/singleM_microbialfrac.tsv"
    threads: 1
    resources: mem_mb=10000, runtime="1h", partition="basic"
    shell:
        """
        for f in intermediates/*spf.tsv; do head -1 $f > intermediates/singleM_microbialfrac.tsv; done
        for f in intermediates/*spf.tsv; do tail -1 $f >> {output} ; done
        """


#nonpareil
rule nonpareil:
    input:
        "intermediates/{sample}.1.phiXclean.fastq.gz"
    output:
        "intermediates/{sample}.npo"
    conda:
        "envs/nonpareil.yaml"
    threads: 16
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        nonpareil -t 16 -s {input} -T kmer -f fastq -b intermediates/{wildcards.sample}
        """


#merge and generate things
rule nonpareil_stats_gen:
    input:
        expand("intermediates/{sample}.npo", sample=NAMES)
    output:
        "intermediates/samples.txt"
    conda:
        "envs/R.yaml"
    threads: 1
    resources: mem_mb=50000, runtime="4h", partition="basic"
    shell:
        """
        echo "FileTABName" > intermediates/samples.txt
        for f in intermediates/*npo; do g=${{f##*/}}; echo "${{f}}TAB${{g%%.npo}}" >> intermediates/samples.txt; done
        sed -i 's/TAB/\t/g' intermediates/samples.txt

        Rscript workflow/scripts/Nonpareil.R
        """
#add R

#merge all things with python for fun

