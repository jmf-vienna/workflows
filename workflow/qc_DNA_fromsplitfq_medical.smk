#This script will take interleave paired-end DNA libraries from .bam files to interleaved .fastq.gz files. It will also remove human, mouse, and rat sequences

#grab names for samples from data directory
FILES = glob_wildcards('data/DNA/{name}.1.fastq.gz')
NAMES = FILES.name


#Request all necessary outputs. For this workflow these are a interleave .fastq.gz and corresponding fastqc files
rule all:
    input:
        expand("data/interleave_DNA/{sample}.interleave.fastq.gz", sample=NAMES),
        expand("intermediates/prefilter_qc/{sample}.1_fastqc.html", sample=NAMES),
        expand("intermediates/prefilter_qc/{sample}.2_fastqc.html", sample=NAMES),
        expand("intermediates/postfilter_qc/{sample}_R1.cleaned_fastqc.html", sample=NAMES),
        expand("intermediates/postfilter_qc/{sample}_R2.cleaned_fastqc.html", sample=NAMES),
        "intermediates/prefilter_qc/multiqc_report.html",
        "intermediates/postfilter_qc/multiqc_report.html"

#Fastqc for each fq in intermediates
rule fastqc_prefilter:
    input:
        "data/DNA/{sample}.1.fastq.gz",
        "data/DNA/{sample}.2.fastq.gz"
    output:
        "intermediates/prefilter_qc/{sample}.1_fastqc.html",
        "intermediates/prefilter_qc/{sample}.2_fastqc.html"
    conda: "envs/fastqc.yaml"
    threads: 16
    resources: mem_mb=10000, time="1-00:00:00"
    shell:
        """
        if [ ! -d "intermediates/prefilter_qc" ]; then mkdir intermediates/prefilter_qc; fi
        fastqc -t {threads} -o intermediates/prefilter_qc {input}
        """


rule multiqc_prefilter:
    input:
        expand("intermediates/prefilter_qc/{sample}.1_fastqc.html", sample=NAMES),
        expand("intermediates/prefilter_qc/{sample}.2_fastqc.html", sample=NAMES)
    output:
        "intermediates/prefilter_qc/multiqc_report.html"
    conda: "envs/multiqc.yaml"
    threads: 1
    resources: mem_mb=10000, time="0-00:10:00"
    shell:
        """
        multiqc intermediates/prefilter_qc --outdir intermediates/prefilter_qc -f
        #rm -r intermediates/prefilter_qc/*fastqc*
        """


#bbduk will trimadapters from reads
#we are keeping reads seperate until the end for some processing reasons (fastqc, mostly)
rule bbduk_trimadapters:
    input:
        "data/DNA/{sample}.1.fastq.gz",
        "data/DNA/{sample}.2.fastq.gz"
    output:
        "intermediates/{sample}_R1.adapterclean.fastq.gz",
        "intermediates/{sample}_R2.adapterclean.fastq.gz"
    threads: 16
    conda:
        "envs/bbmap.yaml"
    resources: mem_mb=10000, time="1-00:00:00"
    shell:
        """
        bbduk.sh threads={threads} -Xmx10g in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} ref=adapters ktrim=r k=21 mink=11 hdist=2 tpe tbo 
        """


#bbduk will remove phiX contamination
rule bbduk_removephiX:
    input:
        "intermediates/{sample}_R1.adapterclean.fastq.gz",
        "intermediates/{sample}_R2.adapterclean.fastq.gz"
    output:
        "intermediates/{sample}_R1.phiXclean.fastq.gz",
        "intermediates/{sample}_R2.phiXclean.fastq.gz"
    threads: 16
    conda:
        "envs/bbmap.yaml"
    resources: mem_mb=10000, time="1-00:00:00"
    shell:
        """
        bbduk.sh threads={threads} -Xmx10g in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} ref=phix ktrim=r k=21 mink=11 hdist=2 minlen=50 qtrim=r trimq=15
        """

#download references for clean up

rule download_ref:
    output:
        "resources/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        "resources/GCF_000001635.27_GRCm39_genomic.fna.gz",
        "resources/GCF_015227675.2_mRatBN7.2_genomic.fna.gz"
    threads: 1
    resources: mem_mb=100000, time="1-00:00:00"
    shell:
        """
        #human
        wget -P resources/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
        #mouse
        wget -P resources/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
        #rat
        wget -P resources/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/227/675/GCF_015227675.2_mRatBN7.2/GCF_015227675.2_mRatBN7.2_genomic.fna.gz
        """

#bbmap will remove human, rat, and mouse sequences to make things easier in the end
#this is done in one step, using 3 commands. It can be adjusted to only remove certain fastas
rule remove_host:
    input:
        "resources/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
        "resources/GCF_000001635.27_GRCm39_genomic.fna.gz",
        "resources/GCF_015227675.2_mRatBN7.2_genomic.fna.gz",
        "intermediates/{sample}_R1.phiXclean.fastq.gz",
        "intermediates/{sample}_R2.phiXclean.fastq.gz"
    output:
        "intermediates/{sample}_R1.cleaned.fastq.gz",
        "intermediates/{sample}_R2.cleaned.fastq.gz",
    threads: 16
    conda:
        "envs/bbmap.yaml"
    resources: mem_mb=60000, time="1-00:00:00"
    shell:
        """
        #human
        bbmap.sh -Xmx50g threads={threads} ref=resources/GCF_000001405.40_GRCh38.p14_genomic.fna.gz nodisk minid=0.98 in1={input[0]} in2={input[1]} outu1=intermediates/{wildcards.sample}_R1.nohuman.fastq.gz outu2=intermediates/{wildcards.sample}_R2.nohuman.fastq.gz
        #mouse
        bbmap.sh -Xmx50g threads={threads} ref=resources/GCF_000001635.27_GRCm39_genomic.fna.gz nodisk minid=0.98 in1=intermediates/{wildcards.sample}_R1.nohuman.fastq.gz in2=intermediates/{wildcards.sample}_R2.nohuman.fastq.gz outu1=intermediates/{wildcards.sample}_R1.nomouse.fastq.gz outu2=intermediates/{wildcards.sample}_R2.nomouse.fastq.gz
        #rat
        bbmap.sh -Xmx50g threads={threads} ref=resources/GCF_015227675.2_mRatBN7.2_genomic.fna.gz nodisk minid=0.98 in1=intermediates/{wildcards.sample}_R1.nomouse.fastq.gz in2=intermediates/{wildcards.sample}_R2.nomouse.fastq.gz outu1={output[0]} outu2={output[1]}
        """

#fastqc after all the filtering
rule fastqc_postfilter:
    input:
        "intermediates/{sample}_R1.cleaned.fastq.gz",
        "intermediates/{sample}_R2.cleaned.fastq.gz"
    output:
        "intermediates/postfilter_qc/{sample}_R1.cleaned_fastqc.html",
        "intermediates/postfilter_qc/{sample}_R2.cleaned_fastqc.html"
    threads: 16
    conda:
        "envs/fastqc.yaml"
    resources: mem_mb=10000, time="1-00:00:00"
    shell:
        """
        if [ ! -d "intermediates/postfilter_qc" ]; then mkdir intermediates/postfilter_qc; fi
        fastqc -t {threads} -o intermediates/postfilter_qc {input}
        """

#convert fastqcs to multiqc
rule multiqc_postfilter:
    input:
        expand("intermediates/postfilter_qc/{sample}_R1.cleaned_fastqc.html", sample=NAMES),
        expand("intermediates/postfilter_qc/{sample}_R2.cleaned_fastqc.html", sample=NAMES)
    output:
        "intermediates/postfilter_qc/multiqc_report.html"
    conda: "envs/multiqc.yaml"
    threads: 1
    resources: mem_mb=10000, time="0-00:10:00"
    shell:
        """
        multiqc intermediates/postfilter_qc --outdir intermediates/postfilter_qc -f
        #rm -r intermediates/postfilter_qc/*fastqc*
        """



#bbmap's reformat.sh will merge the paired end reads into a single interleaved file for ease of processing
rule fastq_merge:
    input:
        "intermediates/{sample}_R1.cleaned.fastq.gz",
        "intermediates/{sample}_R2.cleaned.fastq.gz"
    output:
        "data/interleave_DNA/{sample}.interleave.fastq.gz"
    threads: 16
    conda:
        "envs/bbmap.yaml"
    resources: mem_mb=10000, time="1-00:00:00"
    shell:
        """
        if [ ! -d "data/interleave_DNA" ]; then mkdir data/interleave_DNA; fi
        reformat.sh threads={threads} in1={input[0]} in2={input[1]} out={output}
        """
