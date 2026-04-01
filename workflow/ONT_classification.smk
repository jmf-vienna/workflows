#workflow for non-depleted RNA/transcriptomes libraries

FILES = glob_wildcards('data/NP/{name}.fastq.gz')
NAMES = FILES.name

configfile: "config/general_configs.yaml"

#This rule is just here to "request" the final results and set everything into motion
rule complete:
    input:
        expand("results/{name}/{name}.SSUrRNA.fq.gz", name=NAMES),
        expand("results/{name}/{name}.other.fq.gz", name=NAMES),
        expand("results/{name}/{name}.LSUrRNA.fq.gz", name=NAMES),
        expand("results/{name}/{name}.norRNA.fq.gz", name=NAMES),
        expand("results/{name}/{name}.SSUrRNA.kraken.report.txt", name=NAMES),
        expand("results/{name}/{name}.norRNA.diamond.tsv", name=NAMES),
        expand("results/{name}/{name}.stats.tsv", name=NAMES)


rule chopper_cleanup:
    input:
        "data/NP/{sample}.fastq.gz"
    output:
        "data/cleaned_NP/{sample}.cleaned.fastq.gz"
    threads: 8
    conda:
        "envs/chopper.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        gunzip -c {input} | chopper -q 15 -l 200 -t {threads} | gzip > {output}
        """

#rule wget_SSU:
#    output: "resources/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
#    threads: 1
#    resources: mem_mb=100000, runtime="1d", partition="basic"
#    shell:
#        """
#        wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz -P resources
#        gunzip resources/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
#        """

#rule wget_LSU:
#    output: "resources/SILVA_138.2_LSURef_NR99_tax_silva.fasta"
#    threads: 1
#    resources: mem_mb=100000, runtime="1d", partition="basic"
#    shell:
#        """
#        wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz -P resources
#        gunzip resources/SILVA_138.2_LSURef_NR99_tax_silva.fasta.gz
#        """

rule sortmerna_SSU:
    input:
        "/lisc/data/scratch/jmf/resources/SILVA/SILVA_138.2_SSURef_NR99_tax_silva.fasta",
        "data/cleaned_NP/{sample}.cleaned.fastq.gz"
    output:
        "results/{sample}/{sample}.SSUrRNA.fq.gz",
        "results/{sample}/{sample}.other.fq.gz"
    threads: 16
    conda:
        "envs/sortmerna.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        if [ -d "intermediates/{wildcards.sample}_SSU" ]; then rm -r intermediates/{wildcards.sample}_SSU; fi
        #mkdir intermediates/{wildcards.sample}_SSU
        #mkdir intermediates/{wildcards.sample}_SSU/run
        sortmerna -ref {input[0]} --reads {input[1]} --fastx true --aligned results/{wildcards.sample}/{wildcards.sample}.SSUrRNA --other results/{wildcards.sample}/{wildcards.sample}.other --threads {threads} --workdir intermediates/{wildcards.sample}_SSU/
        """

rule sortmerna_LSU:
    input:
        "/lisc/data/scratch/jmf/resources/SILVA/SILVA_138.2_LSURef_NR99_tax_silva.fasta",
        "results/{sample}/{sample}.other.fq.gz"
    output:
        "results/{sample}/{sample}.LSUrRNA.fq.gz",
        "results/{sample}/{sample}.norRNA.fq.gz"
    threads: 16
    conda:
        "envs/sortmerna.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        if [ -d "intermediates/{wildcards.sample}_LSU" ]; then rm -r intermediates/{wildcards.sample}_LSU; fi
        #mkdir intermediates/{wildcards.sample}_LSU
        #mkdir intermediates/{wildcards.sample}_LSU/run
        sortmerna -ref {input[0]} --reads {input[1]} --fastx true --aligned results/{wildcards.sample}/{wildcards.sample}.LSUrRNA --other results/{wildcards.sample}/{wildcards.sample}.norRNA --threads 16 --workdir intermediates/{wildcards.sample}_LSU/
        """

rule chopper_cleanup16S:
    input:
        "results/{sample}/{sample}.SSUrRNA.fq.gz"
    output:
        "results/{sample}/{sample}.SSUrRNA.800.fq.gz"
    threads: 8
    conda:
        "envs/chopper.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        gunzip -c {input} | chopper -l 800 -t {threads} | gzip > {output}
        """


#rule kraken_download:
#    input:
#    output:
#    threads:
#    conda:
#    resources:
#    shell:
#        """
#        wget https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz -P resources
#        tar -xvzf resources/16S_Silva138_20200326.tgz -C resources
#        """



rule kraken_id:
    input:
        "results/{sample}/{sample}.SSUrRNA.800.fq.gz",
        "/lisc/data/scratch/jmf/resources/Kraken2/16S_SILVA138_k2db"
    output:
        "results/{sample}/{sample}.SSUrRNA.kraken.out.txt",
        "results/{sample}/{sample}.SSUrRNA.kraken.report.txt"
    threads: 16
    conda:
        "envs/kraken2.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        kraken2 --db {input[1]} --threads 16 {input[0]} --output {output[0]} --gzip-compressed --report {output[1]}
        kraken2 --db {input[1]} --threads 16 {input[0]} --output {output[0]}.NAMES --gzip-compressed --report {output[1]}.NAMES --use-mpa-style --use-names
        """

rule fq2fa:
    input:
        "results/{sample}/{sample}.norRNA.fq.gz"
    output:
        "results/{sample}/{sample}.norRNA.fa.gz"
    threads: 16
    conda:
        "envs/seqkit.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        seqkit fq2fa -j {threads} {input} -o {output}
        """


rule diamond:
    input:
        "results/{sample}/{sample}.norRNA.fa.gz",
        "/lisc/data/scratch/jmf/resources/DIAMOND/nr_06020225_diamond.dmnd"
    output:
        "results/{sample}/{sample}.norRNA.diamond.tsv"
    threads: 32
    conda:
        "envs/diamond.yaml"
    resources: mem_mb=100000, runtime="3d", partition="basic"
    shell:
        """
        diamond blastx --threads {threads} --db {input[1]} --query {input[0]} --out {output} --long-reads --outfmt 102 --include-lineage
        """


rule seqkit:
    input:
        "data/cleaned_NP/{sample}.cleaned.fastq.gz",
        "results/{sample}/{sample}.SSUrRNA.fq.gz",
        "results/{sample}/{sample}.SSUrRNA.800.fq.gz",
        "results/{sample}/{sample}.LSUrRNA.fq.gz",
        "results/{sample}/{sample}.norRNA.fq.gz"
    output:
        "results/{sample}/{sample}.stats.tsv"
    threads: 8
    conda:
        "envs/seqkit.yaml"
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        seqkit stats {input[0]} {input[1]} {input[2]} {input[3]} {input[4]} -j 8 -T > {output}
        """

