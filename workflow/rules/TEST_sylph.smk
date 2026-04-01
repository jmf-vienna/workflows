rule sylph:
    input:
        "intermediates/{sample}.cleaned.1.fastq.gz",
        "intermediates/{sample}.cleaned.2.fastq.gz"
    output:
        "intermediates/{sample}.cleaned.1.fastq.gz.sylphmpa"
    conda:
        "envs/sylph.yaml"
    threads: 8
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        sylph profile {config[SYLPHDB]} -o intermediates/{wildcards.sample}.sylph.tsv -u --read-seq-id 99 -t {threads} -1 {input[0]} -2 {input[1]}
        sylph-tax taxprof intermediates/{wildcards.sample}.sylph.tsv -t GlobDB_r226 -o intermediates/
        """

rule sylph_tax:
    input:
        expand("intermediates/{sample}.cleaned.1.fastq.gz.sylphmpa", sample=NAMES)
    output:
        "intermediates/sylph_abundance.tsv"
    conda:
        "envs/sylph.yaml"
    threads: 8
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
     [<35;28;19M]   """
        sylph-tax merge intermediates/*.sylphmpa --column sequence_abundance -o intermediates/sylph_abundance.tsv
        """
