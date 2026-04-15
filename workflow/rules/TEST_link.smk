#this rule links bams to intermediates

rule bam_link:
    input:
        "data/DNA/{sample}.bam"
    output:
        "intermediates/{sample}.bam"
    threads: 48
    resources: mem_mb=750000, runtime="1d", partition="basic"
    shell:
        """
        ln -s {input} -t intermediates -r
        """

rule split_link:
    input:
        "data/DNA/{sample}.1.fastq.gz",
        "data/DNA/{sample}.2.fastq.gz"
    output:
        "intermediates/{sample}.1.fastq.gz",
        "intermediates/{sample}.2.fastq.gz"
    threads: 48
    resources: mem_mb=750000, runtime="1d", partition="basic"
    shell:
        """
        ln -s {input[0]} -t intermediates -r
        ln -s {input[1]} -t intermediates -r
        """