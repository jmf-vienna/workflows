#this rule changes bam to fastq.gz

rule bamtofq:
    input:
        "intermediates/{sample}.bam"
    output:
        "intermediates/{sample}.orig.interleave.fastq.gz"
    conda:
       "../envs/samtools.yaml"
    threads: 8
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        samtools fastq --threads {threads} {input} > intermediates/{wildcards.sample}.orig.interleave.fastq
        pigz -p {threads} intermediates/{wildcards.sample}.orig.interleave.fastq
        """
