#this rule changes bam to fastq.gz

rule mvdna:
    input:
        "intermediates/{sample}.cleaned.interleave.fastq.gz"
    output:
        "data/interleave_DNA/{sample}.cleaned.interleave.fastq.gz"
    threads: 1
    resources: mem_mb=10000, runtime="1h", partition="basic"
    shell:
        """
        mv {input} {output}

        """