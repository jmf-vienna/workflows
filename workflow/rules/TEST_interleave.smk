#this rule changes bam to fastq.gz

rule interleave:
    input:
        "intermediates/{sample}.1.fastq.gz",
        "intermediates/{sample}.2.fastq.gz"
    output:
        "intermediates/{sample}.orig.interleave.fastq.gz"
    conda:
       "../envs/bbmap.yaml"
    threads: 8
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        reformat.sh in1={input[0]} in2={input[1]} threads={threads} out={output}
        """
