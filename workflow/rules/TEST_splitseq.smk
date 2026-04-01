#this rule changes finished interleave fastq.gz to split

rule splitreads:
    input:
        "data/interleave_DNA/{sample}.cleaned.interleave.fastq.gz"
    output:
        "intermediates/{sample}.cleaned.1.fastq.gz",
        "intermediates/{sample}.cleaned.1.fastq.g"
    conda:
       "envs/bbmap.yaml"
    threads: 8
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        reformat.sh threads={threads} in={input} out1={output[0]} out2={output[1]}
        """