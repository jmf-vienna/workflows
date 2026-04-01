#this rule changes bam to fastq.gz

rule fastp:
    input:
        "intermediates/{sample}.orig.interleave.fastq.gz"
    output:
        "intermediates/{sample}.cleaned.interleave.fastq.gz",
        "intermediates/{sample}.json"
    conda:
        "../envs/fastp.yaml"
    threads: 8
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        fastp -i {input} --interleaved_in -j {output[1]} --thread {threads} -l 50 -e 20 --cut_right --stdout > intermediates/{wildcards.sample}.cleaned.interleave.fastq
        pigz -p {threads} intermediates/{wildcards.sample}.cleaned.interleave.fastq
        """
