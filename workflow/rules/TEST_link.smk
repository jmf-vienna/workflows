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

