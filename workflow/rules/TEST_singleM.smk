#globDB singleM
rule singleM:
    input:
        "intermediates/{sample}.cleaned.1.fastq.gz",
        "intermediates/{sample}.cleaned.2.fastq.gz"
    output:
        "intermediates/{sample}.spf.tsv"
    conda:
        "envs/singleM.yaml"
    threads: 16
    resources: mem_mb=100000, runtime="1d", partition="basic"
        shell:
        """
        singlem pipe --threads 16 -1 {input[0]} -2 {input[1]} --metapackage {config[SINGLEMDB]} --taxonomic-profile intermediates/{wildcards.sample}.profile
        touch {output}

        singlem microbial_fraction --metapackage {config[SINGLEMDB]} --forward {input[0]} --reverse {input[1]} -p intermediates/{wildcards.sample}.profile > {output}
        """

#singleM merge
rule singleM_merge:
    input:
        expand("intermediates/{sample}.spf.tsv", sample=NAMES)
    output:
        "intermediates/singleM_microbialfrac.tsv"
    threads: 1
    resources: mem_mb=10000, runtime="1h", partition="basic"
    shell:
        """
        for f in intermediates/*spf.tsv; do head -1 $f > intermediates/singleM_microbialfrac.tsv; done
        for f in intermediates/*spf.tsv; do tail -1 $f >> {output} ; done
        """
