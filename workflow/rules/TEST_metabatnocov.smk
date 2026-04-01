
#binning without coverage input
rule metabat_nocov:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    output:
        directory("results/{sample}_spades/metabat_nocov/")
    threads: 8
    conda:
        "../envs/metabat.yaml"
    resources: mem_mb=50000, runtime="4d", partition="basic"
    log: "log/metabatnocov_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/metabat_nocov" ]; then mkdir results/{wildcards.sample}_spades/metabat_nocov; fi
        metabat2 -m 1500 -t {threads} -i {input} -o results/{wildcards.sample}_spades/metabat_nocov/{wildcards.sample}_metabat_nocov 2> {log}
        """
