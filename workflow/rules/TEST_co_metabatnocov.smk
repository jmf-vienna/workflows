
#binning without coverage input
rule co_metabat_nocov:
    input:
        "results/coassembly/scaffolds_1000bp.fa"
    output:
        directory("results/coassembly/metabat_nocov/")
    threads: 8
    conda:
        "../envs/metabat.yaml"
    resources: mem_mb=50000, runtime="4d", partition="basic"
    log: "log/co_metabatnocov.log"
    shell:
        """
        if [ ! -d "results/coassembly/metabat_nocov" ]; then mkdir results/coassembly/metabat_nocov; fi
        metabat2 -m 1500 -t {threads} -i {input} -o results/coassembly/metabat_nocov/coassembly_metabat_nocov 2> {log}
        """
