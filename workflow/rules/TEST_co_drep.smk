#A quick dREP
rule co_drep:
    input:
        "results/coassembly/metabat_nocov/",
        "results/coassembly/metabat_cov/"

    output:
        "results/coassembly/final_bins/data_tables/genomeInformation.csv"
    threads: 16
    envmodules:
        "dRep"
    #conda:
    #    "../envs/drep.yaml" # drep not working
    resources: mem_mb=300000, runtime="4d", partition="basic"
    log: "log/co_drep.log"
    shell:
        """
        if [ ! -d "results/coassembly/final_bins" ]; then mkdir results/coassembly/final_bins; fi
        cp results/coassembly/metabat_cov/coassembly_metabat_cov*fa results/coassembly/final_bins
        cp results/coassembly/metabat_nocov/coassembly_metabat_nocov*fa results/coassembly/final_bins

        dRep dereplicate results/coassembly/final_bins/ -comp 50 -p {threads} -g results/coassembly/final_bins/*fa 2> {log}
        """
