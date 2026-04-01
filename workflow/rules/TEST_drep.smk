#A quick dREP
rule drep:
    input:
        "results/{sample}_spades/metabat_nocov/",
        "results/{sample}_spades/metabat_cov/"

    output:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv"
    threads: 16
    envmodules:
        "dRep"
    #conda:
    #    "../envs/drep.yaml" # drep not working
    resources: mem_mb=300000, runtime="4d", partition="basic"
    log: "log/drep_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/final_bins" ]; then mkdir results/{wildcards.sample}_spades/final_bins; fi
        cp results/{wildcards.sample}_spades/metabat_cov/{wildcards.sample}_metabat_cov*fa results/{wildcards.sample}_spades/final_bins
        cp results/{wildcards.sample}_spades/metabat_nocov/{wildcards.sample}_metabat_nocov*fa results/{wildcards.sample}_spades/final_bins

        dRep dereplicate results/{wildcards.sample}_spades/final_bins/ -comp 50 -p {threads} -g results/{wildcards.sample}_spades/final_bins/*fa 2> {log}
        """
