rule gtdb_binning_setup:
    input:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv"
    output:
        "results/{sample}_spades/final_bins/GTDB_env.txt"
    threads: 1
    conda:
        "../envs/gtdbtk.yaml"
    resources: mem_mb=1000, runtime="1d", partition="basic"
    shell:
        """
        conda env config vars set GTDBTK_DATA_PATH={config[GTDBPATH]}
        touch {output}
        """

        #at this point lets aim for gtdb and add to the dREP
rule gtdb_binning:
    input:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv",
        "results/{sample}_spades/final_bins/GTDB_env.txt"
    output:
        "results/{sample}_spades/final_bins/dereplicated_genomes/gtdb_classification/gtdbtk.bac120.summary.tsv"
    threads: 16
    conda:
        "../envs/gtdbtk.yaml"
    resources: mem_mb=200000, runtime="1d", partition="basic"
    log: "log/gtdbtk_{sample}.log"
    shell:
        """
        gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes --out_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification --skip_ani_screen 2> {log}
        """
