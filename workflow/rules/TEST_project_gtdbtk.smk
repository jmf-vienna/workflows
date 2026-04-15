rule gtdb_binning_setup:
    input:
        directory("results/FULL_MAGs/dereplicated_MAGs")
    output:
        "results/FULL_MAGs/GTDB_env.txt"
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
        "results/FULL_MAGs/GTDB_env.txt"
    output:
        "results/FULL_MAGs/dereplicated_MAGs/gtdb_classification/gtdbtk.bac120.summary.tsv"
    threads: 16
    conda:
        "../envs/gtdbtk.yaml"
    resources: mem_mb=200000, runtime="1d", partition="basic"
    log: "log/gtdbtk_project.log"
    shell:
        """
        gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/FULL_MAGs/dereplicated_MAGs --out_dir results/FULL_MAGs/dereplicated_MAGs/gtdb_classification 2> {log}
        """
