rule coverM:
    input:
        directory("results/FULL_MAGs/GTDB_env.txt")
    output:
        "results/FULL_MAGs/dereplicated_MAGs/coverM_relative_output.tsv",
        "results/FULL_MAGs/dereplicated_MAGs/coverM_mean_output.tsv"
    threads: 16
    conda:
        "../envs/coverm.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    log:
        "log/coverm_project.log"
    shell:
        """
        files=$(ls data/interleave_DNA/* | tr '\n' ' ')

        coverm genome -d results/FULL_MAGs/dereplicated_MAGs --interleaved $files -x fa --methods relative_abundance -o {output[0]} -t {threads}
        coverm genome -d results/FULL_MAGs/dereplicated_MAGs --interleaved $files -x fa --methods mean -o {output[1]} -t {threads}
        """
