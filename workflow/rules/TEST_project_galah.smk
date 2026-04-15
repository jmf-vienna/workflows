#move all binning sets to final
rule mv_bins_all:
    input:
        expand("results/{sample}_spades/final_bins/dereplicated_MAGs", sample=NAMES),
        "results/coassembly/final_bins/dereplicated_MAGs"
        
    output:
        directory("results/FULL_MAGs/")
    threads: 1
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        if [ ! -d "results/FULL_MAGs" ]; then mkdir "results/FULL_MAGs"; fi
        cp results/*/final_bins/dereplicated_MAGs results/FULL_MAGs/.
        """

#run checkm2
rule checkm2:
    input:
        "results/FULL_MAGs/"
    output:
        "results/FULL_MAGs/checkm2_results/quality_report.tsv"
    threads: 16
    envmodules: "CheckM2"
    resources: mem_mb=100000, runtime="2d", partition="basic"
    log: "log/checkm2_{sample}.log"
    shell:
        """
        checkm2 predict -x fa -t {threads} -i {input} -o results/FULL_MAGs/checkm2_results
        """

#run galah
rule galah:
    input:
        "results/FULL_MAGs/checkm2_results/quality_report.tsv"
    output:
        "results/FULL_MAGs/dereplicated_MAGs"
    threads: 16
    envmodules: "Galah"
    resources: mem_mb=100000, runtime="2d", partition="basic"
    log: "log/checkm2_{sample}.log"
    shell:
        """
        galah cluster --genome-fasta-directory results/FULL_MAGs/final_bins -x fa --output-cluster-definition results/FULL_MAGs/galah_clusters.tsv --output-representative-fasta-directory {output} -t {threads} --checkm2-quality-report {input} --min-completeness 75 --max-contamination 25
        """