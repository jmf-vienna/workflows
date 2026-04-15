#move all binning sets to final
rule mv_bins:
    input:
        "results/{sample}_spades/metabat_nocov/",
        "results/{sample}_spades/metabat_cov/",
        "results/{sample}_spades/binny/"
        
    output:
        directory("results/{sample}_spades/final_bins")
    threads: 1
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/final_bins" ]; then mkdir results/{wildcards.sample}_spades/final_bins; fi
        cp results/{wildcards.sample}_spades/metabat_cov/*fa results/{wildcards.sample}_spades/final_bins
        cp results/{wildcards.sample}_spades/metabat_nocov/*fa results/{wildcards.sample}_spades/final_bins
        cp results/{wildcards.sample}_spades/binny/bins/*fa results/{wildcards.sample}_spades/final_bins
        """

#run checkm2
rule checkm2:
    input:
        "results/{sample}_spades/final_bins"
    output:
        "results/{sample}_spades/final_bins/checkm2_results/quality_report.tsv"
    threads: 16
    envmodules: "CheckM2"
    resources: mem_mb=100000, runtime="2d", partition="basic"
    log: "log/checkm2_{sample}.log"
    shell:
        """
        checkm2 predict -x fa -t {threads} -i {input} -o results/{wildcards.sample}_spades/final_bins/checkm2_results
        """

#run galah
rule galah:
    input:
        "results/{sample}_spades/final_bins/checkm2_results/quality_report.tsv"
    output:
        "results/{sample}_spades/final_bins/dereplicated_MAGs"
    threads: 16
    envmodules: "Galah"
    resources: mem_mb=100000, runtime="2d", partition="basic"
    log: "log/checkm2_{sample}.log"
    shell:
        """
        galah cluster --genome-fasta-directory results/{sample}_spades/final_bins -x fa --output-cluster-definition results/{sample}_spades/final_bins/galah_clusters.tsv --output-representative-fasta-directory {output} -t {threads} --checkm2-quality-report {input} --min-completeness 75 --max-contamination 25
        """