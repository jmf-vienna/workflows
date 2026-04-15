#move all binning sets to final
rule mv_bins:
    input:
        "results/coassembly/metabat_nocov/",
        "results/coassembly/metabat_cov/",
        "results/coassembly/binny/"
        
    output:
        directory("results/coassembly/final_bins")
    threads: 1
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        if [ ! -d "results/coassembly/final_bins" ]; then mkdir results/coassembly/final_bins; fi
        cp results/coassembly/metabat_cov/*fa results/coassembly/final_bins
        cp results/coassembly/metabat_nocov/*fa results/coassembly/final_bins
        cp results/coassembly/binny/bins/*fa results/coassembly/final_bins
        """

#run checkm2
rule checkm2:
    input:
        "results/coassembly/final_bins"
    output:
        "results/coassembly/final_bins/checkm2_results/quality_report.tsv"
    threads: 16
    envmodules: "CheckM2"
    resources: mem_mb=100000, runtime="2d", partition="basic"
    log: "log/checkm2_{sample}.log"
    shell:
        """
        checkm2 predict -x fa -t {threads} -i {input} -o results/coassembly/final_bins/checkm2_results
        """

#run galah
rule galah:
    input:
        "results/coassembly/final_bins/checkm2_results/quality_report.tsv"
    output:
        "results/coassembly/final_bins/dereplicated_MAGs"
    threads: 16
    envmodules: "Galah"
    resources: mem_mb=100000, runtime="2d", partition="basic"
    log: "log/checkm2_{sample}.log"
    shell:
        """
        galah cluster --genome-fasta-directory results/coassembly/final_bins -x fa --output-cluster-definition results/coassembly/final_bins/galah_clusters.tsv --output-representative-fasta-directory {output} -t {threads} --checkm2-quality-report {input} --min-completeness 75 --max-contamination 25
        """