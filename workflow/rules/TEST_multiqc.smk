#this rule generates a multiqc report

rule multiqc:
    input:
        expand("intermediates/{sample}.json", sample=NAMES)
    output:
        "intermediates/multiqc_report.html"
    conda:
        "../envs/multiqc.yaml"
    threads: 8
    resources: mem_mb=10000, runtime="1d", partition="basic"
    shell:
        """
        multiqc intermediates --force --outdir intermediates --clean-up
        """
