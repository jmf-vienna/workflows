#nonpareil
rule nonpareil:
    input:
        "intermediates/{sample}.cleaned.1.fastq.gz"
    output:
        "intermediates/{sample}.npo"
    conda:
        "envs/nonpareil.yaml"
    threads: 16
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        nonpareil -t 16 -s {input} -T kmer -f fastq -b intermediates/{wildcards.sample}
        """


#merge and generate things
rule nonpareil_stats_gen:
    input:
        expand("intermediates/{sample}.npo", sample=NAMES)
    output:
        "intermediates/samples.txt"
    conda:
        "envs/R.yaml"
    threads: 1
    resources: mem_mb=50000, runtime="4h", partition="basic"
    shell:
        """
        echo "FileTABName" > intermediates/samples.txt
        for f in intermediates/*npo; do g=${{f##*/}}; echo "${{f}}TAB${{g%%.npo}}" >> intermediates/samples.txt; done
        sed -i 's/TAB/\t/g' intermediates/samples.txt

        Rscript workflow/scripts/Nonpareil.R
        """
