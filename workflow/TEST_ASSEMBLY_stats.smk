#workflow for spades assembler and binning process afterwards.

FILES = glob_wildcards('results/{name}_spades/scaffolds_1000bp.fa')
NAMES = FILES.name

configfile: "config/general_configs.yaml"

#This rule is just here to "request" the final results and set everything into motion
rule complete:
    input:
        "results/UNMAPPED_merged.tsv"


rule coverM_assembly:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    output:
        "results/{sample}_spades/coverM_relative_output.tsv"
    threads: 16
    conda:
        "envs/coverm.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        coverm genome -f {input[0]} --interleaved data/interleave_DNA/{wildcards.sample}.interleave.fastq.gz --methods relative_abundance -o {output[0]} -t 16

        sed -i 's/unmapped/{wildcards.sample}_unmapped_assembly/' {output[0]}
        """


rule coverM_genomes:
    input:
        directory("results/{sample}_spades/final_bins/dereplicated_genomes/")
    output:
        "results/{sample}_spades/coverMgenomes_relative_output.tsv"
    threads: 16
    conda:
        "envs/coverm.yaml"
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        coverm genome -d results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/ --interleaved data/interleave_DNA/{wildcards.sample}.interleave.fastq.gz -x fa --methods relative_abundance -o {output[0]} -t 16

        sed -i 's/unmapped/{wildcards.sample}_unmapped_genomes/' {output[0]}
        """


rule merge:
    input:
        expand("results/{sample}_spades/coverM_relative_output.tsv", sample=NAMES),
        expand("results/{sample}_spades/coverMgenomes_relative_output.tsv", sample=NAMES)
    output:
        "results/UNMAPPED_merged.tsv"
    threads: 1
    resources: mem_mb=10000, runtime="1h", partition="basic"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in results/*/coverM_relative_output.tsv; do grep "unmapped" $f >> {output}; rm $f; done
        for f in results/*/coverMgenomes_relative_output.tsv; do grep "unmapped" $f >> {output}; rm $f; done
        """
