#binny
rule binny:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa",
        "results/{sample}_spades/bams"
    output:
        directory("results/{sample}_spades/binny")
    threads: 32
    conda:
        "/lisc/data/scratch/jmf/conda_envs/binny"
    resources: mem_mb=250000, runtime="5d", partition="basic"
    log: "log/binny_{sample}.log"
    shell:
        """
        binny --assembly {input[0]} --bam {input[1]}/*sorted.bam --outputdir output --tmp_dir $TMPDIR --threads {threads} --sample {wildcards.sample}_binny
        #should add something to fix the naming on these. It sucks.
        """
