#binny
rule binny:
    input:
        "results/coassembly/scaffolds_1000bp.fa",
        "results/coassembly/bams"
    output:
        directory("results/coassembly/binny")
    threads: 32
    conda:
        "/lisc/data/scratch/jmf/conda_envs/binny"
    resources: mem_mb=250000, runtime="5d", partition="basic"
    log: "log/binny_co.log"
    shell:
        """
        
        binny --assembly {input[0]} --bam {input[1]}/*sorted.bam --outputdir output --tmp_dir $TMPDIR --threads {threads} --sample coassembly_binny
        #should add something to fix the naming on these. It sucks.
        """
