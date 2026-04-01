#mapping prep for the binning programs
rule mapping_prep:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    output:
        directory("results/{sample}_spades/bams/")
    threads: 32
    conda:
        "../envs/bbmap_samtools.yaml"
    resources: mem_mb=250000, runtime="5d", partition="basic"
    log: "log/mapping_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/bams" ]; then mkdir results/{wildcards.sample}_spades/bams; fi
        
        for f in data/interleave_DNA/*.interleave.fastq.gz; do g=${{f##*/}}; bbmap.sh -Xmx250g threads={threads} ref={input} interleaved=true nodisk minid=0.95 in=$f out=results/{wildcards.sample}_spades/bams/${{g%%.interleave.fastq.gz}}.bam; done 2> {log}

        #now we need to sort the bams
        for f in results/{wildcards.sample}_spades/bams/*bam; do samtools sort --write-index -@ {threads} -O BAM -o ${{f%%.bam}}.sorted.bam $f; done
        """

#Metabat with coverage
rule metabat_cov:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa",
        "results/{sample}_spades/bams"
    output:
        directory("results/{sample}_spades/metabat_cov")
    threads: 8
    conda:
        "../envs/metabat.yaml"
    resources: mem_mb=100000, runtime="3d", partition="basic"
    log: "log/metabat_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/metabat_cov" ]; then mkdir results/{wildcards.sample}_spades/metabat_cov; fi
        jgi_summarize_bam_contig_depths --outputDepth results/{wildcards.sample}_spades/metabat_cov/depth.txt results/{wildcards.sample}_spades/bams/*sorted.bam
        metabat2 -i {input[0]} -o results/{wildcards.sample}_spades/metabat_cov/{wildcards.sample}_metabat_cov -m 1500 -t {threads} -a results/{wildcards.sample}_spades/metabat_cov/depth.txt 2> {log}
        """
