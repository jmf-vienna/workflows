#mapping prep for the binning programs
rule co_mapping_prep:
    input:
        "results/coassembly/scaffolds_1000bp.fa"
    output:
        directory("results/coassembly/bams/")
    threads: 32
    conda:
        "../envs/bbmap_samtools.yaml"
    resources: mem_mb=250000, runtime="5d", partition="basic"
    log: "log/co_mapping.log"
    shell:
        """
        if [ ! -d "results/coassembly/bams" ]; then mkdir results/coassembly/bams; fi
        
        for f in data/interleave_DNA/*.interleave.fastq.gz; do g=${{f##*/}}; bbmap.sh -Xmx250g threads={threads} ref={input} interleaved=true nodisk minid=0.95 in=$f out=results/coassembly/bams/${{g%%.interleave.fastq.gz}}.bam; done 2> {log}

        #now we need to sort the bams
        for f in results/coassembly/bams/*bam; do samtools sort --write-index -@ {threads} -O BAM -o ${{f%%.bam}}.sorted.bam $f; done
        """

#Metabat with coverage
rule co_metabat_cov:
    input:
        "results/coassembly/scaffolds_1000bp.fa",
        "results/coassembly/bams"
    output:
        directory("results/coassembly/metabat_cov")
    threads: 8
    conda:
        "../envs/metabat.yaml"
    resources: mem_mb=100000, runtime="3d", partition="basic"
    log: "log/co_metabat.log"
    shell:
        """
        if [ ! -d "results/coassembly/metabat_cov" ]; then mkdir results/coassembly/metabat_cov; fi
        jgi_summarize_bam_contig_depths --outputDepth results/coassembly/metabat_cov/depth.txt results/coassembly/bams/*sorted.bam
        metabat2 -i {input[0]} -o results/coassembly/metabat_cov/coassembly_metabat_cov -m 1500 -t {threads} -a results/coassembly/metabat_cov/depth.txt 2> {log}
        """
