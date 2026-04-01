        
#Assemble with SPAdes using interleaved reads, will only keep >1000bp contigs

rule SPAdes_assembly:
    input:
        "data/interleave_DNA/{sample}.cleaned.interleave.fastq.gz"
    output:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    threads: 64
    conda:
        "../envs/spades.yaml"
    resources: mem_mb=750000, runtime="10d", partition="basic"
    log: "log/SPADES_{sample}.log"
    shell:
        """
        spades.py -t {threads} -m 750 -k 21,31,41,51,61,71,81,91,101,111,121 --meta --tmp-dir $TMPDIR --pe-12 1 {input} -o results/{wildcards.sample}_spades 2> {log}
        reformat.sh in=results/{wildcards.sample}_spades/scaffolds.fasta out={output} minlength=1000
        """
