#Request all necessary outputs

configfile: "config/general_configs.yaml"

rule all:
    input:
        "results/sylph/sylph_merged.tsv"

#sylph
rule sylph:
    output:
        "results/sylph/sylph_merged.tsv"
    conda: "envs/sylph.yaml"
    threads: 8
    resources: mem_mb=100000, runtime="1d", partition="basic"
    shell:
        """
        sylph profile {config[SYLPHDB]} -u --read-seq-id 99 -t {threads} -1 intermediates/test/*R1.phiXclean.fastq.gz -2 intermediates/test/*R2.phiXclean.fastq.gz -o results/sylph/out.tsv


        #sylph profile /lisc/scratch/jmf/resources/sylph/gtdb-r220-c200-dbv1.syldb -u --read-seq-id 99 -t {threads} -1 intermediates/test/*R1.phiXclean.fastq.gz -2 intermediates/test/*R2.phiXclean.fastq.gz -o results/sylph/out.tsv

        # incorporate GTDB-r220 taxonomy into sylph's results
        sylph-tax download --download-to resources/
        sylph-tax taxprof results/sylph/out.tsv -t GTDB_r220 -o results/sylph/

        sylph-tax merge results/sylph/*sylphmpa -o results/sylph/sylph_merged.tsv --column sequence_abundance
        rm results/sylph/*sylphmpa
        rm results/sylph/out.tsv
        
        """
