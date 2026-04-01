#workflow for spades assembler and binning process afterwards.

FILES = glob_wildcards('data/interleave_DNA/{name}.interleave.fastq.gz')
NAMES = FILES.name

configfile: "config/general_configs.yaml"

#This rule is just here to "request" the final results and set everything into motion
rule complete:
    input:
        expand("results/{name}_spades/final_bins/dereplicated_genomes/final_QC.xlsx", name=NAMES)


#Assemble with SPAdes using interleaved reads, will only keep >1000bp contigs
rule SPAdes_assembly:
    input:
        "data/interleave_DNA/{sample}.interleave.fastq.gz"
    output:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    threads: 48
    conda:
        "envs/spades.yaml"
    resources: mem_mb=500000, runtime="1d", partition="basic"
    log: "log/SPADES_{sample}.log"
    shell:
        """
        spades.py -t {threads} -m 500 -k 21,31,41,51,61,71,81,91,101,111,121 --meta --pe-12 1 {input} -o results/{wildcards.sample}_spades 2> {log}
        reformat.sh in=results/{wildcards.sample}_spades/scaffolds.fasta out={output} minlength=1000
        """

#BINNING
#binning without coverage input
rule metabat_nocov:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    output:
        directory("results/{sample}_spades/metabat_nocov/")
    threads: 8
    conda:
        "envs/metabat.yaml"
    resources: mem_mb=50000, runtime="1d", partition="basic"
    log: "log/metabatnocov_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/metabat_nocov" ]; then mkdir results/{wildcards.sample}_spades/metabat_nocov; fi
        metabat2 -m 1500 -t {threads} -i {input} -o results/{wildcards.sample}_spades/metabat_nocov/{wildcards.sample}_metabat_nocov 2> {log}
        """

#mapping prep for the binning programs
rule mapping_prep:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa"
    output:
        directory("results/{sample}_spades/bams/")
    threads: 16
    conda:
        "envs/bbmap_samtools.yaml"
    resources: mem_mb=100000, runtime="3d", partition="basic"
    log: "log/mapping_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/bams" ]; then mkdir results/{wildcards.sample}_spades/bams; fi

        for f in data/interleave_DNA/*.interleave.fastq.gz; do g=${{f##*/}}; bbmap.sh -Xmx50g threads={threads} ref={input} interleaved=true nodisk minid=0.98 in=$f out=results/{wildcards.sample}_spades/bams/${{g%%.interleave.fastq.gz}}.bam; done 2> {log}

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
        "envs/metabat.yaml"
    resources: mem_mb=50000, runtime="1d", partition="basic"
    log: "log/metabat_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/metabat_cov" ]; then mkdir results/{wildcards.sample}_spades/metabat_cov; fi
        jgi_summarize_bam_contig_depths --outputDepth results/{wildcards.sample}_spades/metabat_cov/depth.txt results/{wildcards.sample}_spades/bams/*sorted.bam
        metabat2 -i {input[0]} -o results/{wildcards.sample}_spades/metabat_cov/{wildcards.sample}_metabat_cov -m 1500 -t {threads} -a results/{wildcards.sample}_spades/metabat_cov/depth.txt 2> {log}
        """


#A quick dREP
rule drep:
    input:
        "results/{sample}_spades/metabat_nocov/",
        "results/{sample}_spades/metabat_cov/"

    output:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv"
    threads: 8
#    conda:
#        "envs/drep.yaml"
    resources: mem_mb=200000, runtime="1d", partition="basic"
    log: "log/drep_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/final_bins" ]; then mkdir results/{wildcards.sample}_spades/final_bins; fi
        cp results/{wildcards.sample}_spades/metabat_cov/{wildcards.sample}_metabat_cov*fa results/{wildcards.sample}_spades/final_bins
        cp results/{wildcards.sample}_spades/metabat_nocov/{wildcards.sample}_metabat_nocov*fa results/{wildcards.sample}_spades/final_bins

        module load conda
        conda activate drep-3.5.0

        dRep dereplicate results/{wildcards.sample}_spades/final_bins/ -p {threads} -g results/{wildcards.sample}_spades/final_bins/*fa 2> {log}
        """

rule gtdb_binning_setup:
    input:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv"
    output:
        "results/{sample}_spades/final_bins/GTDB_env.txt"
    threads: 1
    conda:
        "envs/gtdbtk.yaml"
    resources: mem_mb=1000, runtime="1d", partition="basic"
    shell:
        """
        conda env config vars set GTDBTK_DATA_PATH={config[GTDBPATH]}
        touch {output}
        """
#stats for each MAG
rule MAG_stats:
    input:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv"
    output:
        "results/{sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv"
    threads: 1
    conda:
        "envs/bbmap_samtools.yaml"
    resources: mem_mb=10000, runtime="1d", partition="basic"
    log: "log/stats_{sample}.log"
    shell:
        """
        #generate stats file for each MAG
        for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*fa;
        do g=${{f##*/}};
        stats.sh in=$f format=3 > ${{f%%.fa}}.stats.tsv;
        echo 'genome' > ${{f%%.fa}}.name.txt;
        echo "$g" >> ${{f%%.fa}}.name.txt;
        paste -d'\t' ${{f%%.fa}}.name.txt ${{f%%.fa}}.stats.tsv > ${{f%%.fa}}.stats.merged.tsv;
        rm ${{f%%.fa}}.stats.tsv;
        rm ${{f%%.fa}}.name.txt;
        done

        #get headers for stats
        for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*stats.merged.tsv; do head -1 $f > results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv ; done

        #merge all stats into 1 file
        for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*stats.merged.tsv;
        do tail -1 $f >> results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv;
        rm $f;
        done
        """

#at this point lets aim for gtdb and add to the dREP
rule gtdb_binning:
    input:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv",
        "results/{sample}_spades/final_bins/GTDB_env.txt"
    output:
        "results/{sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv"
    threads: 16
    conda:
        "envs/gtdbtk.yaml"
    resources: mem_mb=200000, runtime="1d", partition="basic"
    log: "log/gtdbtk_{sample}.log"
    shell:
        """
        gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes --out_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification --skip_ani_screen 2> {log}
        mv results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification/classify/gtdbtk.bac120.summary.tsv results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv
        """


rule drep_stats_gtbd_merge:
    input:
        "results/{sample}_spades/final_bins/data_tables/genomeInformation.csv",
        "results/{sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv",
        "results/{sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv"
    output:
        "results/{sample}_spades/final_bins/dereplicated_genomes/final_QC.xlsx"
    conda:
        "envs/python3_modules.yaml"
    threads: 1
    resources: mem_mb=10000, runtime="1h", partition="basic"
    shell:
        """
        python3 workflow/scripts/derep_gtdb_stats_merge.py {input[0]} {input[1]} results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification/gtdbtk.ar53.summary.tsv {input[2]} {output}
        """




#dREP all
rule drep_all:
    input:
    expand("results/{name}_spades/metabat_nocov/", name=NAMES),
    expand("results/{name}_spades/metabat_cov/", name=NAMES)

    output:
        "results/final_bins/data_tables/genomeInformation.csv"
    threads: 16
#    conda:
#        "envs/drep.yaml"
    resources: mem_mb=200000, runtime="2d", partition="basic"
    log: "log/drepall_all.log"
    shell:
        """
        if [ ! -d "results/final_bins" ]; then mkdir results/final_bins; fi
        for f in results/*_spades/metabat_cov/*metabat_cov*fa; ln -s $f results/final_bins/${f##*/} ; done
        for f in results/*_spades/metabat_nocov/*metabat_nocov*fa; ln -s $f results/final_bins/${f##*/} ; done

        module load conda
        conda activate drep-3.5.0

        dRep dereplicate results/final_bins/ -p {threads} -g results/final_bins/*fa 2> {log}
        """


#stats for each MAG
rule MAG_stats_all:
    input:
	"results/final_bins/data_tables/genomeInformation.csv"
    output:
	"results/final_bins/dereplicated_genomes/complete_stats.tsv"
    threads: 1
    conda:
	"envs/bbmap_samtools.yaml"
    resources: mem_mb=10000, runtime="1d", partition="basic"
    log: "log/stats_all.log"
    shell:
	"""
        #generate stats file for each MAG
        for f in results/final_bins/dereplicated_genomes/*fa;
        do g=${{f##*/}};
        stats.sh in=$f format=3 > ${{f%%.fa}}.stats.tsv;
        echo 'genome' > ${{f%%.fa}}.name.txt;
        echo "$g" >> ${{f%%.fa}}.name.txt;
        paste -d'\t' ${{f%%.fa}}.name.txt ${{f%%.fa}}.stats.tsv > ${{f%%.fa}}.stats.merged.tsv;
        rm ${{f%%.fa}}.stats.tsv;
        rm ${{f%%.fa}}.name.txt;
        done

        #get headers for stats
        for f in results/final_bins/dereplicated_genomes/*stats.merged.tsv; do head -1 $f > results/final_bins/dereplicated_genomes/complete_stats.tsv ; done

        #merge all stats into 1 file
        for f in results/final_bins/dereplicated_genomes/*stats.merged.tsv;
        do tail -1 $f >> results/final_bins/dereplicated_genomes/complete_stats.tsv;
        rm $f;
        done
        """

#at this point lets aim for gtdb and add to the dREP
rule gtdb_binning_all:
    input:
        "results/final_bins/data_tables/genomeInformation.csv",
        expand("results/{name}_spades/final_bins/GTDB_env.txt", name=NAMES)
    output:
        "results/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv"
    threads: 16
    conda:
        "envs/gtdbtk.yaml"
    resources: mem_mb=200000, runtime="1d", partition="basic"
    log: "log/gtdbtkall.log"
    shell:
        """
        gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/final_bins/dereplicated_genomes --out_dir results/final_bins/dereplicated_genomes/gtdb_classification --skip_ani_screen 2> {log}
        mv results/final_bins/dereplicated_genomes/gtdb_classification/classify/gtdbtk.bac120.summary.tsv results/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv
        """


rule drep_stats_gtbd_merge_all:
    input:
        "results/final_bins/data_tables/genomeInformation.csv",
        "results/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv",
        "results/final_bins/dereplicated_genomes/complete_stats.tsv"
    output:
	"results/final_bins/dereplicated_genomes/final_QC.xlsx"
    conda:
	"envs/python3_modules.yaml"
    threads: 1
    resources: mem_mb=10000, runtime="1h", partition="basic"
    shell:
        """
        python3 workflow/scripts/derep_gtdb_stats_merge.py {input[0]} {input[1]} results/final_bins/dereplicated_genomes/gtdb_classification/gtdbtk.ar53.summary.tsv {input[2]} {output}
        """
