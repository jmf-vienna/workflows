#workflow to rebin samples using selected readset (done after a larger binning using all samples)
#must delete undesired mapping data for each sample from "results/{sample}_spades/bams"

#make sure that the original binnings were as successful as desired as this will remove sample dataset for that


FILES = glob_wildcards('data/interleave_DNA/{name}.interleave.fastq.gz')
NAMES = FILES.name

configfile: "config/general_configs.yaml"

#This rule is just here to "request" the final results and set everything into motion
rule complete:
    input:
        expand("results/{name}_spades/final_bins_rebin/dereplicated_genomes/final_rebin_QC.xlsx", name=NAMES)


#Metabat with coverage
rule metabat_cov:
    input:
        "results/{sample}_spades/scaffolds_1000bp.fa",
        "results/{sample}_spades/bams_rebin"
    output:
        directory("results/{sample}_spades/metabat_cov_rebin")
    threads: 8
    conda:
        "envs/metabat.yaml"
    resources: mem_mb=50000, time="1-00:00:00"
    log: "log/metabat_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/metabat_cov_rebin" ]; then mkdir results/{wildcards.sample}_spades/metabat_cov_rebin; fi
        jgi_summarize_bam_contig_depths --outputDepth results/{wildcards.sample}_spades/metabat_cov_rebin/depth.txt results/{wildcards.sample}_spades/bams_rebin/*sorted.bam
        metabat2 -i {input[0]} -o results/{wildcards.sample}_spades/metabat_cov_rebin/{wildcards.sample}_metabat_cov -m 1500 -t {threads} -a results/{wildcards.sample}_spades/metabat_cov_rebin/depth.txt 2> {log}
        """


#A quick dREP
rule drep:
    input:
        "results/{sample}_spades/metabat_nocov/",
        "results/{sample}_spades/metabat_cov_rebin/"

    output:
        "results/{sample}_spades/final_bins_rebin/data_tables/genomeInformation.csv"
    threads: 8
    conda:
        "envs/drep.yaml"
    resources: mem_mb=200000, time="1-00:00:00"
    log: "log/drep_{sample}.log"
    shell:
        """
        if [ ! -d "results/{wildcards.sample}_spades/final_bins_rebin" ]; then mkdir results/{wildcards.sample}_spades/final_bins_rebin; fi
        cp results/{wildcards.sample}_spades/metabat_cov_rebin/{wildcards.sample}_metabat_cov*fa results/{wildcards.sample}_spades/final_bins_rebin
        cp results/{wildcards.sample}_spades/metabat_nocov/{wildcards.sample}_metabat_nocov*fa results/{wildcards.sample}_spades/final_bins_rebin

        dRep dereplicate results/{wildcards.sample}_spades/final_bins_rebin/ -p {threads} -g results/{wildcards.sample}_spades/final_bins_rebin/*fa 2> {log}
        """

#stats for each MAG
rule MAG_stats:
    input:
        "results/{sample}_spades/final_bins_rebin/data_tables/genomeInformation.csv"
    output:
        "results/{sample}_spades/final_bins_rebin/dereplicated_genomes/complete_stats.tsv"
    threads: 1
    conda:
        "envs/bbmap_samtools.yaml"
    resources: mem_mb=10000, time="0-01:00:00"
    log: "log/stats_{sample}.log"
    shell:
        """
        #generate stats file for each MAG
        for f in results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/*fa;
        do g=${{f##*/}};
        stats.sh in=$f format=3 > ${{f%%.fa}}.stats.tsv;
        echo 'genome' > ${{f%%.fa}}.name.txt;
        echo "$g" >> ${{f%%.fa}}.name.txt;
        paste -d'\t' ${{f%%.fa}}.name.txt ${{f%%.fa}}.stats.tsv > ${{f%%.fa}}.stats.merged.tsv;
        rm ${{f%%.fa}}.stats.tsv;
        rm ${{f%%.fa}}.name.txt;
        done

        #get headers for stats
        for f in results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/*stats.merged.tsv; do head -1 $f > results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/complete_stats.tsv ; done

        #merge all stats into 1 file
        for f in results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/*stats.merged.tsv;
        do tail -1 $f >> results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/complete_stats.tsv;
        rm $f;
        done
        """

rule gtdb_binning_setup:
    input:
        "results/{sample}_spades/final_bins_rebin/data_tables/genomeInformation.csv"
    output:
        "results/{sample}_spades/final_bins_rebin/GTDB_env.txt"
    threads: 1
    conda:
        "envs/gtdbtk.yaml"
    resources: mem_mb=1000, time="1-00:00:00"
    shell:
        """
        conda env config vars set GTDBTK_DATA_PATH={config[GTDBPATH]}
        touch {output}
        """

#reanalyze things
#at this point lets aim for gtdb and add to the dREP
rule gtdb_binning:
    input:
        "results/{sample}_spades/final_bins_rebin/data_tables/genomeInformation.csv",
        "results/{sample}_spades/final_bins_rebin/GTDB_env.txt"
    output:
        "results/{sample}_spades/final_bins_rebin/dereplicated_genomes/gtdbtk.bac120.summary.tsv"
    threads: 16
    conda:
        "envs/gtdbtk.yaml"
    resources: mem_mb=200000, time="2-00:00:00"
    log: "log/gtdbtk_{sample}.log"
    shell:
        """
        gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes --out_dir results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/gtdb_classification --skip_ani_screen 2> {log}
        mv results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/gtdb_classification/classify/gtdbtk.bac120.summary.tsv results/{wildcards.sample}_spades/final_bins_rebin/dereplicated_genomes/gtdbtk.bac120.summary.tsv
        """


rule drep_gtbd_merge:
    input:
        "results/{sample}_spades/final_bins_rebin/data_tables/genomeInformation.csv",
        "results/{sample}_spades/final_bins_rebin/dereplicated_genomes/gtdbtk.bac120.summary.tsv",
        "results/{sample}_spades/final_bins_rebin/dereplicated_genomes/complete_stats.tsv"
    output:
        "results/{sample}_spades/final_bins_rebin/dereplicated_genomes/final_rebin_QC.xlsx"
    conda:
        "envs/python3_modules.yaml"
    threads: 1
    resources: mem_mb=10000, time="1-00:00:00"
    shell:
        """
        python3 workflow/scripts/derep_gtdb_stats_merge.py {input[0]} {input[1]} {input[2]} {output}
        """
