#workflow for spades assembler and binning process afterwards. It will 

FILES = glob_wildcards('data/interleave/{name}.interleave.fastq.gz')
NAMES = FILES.name

configfile: "config/general_configs.yaml"

#This rule is just here to "request" the final results and set everything into motion
rule complete:
	input:
		expand("results/{name}_spades/final_bins/dereplicated_genomes/final_QC.xlsx", name=NAMES)


#Assemble with SPAdes using interleaved reads, will only keep >1000bp contigs
rule SPAdes_assembly:
	input:
		"data/interleave/{sample}.interleave.fastq.gz"
	output:
		"results/{sample}_spades/scaffolds_1000bp.fa"
	threads: 48
	conda:
		"envs/spades.yaml"
	resources: mem_mb=500000, time="4-00:00:00"
	log: "log/SPADES_{sample}.log"
	shell:
		"""
		spades.py -t {threads} -m 500 -k 21,31,41,51,61,71,81,91 --meta --pe-12 1 {input} -o results/{wildcards.sample}_spades 2> {log}
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
	resources: mem_mb=50000, time="1-00:00:00"
	log: "log/metabatnocov_{sample}.log"
	shell:
                """
                if [ ! -d "results/{wildcards.sample}_spades/metabat_nocov" ]; then mkdir results/{wildcards.sample}_spades/metabat_nocov; fi
                metabat2 -m 1500 -t {threads} -i {input} -o results/{wildcards.sample}_spades/metabat_nocov/metabat_nocov 2> {log}
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
	resources: mem_mb=50000, time="1-00:00:00"
	shell:
                """
                if [ ! -d "results/{wildcards.sample}_spades/bams" ]; then mkdir results/{wildcards.sample}_spades/bams; fi
                
                bbmap.sh -Xmx50g threads={threads} ref={input} interleaved=true nodisk minid=0.98 in=data/interleave/{wildcards.sample}.interleave.fastq.gz out=results/{wildcards.sample}_spades/bams/{wildcards.sample}.bam

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
	resources: mem_mb=50000, time="1-00:00:00"
	log: "log/metabat_{sample}.log"
	shell:
		"""
		if [ ! -d "results/{wildcards.sample}_spades/metabat_cov" ]; then mkdir results/{wildcards.sample}_spades/metabat_cov; fi
		jgi_summarize_bam_contig_depths --outputDepth results/{wildcards.sample}_spades/metabat_cov/depth.txt results/{wildcards.sample}_spades/bams/*sorted.bam
		metabat2 -i {input[0]} -o results/{wildcards.sample}_spades/metabat_cov/metabat_cov -m 1500 -t {threads} -a results/{wildcards.sample}_spades/metabat_cov/depth.txt 2> {log}
		"""


#A quick dREP
rule drep:
	input:
		"results/{sample}_spades/metabat_nocov/",
		"results/{sample}_spades/metabat_cov/"

	output:
		"results/{sample}_spades/final_bins/data_tables/genomeInfo.csv"
	threads: 8
	conda:
		"envs/drep.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	log: "log/drep_{sample}.log"
	shell:
                """
                if [ ! -d "results/{wildcards.sample}_spades/final_bins" ]; then mkdir results/{wildcards.sample}_spades/final_bins; fi
                cp results/{wildcards.sample}_spades/metabat_cov/metabat_cov*fa results/{wildcards.sample}_spades/final_bins
                cp results/{wildcards.sample}_spades/metabat_nocov/metabat_nocov*fa results/{wildcards.sample}_spades/final_bins

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
		"results/{sample}_spades/final_bins/data_tables/genomeInformation.csv",
		"results/{sample}_spades/final_bins/GTDB_env.txt"
	output:
		"results/{sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv"
	threads: 16
	conda:
		"envs/gtdbtk.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	log: "log/gtdbtk_{sample}.log"
	shell:
                """
		gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes --out_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification --skip_ani_screen 2> {log}
                mv results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification/classify/gtdbtk.bac120.summary.tsv results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv
                """


rule drep_gtbd_merge:
	input:
		"results/{sample}_spades/final_bins/data_tables/genomeInformation.csv",
		"results/{sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv"
	output:
		"results/{sample}_spades/final_bins/dereplicated_genomes/final_QC.xlsx"
	conda:
		"workflow/envs/python3_modules.yaml"
	threads: 1
	resources: mem_mb=10000, time="1-00:00:00"
	shell:
	"""
	python3 workflow/scripts/derep_gtdb_merge.py {input[0]} {input[1]} {output}
	"""

