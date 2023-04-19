#workflow for spades assembler and binning process afterwards. It will 

FILES = glob_wildcards('data/interleave/{name}.interleave.fastq.gz')
NAMES = FILES.name

configfile: "config/general_configs.yaml"

#This rule is just here to "request" the final results and set everything into motion
rule complete:
	input:
		expand("results/{name}_spades/final_bins/dereplicated_genomes/final_QC.csv", name=NAMES)


#Assemble with SPAdes using interleaved reads, will only keep >1000bp contigs
rule SPAdes_assembly:
	input:
		"data/interleave/{sample}.interleave.fastq.gz"
	output:
		"results/{sample}_spades/scaffolds_1000bp.fa"
	threads: 48
	conda:
		"envs/spades.yaml"
	resources: mem_mb=750000, time="4-00:00:00", partition="himem"
	log: "log/SPADES_{sample}.log"
	shell:
		"""
		spades.py -t {threads} -m 750 -k 21,31,41,51,61,71,81,91 --meta --pe-12 1 {input} -o results/{wildcards.sample}_spades 2> {log}
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
	log: "log/mapping_{sample}.log"
	shell:
                """
                if [ ! -d "results/{wildcards.sample}_spades/bams" ]; then mkdir results/{wildcards.sample}_spades/bams; fi
                
                for f in data/interleave/*.interleave.fastq.gz; do g=${{f##*/}}; bbmap.sh -Xmx50g threads={threads} ref={input} interleaved=true nodisk minid=0.98 in=$f out=results/{wildcards.sample}/bams/${{g%%.interleave.fastq.gz}}.bam; done 2> {log}

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

#Concoct with coverage
rule concoct:
	input:
		"results/{sample}_spades/scaffolds_1000bp.fa",
		"results/{sample}_spades/bams/"
	output:
		directory("results/{sample}_spades/concoct/fasta_bins/")
	threads: 1
	conda:
		"envs/concoct.yaml"
	resources: mem_mb=50000, time="1-00:00:00"
	log: "log/concoct_{sample}.log"
	shell:
                """
                if [ ! -d "results/{wildcards.sample}_spades/concoct" ]; then mkdir results/{wildcards.sample}_spades/concoct; fi

                cut_up_fasta.py {input[0]} -c 10000 -o 0 --merge_last -b results/{wildcards.sample}_spades/concoct/contigs_10K.bed > results/{wildcards.sample}_spades/concoct/contigs_10K.fa
                concoct_coverage_table.py results/{wildcards.sample}_spades/concoct/contigs_10K.bed results/{wildcards.sample}_spades/bams/*.sorted.bam > results/{wildcards.sample}_spades/concoct/coverage_table.tsv 2> {log}

                #concoct only likes to run on 1 thread
                concoct --composition_file results/{wildcards.sample}_spades/concoct/contigs_10K.fa --coverage_file results/{wildcards.sample}_spades/concoct/coverage_table.tsv -b results/spades_assembly/concoct --threads 1 2>> {log}
                merge_cutup_clustering.py results/{wildcards.sample}_spades/concoct/clustering_gt1000.csv > results/{wildcards.sample}_spades/concoct/clustering_merged.csv 2>> {log}

                if [ ! -d "results/{wildcards.sample}_spades/concoct/fasta_bins" ]; then mkdir results/{wildcards.sample}_spades/concoct/fasta_bins; fi
                extract_fasta_bins.py {input[0]} results/{wildcards.sample}_spades/concoct/clustering_merged.csv --output_path results/{wildcards.sample}_spades/concoct/fasta_bins 2>> {log}
                """

#A quick dREP
rule drep:
	input:
		"results/{sample}_spades/metabat_nocov/",
		"results/{sample}_spades/metabat_cov/",
		"results/{sample}_spades/concoct/fasta_bins/"

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
                cp results/{wildcards.sample}_spades/concoct/fasta_bins/*fa results/{wildcards.sample}_spades/final_bins

                dRep dereplicate results/{wildcards.sample}_spades/final_bins/ -p {threads} -g results/{wildcards.sample}_spades/final_bins/*fa 2> {log}
                
                """


#reanalyze things
#at this point lets aim for gtdb and add to the dREP
rule gtdb_binning:
	input:
		"results/{sample}_spades/final_bins/data_tables/genomeInfo.csv"
	output:
		"results/{sample}_spades/final_bins/dereplicated_genomes/final_QC.csv"
	threads: 16
	conda:
		"envs/gtdbtk.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	log: "log/gtdbtk_{sample}.log"
	shell:
                """
		conda env config vars set GTDBTK_DATA_PATH={config[GTDBPATH]}
		gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes --out_dir results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification --skip_ani_screen 2> {log}
                mv results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdb_classification/classify/gtdbtk.bac120.summary.tsv results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv

                #make final csv
                echo "genome,completeness,contamination,strain_hetero,length,N50,GTDB_id" > results/spades_assembly/final_bins/dereplicated_genomes/final_QC.csv

                for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*fa; do g=${{f##*/}}; DREP="$(grep "^$g" results/{wildcards.sample}_spades/final_bins/data_tables/genomeInfo.csv)"; GTD="$(grep "^${{g%%.fa}}" results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv | cut -f2)"; echo "${{DREP}},${{GTD}}" >> results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/final_QC.csv; done
                """
