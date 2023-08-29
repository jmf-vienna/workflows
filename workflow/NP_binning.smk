#This is a generalized nanopore binning script

configfile: "config/general_configs.yaml"


FILES = glob_wildcards('data/NP/{name}.fastq.gz')
NAMES = FILES.name


rule complete:
	input:
		expand("results/{name}_flye/mapping/checkM_stats.txt", name=NAMES)


rule minimap:
	input:
		"results/{sample}_flye/assembly_1000bp.fa",
		"data/NP/{readset}.fastq.gz"
	output:
		"results/{sample}_flye/mapping/{readset}.np.cov.bam"
	threads: 16
	conda:
		"envs/minimap.yaml"
	resources: mem_mb=100000, time="4-00:00:00", partition="basic"
	shell:
		"""
		if [ ! -d "results/{wildcards.sample}_flye/mapping" ]; then mkdir results/{wildcards.sample}_flye/mapping; fi
		minimap2 -ax sr -t {threads} {input[0]} {input[1]} | samtools view --threads {threads} -Sb -F 0x104 - | samtools sort --threads {threads} - > results/{wildcards.sample}_flye/mapping/{wildcards.readset}.np.cov.bam
		"""


rule metabat:
	input:
		"results/{sample}_flye/assembly_1000bp.fa", 
		"results/{sample}_flye/mapping/{sample}.np.cov.bam"
	output:
		directory("results/{sample}_flye/mapping/metabat_cov/")
	threads: 16
	conda:
		"envs/metabat.yaml"
	resources: mem_mb=100000, time="4-00:00:00", partition="basic"
	shell:
		"""
                jgi_summarize_bam_contig_depths --percentIdentity 80 --outputDepth results/{wildcards.sample}_flye/mapping/depth.txt results/{wildcards.sample}_flye/mapping/*.np.cov.bam
                metabat2 -i {input[0]} -o results/{wildcards.sample}_flye/mapping/metabat_cov/metabat_cov -m 1500 -t {threads} -a results/{wildcards.sample}_flye/mapping/depth.txt
		"""

rule checkm:
	input:
		"results/{sample}_flye/mapping/metabat_cov/"
	output:
		"results/{sample}_flye/mapping/checkM_stats.txt"
	threads: 16
	conda:
		"envs/checkm_genome.yaml"
	resources: mem_mb=50000, time="1-00:00:00", partition="basic"
	shell:
		"""
		checkm lineage_wf -t {threads} --pplacer_threads {threads} -f {output} -x fa {input} {input}/checkm
		rm -r {input}/checkm
		"""
