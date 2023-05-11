#This is a generalized nanopore assembly workflow that assembles and polishes the data

configfile: "config/general_configs.yaml"


FILES = glob_wildcards('data/NP/{name}.fastq.gz')
NAMES = FILES.name


rule complete:
	input:
		expand("results/{name}_flye/assembly_medaka1x.fa.gz", name=NAMES)


rule flye_assembly:
	input:
		"data/NP/{sample}.fastq.gz"
	output:
		"results/{sample}_flye/assembly_1000bp.fa"
	threads: 32
	conda:
		"envs/flye.yaml"
	resources: mem_mb=100000, time="4-00:00:00"
	shell:
		"""
		flye --nano-hq {input} --threads {threads} --meta --out-dir results/{wildcards.sample}_flye
		seqtk seq -L 1000 results/{wildcards.sample}_flye/assembly.fasta > {output}
		"""


rule medaka_polish:
	input:
		"results/{sample}_flye/assembly_1000bp.fa"
	output:
		"results/{sample}_flye/assembly_medaka1x.fa.gz"
	threads: 8
	conda:
		"envs/medaka.yaml"
	resources: mem_mb=100000, time="4-00:00:00", partition="basic"
	shell:
		"""
		MEDAKAMODEL=r1041_e82_400bps_sup_g615
		medaka_consensus -i data/NP/{wildcards.sample}.fastq.gz -b 20 -d {input} -o results/{wildcards.sample}_flye -t {threads} -m $MEDAKAMODEL
		cat results/{wildcards.sample}_flye/consensus.fasta | gzip > {output}
		rm results/{wildcards.sample}_flye/consensus.fasta
		"""
