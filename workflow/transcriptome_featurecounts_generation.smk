#This workflow will generate annotations for a reference, align reads to reference, make featurecounts table, and use DESeq for analysis 
#A metadata csv is also needed
#THIS currently dislikes "_" in the fastq.gz names

#this should only be 1 file
READS = glob_wildcards('data/interleave_RNA/{library}.interleave.fastq.gz')
LIBS = READS.library

GENOMES = glob_wildcards('data/ref/{genome}.fa')
MAGS = GENOMES.genome


#This rule is just here to "request" the final results and set everything into motion
rule complete:
	input:
		expand("results/FeatureCounts/{genome}.featurecounts.tsv", genome=MAGS)

#PROKKA
rule PROKKA:
	input:
		"data/ref/{genome}.fa"
	output:
		"data/ref/{genome}/{genome}.gff"
	threads: 8
	conda:
		"envs/prokka.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	shell:
		"""
		prokka --cpus {threads} --outdir data/ref/{wildcards.genome} --prefix {wildcards.genome} {input} --force
		"""


#ALIGN RNA to ref
rule BBMAP_align:
	input:
		"data/ref/{genome}.fa",
		"data/interleave_RNA/{library}.interleave.fastq.gz"
	output:
		"data/ref/{genome}_{library}.bam"
	threads: 16
	conda:
		"envs/bbmap.yaml"
	resources: mem_mb=50000, time="1-00:00:00"
	shell:
		"""
		bbmap.sh -Xmx50g threads={threads} minid=0.99 ambiguous=best interleaved=true ref={input[0]} nodisk in={input[1]} out={output}
		"""
#FEATUREC
rule FeatureCounts:
	input:
		"data/ref/{genome}/{genome}.gff", 
		expand("data/ref/{genome}_{library}.bam", library=LIBS, genome=MAGS)
	output:
		"results/FeatureCounts/{genome}.featurecounts.tsv"
	threads: 8
	conda:
		"envs/subread.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	shell:
		"""
		featureCounts -s 2 -p -T 8 -O -g ID -t CDS -a {input[0]} -o {output} data/ref/{wildcards.genome}_*.bam 
		"""
