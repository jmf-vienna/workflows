rule megahit_assembly:
	input:
		"data/interleave/{sample}.interleave.fastq.gz"
	output:
		"results/{sample}_megahit/scaffolds_1000bp.fa"
	threads: 48
	conda:
		"../envs/megahit.yaml"
	resources: mem_mb=500000, time="4-00:00:00", slurm_partition="himem"
	shell:
		"""
		#Output directory cannot exist for megahit
		rm -r results/{sample}_megahit
		
		#assemble individual fastqs with megahit
		megahit --k-min 21 --k-max 121 --k-step 10 -t {threads} -m 500000000000 -o results/{wildcards.sample}_megahit --12 {input}
		reformat.sh in=results/{wildcards.sample}_megahit/final_contigs.fa out={output} minlength=1000
		"""

