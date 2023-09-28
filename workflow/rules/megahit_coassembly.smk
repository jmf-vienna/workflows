rule megahit_assembly:
	input:
		"data/interleave/"
	output:
		"results/coassembly_megahit/scaffolds_1000bp.fa"
	threads: 48
	conda:
		"../envs/megahit.yaml"
	resources: mem_mb=500000, time="4-00:00:00", slurm_partition="himem"
	shell:
		"""
		#Output directory cannot exist for megahit to run
		rm -r results/coassembly_megahit
		
		#Generate a list of all interleave fastqs for megahit
		echo "$(ls data/interleave/*gz | tr '\n' ',')" > intermediates/list.txt
		LIST="$(sed 's/\(.*\),/\1 /' intermediates/list.txt)"
		
		#assemble with megahit
		megahit --k-min 21 --k-max 121 --k-step 10 -t {threads} -m 500000000000 -o results/coassembly_megahit --12 LIST
		reformat.sh in=results/coassembly_megahit/final_contigs.fa out={output} minlength=1000
		"""
