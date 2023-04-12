#This rule will run checkm on a directory (results/bins) and assign taxonomy to all .fa files

rule checkM:
	input:
		"results/bins/"
	output:
		"results/bins/checkM_stats.txt"
	threads: 16
	resources: mem_mb=10000, time="1-00:00:00"
	conda:
		"../envs/checkm.yaml"
	shell:
		"""
		checkm lineage_wf -t {threads} --pplacer_threads {threads} -f {output} -x fa results/bins results/bins/checkm
		rm -r results/spades_assembly/checkm
		"""

