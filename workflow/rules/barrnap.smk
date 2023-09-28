#This rule will run barrnap on all .fa files in a directory (results/bins)

rule barrnap:
	input:
		"results/bins/{sample}.fa"
	output:
		"results/bins/{sample}_barrnap.txt"
	threads: 1
	resources: mem_mb=10000, time="1-00:00:00"
	conda:
		"../envs/barrnap.yaml"
	shell:
		"""
		barrnap {input} --quiet > {output}
		"""
