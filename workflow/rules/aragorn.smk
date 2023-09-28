#This rule will run aragorn on all .fa files in a directory (results/bins)
rule aragorn:
	input:
		"results/bins/{sample}.fa"
	output: 
		"results/bins/{sample}_aragorn.txt"
	threads: 1
	resources: mem_mb=10000, time="1-00:00:00"
	conda:
		"../envs/aragorn.yaml"
	shell:
		"""
		aragorn {input} -t -d > {output}
		"""

