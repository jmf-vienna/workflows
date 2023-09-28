#This rule will run gtdbtk on a directory (results/bins) and assign taxonomy to all .fa files

#Requires the GTDBPATH variable from configs (and a download)
configfile: "config/general_configs.yaml"

rule gtdb_binning:
	input:
		"results/bins"
	output:
		"results/bins"
	threads: 16
	conda:
		"../envs/gtdbtk.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	shell:
		"""
		GTDBTK_DATA_PATH="{config[GTDBPATH]}"
		gtdbtk classify_wf --extension fa --cpus {threads} --genome_dir results/bins/ --out_dir results/bins/gtdb_classification --skip_ani_screen
		"""

