
#mapping prep for the binning programs
rule mag_stats:
	input:
		"results/{sample}_spades/final_bins/dereplicated_genomes"
	output:
		"results/{sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv"
	threads: 1
	conda:
		"envs/bbmap_samtools.yaml"
	resources: mem_mb=50000, time="0-01:00:00"
	shell:
                """
		for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*fa; do stats.sh in=$f format=3 > ${f%%.fa}.stats.tsv; done

		touch results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/all_stats.tsv
		cat results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*.stats.tsv >> results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/all_stats.tsv
		for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*.stats.tsv; do tail -1 $f; done >> results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/all_stats.tsv

		touch results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/names.txt
		echo 'genome' >> results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/names.txt
		for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*.stats.tsv; do echo ${f%%.stats.tsv}.fa; done >> results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/names.txt

		paste -d'\t' results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/names.txt results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/all_stats.tsv > results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv
		
		rm results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/names.txt
		rm results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/all_stats.txt
		rm results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*.stats.tsv
                """

rule drep_gtbd_stats_merge:
	input:
		"results/{sample}_spades/final_bins/data_tables/genomeInformation.csv",
		"results/{sample}_spades/final_bins/dereplicated_genomes/gtdbtk.bac120.summary.tsv",
		"results/{sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv"
	output:
		"results/{sample}_spades/final_bins/dereplicated_genomes/final_QC.xlsx"
	conda:
		"envs/python3_modules.yaml"
	threads: 1
	resources: mem_mb=10000, time="0-00:10:00"
	shell:
		"""
		python3 workflow/scripts/derep_gtdb_stats_merge.py {input[0]} {input[1]} {input[2]} {output}
		"""
