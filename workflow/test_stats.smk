
#generating mag stats for all things
rule mag_stats:
	input:
		"results/{sample}_spades/final_bins/dereplicated_genomes"
	output:
		"results/{sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv"
	threads: 1
	conda:
		"envs/bbmap_samtools.yaml"
	resources: mem_mb=50000, time="0-01:00:00", partition="basic"
	shell:
                """
		#generate stats file for each MAG
		for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*fa;
		do stats.sh in=$f format=3 > ${f%%.fa}.stats.tsv;
		echo 'genome' > ${f%%.fa}.name.txt;
		echo "$f" >> ${f%%.fa}.name.txt;
		paste -d'\t' ${f%%.fa}.name.txt ${f%%.fa}.stats.tsv > ${f%%.fa}.stats.merged.tsv;
		rm ${f%%.fa}.stats.tsv;
		rm ${f%%.fa}.name.txt;
		done
		
		#get headers for stats
		cat results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*stats.merged.tsv | head -1 > results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/complete_stats.tsv

		#merge all stats into 1 file
		for f in results/{wildcards.sample}_spades/final_bins/dereplicated_genomes/*stats.merged.tsv;    
		do tail	-1 $f >> complete_stats.tsv;
		rm $f;
		done

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
	resources: mem_mb=10000, time="0-00:10:00", partition="basic"
	shell:
		"""
		python3 workflow/scripts/derep_gtdb_stats_merge.py {input[0]} {input[1]} {input[2]} {output}
		"""
