#This rule will run phyloflash on interleaved fastqs


#remember to change the read length!!!
rule phyloflash:
	input:
		"data/interleave/{sample}.interleave.fastq.gz
	output:
		directory("results/{sample}_phyloflash")
	threads: 16
	conda:
		"../envs/phyloflash.yaml"
	resources: mem_mb=10000, time="1-00:00:00"
	shell:
		"""
		cd results/
		phyloFlash.pl -lib {wildcards.sample}_phyloflash -CPUs 16 -everything -read1 ../{input} --interleaved -readlength 100 -clusterid 99
		cd ../
		"""
