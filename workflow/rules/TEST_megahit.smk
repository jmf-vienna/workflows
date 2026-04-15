#coassembly
rule megahit_coassembly:
        input:
                "data/interleave_DNA/"
        output:
                "results/coassembly/scaffolds_1000bp.fa"
        threads: 32
        conda:
                "../envs/megahit.yaml"
        resources: mem_mb=500000, runtime="4d", slurm_partition="basic"
        shell:
                """
                if [ -d "results/coassembly" ]; then rm -rf results/coassembly; fi
                
                #Generate a list of all interleave fastqs for megahit
                echo "$(ls data/interleave_DNA/*gz | tr '\n' ',')" > intermediates/list.txt
                LIST="$(cat intermediates/list.txt)"
                
                megahit --k-min 21 --k-max 121 --k-step 10 -t {threads} -m 500000000000 -o results/coassembly --12 $LIST --tmp-dir $TMPDIR
                
                sed -i 's/ .*//' final_contigs.fa

                #remove short contigs
                reformat.sh in=results/coassembly/final_contigs.fa out={output} minlength=1000
                """