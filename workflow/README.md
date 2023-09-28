# JMF independent workflows 
This repository contains dependency-independent workflows that can be shared and used on Joint Microbiome Facility (University of Vienna and Medical University of Vienna) generated data.
It is regularly changed (for now) and that frequently results in some broken code... sorry. 

These workflows are designed to be run on a SLURM-based server system with snakemake but should be able to run on any linux system as all of the dependencies are installed using bioconda through snakemake.

## Data structure
For reproducibility these directories are designed around snakemake's directory structure (found here: https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) 

That looks like this:
```
├── .gitignore
├── README.md
├── LICENSE.md
├── workflow
│   ├── rules
|   │   ├── module1.smk
|   │   └── module2.smk
│   ├── envs
|   │   ├── tool1.yaml
|   │   └── tool2.yaml
│   ├── scripts
|   │   ├── script1.py
|   │   └── script2.R
│   ├── notebooks
|   │   ├── notebook1.py.ipynb
|   │   └── notebook2.r.ipynb
│   ├── report
|   │   ├── plot1.rst
|   │   └── plot2.rst
|   └── Snakefile
├── log
├── intermediates
├── config
│   ├── config.yaml
│   └── some-sheet.tsv
├── data
│   ├── DNA
│   ├── RNA
│   ├── NP
│   └── ref
├── results
└── resources
```

But don't worry, you can set this entire thing up like this:
```
mkdir <desired directory name>
cd <desired directory name>

#clone this repo
git clone https://github.com/osvatic/JMF_independent_workflows.git .

#make the rest of the directories
mkdir data
mkdir data/DNA
mkdir data/RNA
mkdir data/ref
mkdir data/NP
mkdir log
mkdir intermediates
mkdir results
mkdir resources

```

## How to run
In general snakemake workflows/rules can be run from the main directory ("desired directory name" from above) like this:

**For a system with modules**
```
module load conda
module load snakemake
snakemake --cores <total available threads> -s <snakemake workflow/rule> -j10 --use-conda
```
If you are not using this on a using with modules, please ignore the "module load" commands.
In that case you need to have conda and snakemake installed and accessible. 


**To submit the entire workflow or rule to a slurm-based workflow manager**
```
sbatch workflows/scripts/runSnakemake.sh <workflow/rule>
```
The `runSnakemake.sh` will submit all workflow jobs to the servers as needed. 


## Details on individual workflows

### Quality checking, adapter and phiX removal from bam files
These workflows can be used with paired-end (1 bam file) illumina reads in the data/DNA directory. The workflow will clean all .bam files and use the sample name from data/DNA/sample.bam to label the output interleaved fastq.gz. The final readsets will be in data/interleave_DNA.

For environmental samples (non-medical):
`qc_DNA_frombam_nonmedical.smk`

For medical samples (removes reads aligned to mouse, rat and human genomes):
`qc_DNA_frombam_medical.smk`


Similarly, RNA files can be processed from data/RNA. This workflow is similar to the DNA version but has a stricter q-value cut off (q=28 instead of q=15).

For environmental samples (non-medical):
`qc_RNA_frombam_nonmedical.smk`

For medical samples (removes reads aligned to mouse, rat and human genomes):
`qc_RNA_frombam_medical.smk`


### Quality checking, adapter and phiX removal from split fastq files
These workflows can be used with paired-end (split files) illumina reads in the data/DNA directory. The workflow will clean all .fastq.gz files and use the sample name from data/DNA/sample.1.fastq.gz to label the output interleaved fastq.gz. Readsets should named like this "samplename.1.fastq.gz" and "samplename.2.fastq.gz" for the left and right read files respectively. The final readsets will be in data/interleave_DNA.

For environmental samples (non-medical):
`qc_DNA_fromsplitfq_nonmedical.smk`

For medical samples (removes reads aligned to mouse, rat and human genomes):
`qc_DNA_fromsplitfq_medical.smk`


### Metagenome assembly and binning
The workflows all use SPAdes (https://cab.spbu.ru/files/release3.15.2/manual.html) as the assembler and then attempt to bin the metagenome into putative MAGs using Metabat2 (https://bitbucket.org/berkeleylab/metabat/src/master/), both with coverage and without. All workflows will attempt to assemble all data/interleave_DNA/sample.interleave.fastq.gz files and use some combination of those same files for coverage binning. 

For metagenomic assembly and binning with all readsets in data/interleave_DNA/:
`spades_metagenome_assembly_binning_noconcoct.smk`

For metagenomic assembly and binning with only the sample's readset in data/interleave_DNA/:
`spades_metagenome_assembly_binning_noconcoct_individualsample_binning.smk`

For isolate (or very low diversity) assembly and binning with only the sample's readset in data/interleave_DNA/:
`spades_isolate_assembly_binning_noconcoct_individualsample_binning.smk`


#### Output
Final MAGs (high quality ones) can be found here: **results/{sample}_spades/final_bins_rebin/dereplicated_genomes/**

A quality and taxonomy excel file can be found here for all of these final MAGs: **results/{sample}_spades/final_bins_rebin/dereplicated_genomes/final_rebin_QC.xlsx**


The workflows `megahit_metagenome_assembly_binning.smk` and `spades_metagenome_assembly_binning.smk` are not currently in use because they use concoct, which I cannot get to work at the moment. 

### Binning metagenomes using a different combination of readsets for coverage
Many projects have multiple treatment groups or subsets with the project, which could make binning with only a subset more successful (or less) than binning with all data/interleave_DNA/sample.fastq.gz. 

A subset can be used after running `spades_metagenome_assembly_binning_noconcoct.smk`. Upon completion, all bams for the mapping will be in **results/{sample}_spades/bams**. Just move any mapping files (named after samples) to a newly created directory called **results/{sample}_spades/bams_rebin**. Then you can run `rebinner_spades.smk`.

This is not the perfect set up as it wastes some time and computing power aligning more samples to the metagenome than needed.... but it works. 

#### Output
Assemblies can be found here upon completion: **results/{sample}_spades/** 

Final MAGs (high quality ones) can be found here: **results/{sample}_spades/final_bins/dereplicated_genomes/**

A quality and taxonomy excel file can be found here for all of these final MAGs: **results/{sample}_spades/final_bins/dereplicated_genomes/final_QC.xlsx**


### Microbial Transcriptomic analysis
Files for RNA DESeq/edgeR analysis can be prepared using:
`transcriptome_featurecounts_generation.smk`

Bacterial MAGs or Assemblies can be added to **data/ref** with the ending .fa and RNA readsets can be added to **/data/interleave_RNA** (should be in there if you used the above RNA workflows). 
This workflow will annotation all MAGs in data/ref using PROKKA and align all readsets to them using bbmap and generate the appropriate FeatureCounts tables. 

#### Output: 
The FeatureCounts tables can be found here: "results/FeatureCounts/{genome}.featurecounts.tsv"

The PROKKA annotations can be found here: data/ref/{genome name}

## Coming soon
Nanopore workflows

## Comments and questions

Please feel free to share comments/issues/suggestions with me.


