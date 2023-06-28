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
### Starting from bam files
qc_DNA_frombam_nonmedical.smk
qc_DNA_frombam_medical.smk


qc_RNA_frombam_medical.smk
qc_RNA_frombam_nonmedical.smk

qc_DNA_fromsplitfq_medical.smk
qc_DNA_fromsplitfq_nonmedical.smk




## Comments and questions

Please feel free to share comments/issues/suggestions with me.
