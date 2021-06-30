# MeChIP2: A Comprehensive ChIP-Seq Data Analysis Bioinformatics Pipeline

## What is MeCHIP2?
MeCHIP2 is an automated bioinformatics pipeline that downloads raw paired-end ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) and uses it to generate the files required for data visualization and peak delineation. ChIP-seq data is downloaded directly from the SRA using the SRA Toolkit (sra-tools, v2.10.9). The raw sequence read file is then split into its respective paired-end read files by the SRA Toolkit. Paired-end reads are aligned to the GRCm39/mm39 mouse genome assembly using the Burrows-Wheeler Algorithm (BWA) Maximal Exact Match (MEM) software (bwa, v0.7.17). Samtools (samtools, v1.12) is used for file format conversion and deduplication of sequence data. Deeptools (deeptools, v3.5.1) is used to convert data to the bigwig file format, which can be used for visualization of ChIP-seq data in the UCSC Genome Browser or other visualization tools. Peaks are called using MACS2 (macs2, v2.2.7.1). FastQC (fastqc, v0.11.9) carries out a sequence quality control analysis pre-processing, on the raw sequence data, and post-processing, after sequence alignment takes place. Overall, this pipeline carries out the major steps of ChIP-seq data analysis and generates output files to be used as figures and to be used in further analysis.

## Using the Pipeline
In order to use the pipeline, clone the [GitHub repository](https://github.com/vhaghani26/MeChIP2.git) into the directory of your choice. The contents of the repository are broken down below for your reference. See the "Snakefile" section under "Repository Contents" for more specifics on file execution. It is  recommended to try the tutorial outlined at the end of this document to better understand specific usage of the program. Additionally, Snakemake should be installed in order to run the pipeline. See [here](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) for instructions on Snakemake installation.

### Repository Contents

#### samples.yaml
Open the `samples.yml` file and add each SRA ID you are interested in visualizing on its own line using the following format:
```
---
samples: 
   - SRAID1
   - SRAID2
```
This file will be used as a configuration file to pipe in the samples used in the analysis. The syntax in this file is important, so be aware that errors relating to `wildcards` are likely caused by the entries in this file.

#### Snakefile
The Snakefile contains the code required to execute the entire pipeline. Here, I will go through the basics of Snakemake usage for this program. For additional information, see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). For specifics regarding command-line interface, [this](https://snakemake.readthedocs.io/en/stable/executing/cli.html) may be helpful.

Briefly, the whole program can be executed by running `snakemake -j 4 -p`, where the `-j` or `--jobs` option specifies the highest number of jobs (highest number of CPU) that can be run in parallel. The `-p` or `--prioritize` option tells the job scheduler to prioritize certain jobs given their dependencies. Without specifying anything else, the whole workflow will run.

To run a rule in particular, you can use `snakemake -j 4 -p rulename_wc`. The `_wc` tag at the end of the rule name is required when running an individual rule to ensure that your `wildcards` (samples) are properly incorporated. If the whole workflow is being run, the `_wc` tag is not necessary.

The list of rules and their descriptions goes as follows:
- `rule all`: runs the entire workflow
- `rule make_directories`: makes directories for data organization. See "Outputs" section for more information.
- `rule download_data`: downloads raw sequence data files (from SRA IDs)
- `rule split_paired_reads`: splits paired-end read data into separate files
- `rule fastqc_precheck`: run quality control on samples pre-processing
- `rule download_genome`: download the GRCm39/mm39 mouse genome from the UCSC Genome Browser
- `rule process_genome`: decompress the genome and concatenate individual chromosome files to create the full assembly
- `rule set_alignment_reference`: set the GRCm39/mm39 mouse genome assembly as the reference genome for alignment
- `rule align_reads`: align paired end reads to GRCm39/mm39 reference genome
- `rule sam_to_bam`: convert SAM to BAM file format
- `rule sam_fixmate`: remove secondary and unmapped reads and add tags to reads for deduplication
- `rule sam_sort`: sort reads by chromosome coordinates
- `rule sam_markdup`: mark and remove duplicates
- `rule sam_index`: index deduplicated BAM file
- `rule bam_to_bigwig`: convert BAM file format to bigwig file format for visualization
- `rule call_peaks`: call ChIP-seq peaks
- `rule fastqc_postprocessing`: run quality control on samples post-processing

These rules are loosely written in the order they get executed. If you want to run the whole pipeline up to a certain point, using the `snakemake -j 4 -p rulename_wc` command will determine all the dependencies of that rule and generate the output required by the specified rule. This means everything needed beforehand will be run and their outputs saved, but anything beyond the specified rule will not be run.

#### 00_conda_software
The repository contains a folder called `00_conda_software` that contains individual `yml` files with each specified software, including its version number, to be used for each rule in the pipeline. This ensures replicability and reproducibility of analysis results.

#### .gitignore
The `.gitignore` file contains the directories and consequently the data generated by each rule. Because the data is not pushed to git, it is strongly recommended to back up your data somewhere that it is protected.
 
#### dag_mechip2.pdf
`dag_mechip2.pdf` was generated using the command `snakemake --dag | dot -Tpdf > dag_mechip2.pdf`. The file is a visual depiction of the workflow that is generated based on the dependencies of each rule.

#### mechip2.slurm
The `mechip2.slurm` file is a sample SLURM file that can be used as a template for those wishing to run the workflow on the cluster/HPC.

## Outputs
There are several output directories, each containing a component of the pipeline. These directories are automatically generated when the analysis is run and outputs are automatically sorted into each directory.

### 00_conda_software
See the "00_conda_software" section under "Repository Contents" above.

### 00_logs
The `00_logs` directory contains output logs. Each log is labeled based on the sample name and rule. The logs can be referenced if the pipeline fails at a specific rule. It will contain the run information to be referenced. 

### 01_raw_data
All sequence data for both the samples and reference genome, including reference genome alignment files, are stored in this directory.

### 02_fastqc_analysis
FastQC analysis (quality control) is carried out on raw sequence data, specifically after paired end read data are split, and again after sequence alignment and processing. 

### 03_sam_files
`03_SAM_files` contains the SAM files generated for each sample by the reference genome alignment.

### 04_bam_files
All BAM files are stored in this folder, including intermediates of samtools flagging, sorting, and deduplication. Steps are labeled using tags in the file name.

### 05_bigwig_files
Bigwig files are used for visualization of ChIP-seq data and are the final product of the pipeline. Although bedGraphToBigWig can be used to convert the {sample}_treat_pileup.bdg files generated by MAS2 to bigwig format for visualization, deeptools was chosen to generate bigwig files because it allows for greater flexibility in the peak-calling software used if a user decides to change the peak-caller. 

### 06_macs2_peaks

## Tutorial

