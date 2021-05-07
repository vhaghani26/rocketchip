# MeChIP2

## What is MeCHIP2?
MeCHIP2 is an automated bioinformatics pipeline that downloads raw paired-end ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA). ChIP-seq data is downloaded directly from the SRA using the SRA Toolkit (sra-tools, v2.10.9). The raw sequence read file is then split into its respective paired-end read files by the SRA Toolkit. Paired-end reads are aligned to the GRCm39/mm39 mouse genome assembly using the Burrows-Wheeler Algorithm (BWA) Maximal Exact Match (MEM) software (bwa, v0.7.17). Samtools (samtools, v1.12) is used for file format conversion and deduplication of sequence data. Finally, deeptools (deeptools, v3.5.1) is used to convert data to the bigwig file format, which can be used for visualization of ChIP-seq data in the UCSC Genome Browser or other visualization tools.

## Using the Pipeline
In order to use the pipeline, clone the GitHub repository (https://github.com/vhaghani26/MeChIP2.git). Open the ```samples.txt``` file and add each SRA ID you are interested in visualizing on its own line. 

[insert information about snakefile use here]

## Outputs
There will be four output directories, each containing a component of the pipeline.

[add more in-depth info here about what the file is]

### 01_raw_data
All sequence data for both the samples and reference genome, including reference genome alignment files, are stored in this directory.

### 02_sam_files
SAM files contains the SAM files generated for each sample by the reference genome alignment.

### 03_bam_files
All BAM files are stored in this folder, including intermediates of samtools flagging, sorting, and deduplication.

### 04_bigwig_files
Bigwig files are used for visualization of ChIP-seq data and are the final product of the pipeline.
