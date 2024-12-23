# Rocketchip: A Comprehensive Bioinformatics Workflow for ChIP-Seq Data Analysis

Rocketchip (v1.0.0) is an automated bioinformatics workflow that is capable of analyzing local ChIP-seq data or ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA), the largest publicly available sequence data database. `rocketchip` takes raw data inputs and generates the files required for data visualization and peak delineation.

![Workflow](https://github.com/vhaghani26/rocketchip/blob/main/rocketchip_flowchart.png)

## Table of Contents

* [Installation](#installation)
    * [Setting Up Your Project Directory](#setting-up-your-project-directory)
    * [Creating the Rocketchip Environment](#creating-the-rocketchip-environment)
 	* [Making a Project File](#making-a-project-file)
 * [Running Rocketchip](#running-rocketchip)
    * [Data Storage and Source Code Storage](#data-storage-and-source-code-storage)
    * [Executing Rocketchip](#executing-rocketchip)
* [Interpretting Outputs](#Interpretting-Outputs)
* [Installing Cisgenome](#installing-cisgenome)

## Installation

### Setting Up Your Project Directory

Clone the repository using your project name in the directory you plan to host the project:

```
git clone https://github.com/vhaghani26/rocketchip {project_name}
```

Enter the directory

```
cd {project_name}
```

### Creating the Rocketchip Environment

Prior to starting this, make sure you have Conda installed. Run the following command in your project directory. It will clone the Conda environment with all dependencies needed in order to run the workflow outlined here. This creates an environment called `rocketchip`. If you would like to change the name, feel free to do so where the command says `rocketchip`. Please note that this may take quite a few minutes to run.

```
conda env create -f environment.yml
```

Activate your environment using

```
conda activate rocketchip
```

Run everything downstream of this point in this Conda environment. Note that you must activate this environment every time you restart your terminal.

### Making a Project File

In order to run Rocketchip, you will need to create a project file. A template, `project_file.yaml`, is included in this repository. The contents should look like this:

```
Author:
Project:
Genome:
  Name:
  Location:
Reads:
  Samples:
    grp1:
  Controls:
    ctl1:
Readtype:
Peaktype:
Aligner:
Deduplicator:
Peakcaller:
Threads:
```

* Author - write your name and collaborators' names (if any), but do not exceed one line
* Project - write the name of the project
* Genome - leave blank
    * Name - write the name of the genome you are using (see examples below)
    * Location - if you have a local copy of the genome, put the path to the genome here (e.g. absolute/path/to/my/genome.fa), otherwise put the link corresponding to whatever genome you are using. Here are some commonly used genomes and links that have been proven to work with Rocketchip:

| Organism  | Genome   | Link                                                                        |
| :-------: | :------: | :------------------------------------------------------------------------:  |
| Fruitfly  | dm6      | https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz            |
| Human     | hg19     | https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz      |
| Human     | hg38     | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz |
| Mouse     | mm9      | https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz       |
| Mouse     | mm10     | https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz      |
| Rat       | rn6      | https://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz            |
| Worm      | ce11     | https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz      |
| Yeast     | sacCer3  | https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz   |
| Zebrafish | danRer11 | https://hgdownload.cse.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz  |

* Reads - leave blank
    * Samples - leave blank
        * Sample Groups - put all replicates of a sample in one group, separating samples by group (grp1, grp2, grp3, ...). See note
    * Controls - leave blank
        * Control Groups - leave blank if you are not using a control; if you are using a control, put the replicates of the control in one group. See note
* Readtype - the endedness of the data; options include `single` or `paired`
* Peaktype - this is determined based on whatever element your antibody targets; options include `narrow` or `broad`
    * Note that the only peak-caller explicitly written to handle broad-peak calling is MACS3
* Aligner - software to be used for alignment; options include `bwa_mem`, `bowtie2`, or `STAR`
* Deduplicator - software to be used for deduplication; options include `samtools`, `picard`, `sambamba`, or `no_deduplication`
* Peakcaller - software to be used for peak-calling; options include `macs3`, `genrich`, `pepr`, or `cisgenome`
    * Note that if you are using Cisgenome, it will need to be installed separately (see provided instructions titled "Installing Cisgenome")
* Threads - the number of threads to be used in subsequent analysis steps

Note: For local read data, use the absolute paths as entries. For example, the read path can be `absolute/path/to/my/sample_id`. Do not put the `fq.gz` extension or the read direction (forward vs. reverse) in the path to your read names. Just end the path at the sample name. Based on if you tell Rocketchip whether the data is single- or paired-end, it will match the appropriate files automatically. For read data from the SRA, use the SRA ID for the sample instead of a path. This goes for both the samples and controls.

Here are various examples of project yaml files:
* [One sample with one replicate, no control](https://github.com/vhaghani26/rocketchip_tests/blob/main/run_times/human/human_SRR17409984.yaml)
* [Two samples with three replicates each, no control](https://github.com/vhaghani26/rocketchip_tests/blob/main/cut_and_tag_run/cut_and_tag.yaml)
* [Three samples with two replicates each, no control](https://github.com/vhaghani26/rocketchip_tests/blob/main/cut_and_tag_run/cut_and_run_se.yaml)
* [One sample with two replicates, one control with one replicate](https://github.com/vhaghani26/rocketchip_tests/blob/main/cut_and_tag_run/cut_and_run_pe.yaml)
* [One sample with two replicates, one control with two replicates](https://github.com/vhaghani26/rocketchip_tests/blob/main/replicability/replicability.yaml)

## Running Rocketchip

### Data Storage and Source Code Storage

There are a few considerations to make regarding data storage and source code storage. Briefly, you can either choose to export the `ROCKETCHIP_DATA` and `ROCKETCHIP_SRC` variables, or you can simply use the flags at the command line (the former being more useful, but the latter being easier).

1. `--data`, `ROCKETCHIP_DATA`

If you are using your own data, then use the `--data` option at the command line when you run Rocketchip and set the path to the raw data as the argument (i.e. `--data path/to/your/raw/data/`). If you are planning to use data from the NCBI SRA, then you can either set the `--data` option to whatever directory you want the data to get stored in or you can export `ROCKETCHIP_DATA`. The benefit of exporting `ROCKETCHIP_DATA` is that the data will be stored in the location you designate for whatever analyses you use Rocketchip for. For individual analyses, these files get aliased into the project directory. This means that if you are using the same genome or sample data for multiple analyses, the data is only stored once and not duplicated. This can save both time and storage. To do so, add the `ROCKETCHIP_DATA` variable to your configuration file (`.bashrc`, `.profile`, or whatever file your system uses) like so, ensuring that you change the path to wherever you want to store your data:

```
export ROCKETCHIP_DATA="/share/mylab/raw_data/"
```

If this sounds too complicated, no worries! To keep it simple, just enter your project directory and run Rocketchip with the flag `--data .`, which will store the data in the project directory that you are working in.

2. `--src`, `ROCKETCHIP_SRC`

`ROCKETCHIP_SRC` refers to where the Rocketchip source code is maintained. Essentially, this should just be the path to the `rocketchip` script. You can use the `--src` flag and set the path yourself (i.e. `--src path/to/the/rocketchip/source/code/`). Alternatively, if you plan to use Rocketchip for several projects, it can be helpful to put this in a designated location and export the source code path so you don't have to use the command line argument or reinstall Rocketchip in the future. To do so, add the `ROCKETCHIP_SRC` variable to your configuration file (`.bashrc`, `.profile`, or whatever file your system uses) like so, ensuring that you change the path to wherever you want to store your data:

```
export ROCKETCHIP_SRC="path/to/the/rocketchip/source/code/"
```

In addition to adding or defining a source code location for Rocketchip, you should also export the path:

```
export PATH="$PATH:$ROCKETCHIP_SRC"
```

The path should lead to the directory containing the `rocketchip` script, not the script itself. This will allow Rocketchip to be executed from anywhere at your terminal.

To test if you have set up `ROCKETCHIP_SRC` correctly, restart your terminal or source `.bashrc`/`.profile` or whatever file you typically source from, then run `rocketchip --help`. This should show something like this:

```
usage: rocketchip [-h] [--data <str>] [--src <str>] [--output_file <str>] <path>

Make Snakefiles

positional arguments:
  <path>               Path to configuration file. See README for details

options:
  -h, --help           show this help message and exit
  --data <str>         override/set current ROCKETCHIP_DATA environment variable
  --src <str>          override/set current ROCKETCHIP_SRC environment variable
  --output_file <str>  output snakefile name (default: STDOUT)
```

### Executing Rocketchip

1. **Run Rocketchip**

Enter the directory containing the `project_file.yaml` that you have set up (you can rename this, just make sure to change the name in the command below). Assuming you have set `ROCKETCHIP_DATA` and `ROCKETCHIP_SRC`, all you need to do is run the following:

```
rocketchip project_file.yaml --output_file {output_file_name}
```

If you have not set `ROCKETCHIP_DATA` and `ROCKETCHIP_SRC`, you will need to set them at the command line:

```
rocketchip project_file.yaml --output_file {output_file_name} --data {directory_to_store_the_data} --src {directory_containing_source_code}
```

This will generate the Snakefile you have named `{output_file_name}` that we will run in the next step.

2. **Run Snakemake**

Now, you will run Snakemake. This follows Snakemake's command line usage, but at it's simplest, you can run:

```
snakemake -j 1 -s {output_file_name}
```

Increase `-j` to match the number of jobs you would like to parallelize.

## Interpretting Outputs

There are several output directories, each containing a component of the analysis. These directories are automatically generated when the analysis is run and outputs are automatically sorted into each directory.

### 00_logs
The `00_logs` directory contains output logs. Each log is labeled based on the sample name and rule. The logs can be referenced if the analysis fails at a specific rule. It will contain the run information to be referenced.

### 01_raw_data
All sequence data for both the samples and reference genome, including reference genome alignment files, are aliased in this directory. The files are aliased to the files downloaded in `ROCKETCHIP_DATA` so that downloading and processing only occurs once per genome/sample. Aliases in the local directly allow the user to see the samples and genome used in a specific project's analysis. These aliases are required for Rocketchip to correctly carry out the workflow.

### 02_fastqc_analysis
FastQC analysis (quality control) is carried out on raw sequence data, specifically after conversion from an SRA file to FASTQ file, and again after sequence alignment and processing.

### 03_sam_files
`03_SAM_files` contains the SAM files generated for each sample by the reference genome alignment.

### 04_bam_files
All BAM files are stored in this folder, including intermediates of samtools flagging, sorting, and deduplication. Steps are labeled using tags in the file name.

### 05_bigwig_files
Bigwig files are used for visualization of ChIP-seq data and are one of the final products of the analysis.

### 06_{peakcaller}_peaks
This directory contains the files delineating the peaks. These peaks will be used in answering the biological question you are asking using the data. In many instances, the peaks correspond to the binding sites of a protein of interest.

## Installing Cisgenome

If you are using Cisgenome as your peak caller, you will need to install it [separately](https://www.biostat.jhsph.edu/~hji/cisgenome/index_files/download.htm), as it is not available through Conda. Rocketchip is compatible with version 2.0.

To install it, enter the `tools/` directory of this repository or the directory of your choice (if so, just change the paths appropriately) and carry out the following commands.

1. Download Cisgenome v2.0

```
wget http://jilab.biostat.jhsph.edu/software/cisgenome/executables/cisgenome_v2.0_linux.tar.gz
```

2. Unzip and untar the file

```
tar zvfx cisgenome_v2.0_linux.tar.gz
```

3. Enter the Cisgenome folder

```
cd cisgenome_project/
```

4. Run

```
./makefile
```

Fortunately, the executables work after unzipping and untarring, so if this last step fails, then you can instead add the `bin` directory to your configuration file (e.g. `.bashrc`, `.bash_profile`, `.profile`) like so, making sure to edit the part that says {your_directory} to reflect your directory structure:

```
export PATH=$PATH:{your_directory}/rocketchip/tools/cisgenome_project/bin
```

Now, either source your configuration file or restart your terminal.  To confirm proper installation, you can run:

```
seqpeak
```

This will display the options available for use, indicating that you are able to execute it.

## Questions, Comments, or Concerns?

Feel free to contact me at [vhaghani@ucdavis.edu](vhaghani@ucdavis.edu), and I will get back to you as soon as possible.
