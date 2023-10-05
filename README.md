# Rocketchip: A Comprehensive Bioinformatics Workflow for ChIP-Seq Data Analysis

v0.0.3

**Rocketchip is currently being redone. The README.md will be redone shortly to include updated usage instructions and a tutorial for getting started.**

_______________________________________________________________________________________________

## Table of Contents

* [Project Set-Up](#project-set-up)
	* [Setting Up Your Project Directory](#setting-up-your-project-directory)
	* [Installation](#installation)
* [Running `rocketchip`](#running-rocketchip)
    * [How Does `rocketchip` Work?](#how-does-rocketchip-work)
	* [Making a Project File](#making-a-project-file)
	* [Running `rocketchip` Locally](#running-rocketchip-locally)
	* [Running `rocketchip` on SLURM (Recommended)](#Running-rocketchip-on-SLURM-Recommended)
* [Tutorial](#tutorial)
* [Interpretting Outputs](#Interpretting-Outputs)

## Project Set-Up

### Setting Up Your Project Directory

Clone the repository using your project name in the directory you plan to host the project:

```
git clone https://github.com/vhaghani26/rocketchip {project_name}
```

Enter the directory

```
cd {project_name}
```

### Installation

#### Creating the `rocketchip` Environment

Run the following command in your project directory. It will clone the conda environment with all dependencies needed in order to run the workflow outlined here. This creates an environment called `rocketchip`. If you would like to change the name, feel free to do so where the command says `rocketchip`. Please note that this may take quite a few minutes to run.

```
conda env create -f environment.yml --name rocketchip
```

Activate your environment using

```
conda activate rocketchip
```

Run everything downstream of this point in this conda environment. Note that you must activate this environment every time you restart your terminal.

#### Setting Up `rocketchip` Variables

You should do the following:

* Define `ROCKETCHIP_DATA`
* Define `ROCKETCHIP_SRC`
* Add `rocketchip` to your `PATH`

Modify your `.profile`, `.bash_profile`, `.bashrc`, `.zshrc`, or whatever your shell reads upon login with something like the following:

```
export ROCKETCHIP_DATA="/share/mylab/data/rocketchip
export ROCKETCHIP_SRC="/share/mylab/pkg/rocketchip
PATH="$PATH:$ROCKETCHIP_SRC"
```

(MAKE INSTRUCTIONS MORE CLEAR)

#### Installing Cisgenome (optional)

If you are using Cisgenome as your peak caller, you will need to install it [separately](https://www.biostat.jhsph.edu/~hji/cisgenome/index_files/download.htm), as it is not available through Conda. `rocketchip` is compatible with version 2.0.

To install it, enter the `tools/` directory or the directory of your choice (if so, just change the paths appropriately) and carry out the following commands.

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

## Running `rocketchip`

### How Does `rocketchip` Work?

`rocketchip` is an automated bioinformatics workflow that is capable of analyzing local ChIP-seq data or ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA), the largest publicly available sequence data database. `rocketchip` takes raw data inputs and generates the files required for data visualization and peak delineation.

![Workflow](https://github.com/vhaghani26/rocketchip/blob/main/rocketchip_flowchart.png)

### Making a Project File

In order to run `rocketchip`, you will need to create a project file. A template, `project_file.yaml`, is included in this repository. The contents should look like this:

```
Author: 
Project: 
Genome:
  Name: 
  Location: 
Reads:
  Samples:
    grp1: 
    grp2:
  Controls:
    ctl1: 
Readtype: 
Peaktype: 
Aligner: 
Deduplicator: 
Peakcaller: 
Threads: 
```

#### Author

List the authors associated with the project. This may be only you or include others as well. (Spaces allowed?)

#### Project

Write the name of your project. (Spaces allowed?)

#### Genome

| Organism  | Genome   | Link                                                                       |
| :-------: | :------: | :------------------------------------------------------------------------: |
| Fruitfly  | dm6      | https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz           |
| Human     | hg19     | https://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz     |
| Mouse     | mm10     | https://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz     |
| Rat       | rn6      | https://hgdownload.cse.ucsc.edu/goldenPath/rn6/bigZips/rn6.fa.gz           |
| Worm      | ce11     | https://hgdownload.soe.ucsc.edu/goldenPath/ce11/bigZips/chromFa.tar.gz     |
| Yeast     | sacCer3  | https://hgdownload.cse.ucsc.edu/goldenPath/sacCer3/bigZips/chromFa.tar.gz  |
| Zebrafish | danRer11 | https://hgdownload.cse.ucsc.edu/goldenPath/danRer11/bigZips/danRer11.fa.gz |

#### Name

#### Location

#### Reads

#### Samples

##### Groups

#### Controls

#### Readtype

#### Peaktype

#### Aligner

#### Deduplicator

#### Peakcaller

#### Threads

### Running `rocketchip` Locally

### Running `rocketchip` on SLURM (Recommended)

## Tutorial

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
Bigwig files are used for visualization of ChIP-seq data and are the final product of the analysis. Although bedGraphToBigWig can be used to convert the {sample}_treat_pileup.bdg files generated by MAS2 to bigwig format for visualization, deeptools was chosen to generate bigwig files because it allows for greater flexibility in the peak-calling software used if a user decides to change the peak-caller.

### 06_macs2_peaks
MACS2 is used to call peaks in the data. These peaks will be used in answering the biological question you are asking using the data. In many instances, the peaks correspond to the binding sites of a protein of interest.

(REWRITE AS NEEDED)

## Questions, Comments, or Concerns?

Feel free to contact me at [vhaghani@ucdavis.edu](vhaghani@ucdavis.edu), and I will get back to you as soon as possible.