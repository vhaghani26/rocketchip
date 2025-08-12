# Rocketchip: A Comprehensive Bioinformatics Workflow for ChIP-Seq Data Analysis

Rocketchip is an automated bioinformatics workflow that is capable of analyzing local ChIP-seq data or ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA), the largest publicly available sequence data database. `rocketchip` takes raw data inputs and generates the files required for data visualization and peak delineation.

![Workflow](https://github.com/vhaghani26/rocketchip/blob/main/rocketchip_flowchart.png)

## Table of Contents

* [Citation](#citation)
* [Installation](#installation)
    * [Dependencies](#dependencies)
 * [Running Rocketchip](#running-rocketchip)
	* [Data Storage](#data-storage)
  	* [Making a Project File](#making-a-project-file)
    * [Executing Rocketchip Locally](#executing-rocketchip-locally)
	* [Executing Rocketchip via SLURM Submission](#executing-rocketchip-via-slurm-submission)
* [Interpretting Outputs](#Interpretting-Outputs)
* [Installing Cisgenome](#installing-cisgenome)

## Citation

If you use Rocketchip, please cite it using the following:

Haghani V, Goyal A, Zhang A et al. Improving rigor and reproducibility in chromatin immunoprecipitation assay data analysis workflows with Rocketchip [version 1; peer review: 1 approved, 1 approved with reservations]. F1000Research 2025, 14:625 (https://doi.org/10.12688/f1000research.164319.1) 
 
## Installation

In order to install the `rocketchip` source code, please run:

```
pip install rocketchip
```

To confirm successful installation, confirm by running:

```
rocketchip --version
```

### Dependencies

Prior to installing the necessary dependencies, make sure you have Conda installed. Run the following command in your project directory. It will clone the Conda environment with all dependencies needed in order to run the workflow using all software options available. This creates an environment called `rocketchip` (you may change the name). Please note that this may take quite a few minutes to run.

```
conda env create -f environment.yml
```

Activate your environment using

```
conda activate rocketchip
```

Run everything downstream of this point in this Conda environment. Note that you must activate this environment every time you run the workflow. Please note that you can also modify the environment file as well. For instance, if you are only using `samtools` for deduplication, then you can remove options to install `sambamba` and `picard` to save space and make installation quicker.

## Running Rocketchip

### Data Storage

You have two options for managing data storate locations.

#### Option 1: Export `ROCKETCHIP_DATA`

**Purpose**: This method ensures that your raw data is consistently stored in a single location. It is particularly useful for labs that utilize public data or share the same raw data across multiple analyses.

**Benefits**

* Data is stored in a designated location for all analyses using Rocketchip
* Files are aliased into the project directory, preventing duplication. This is advantageous if you are using the same genome or sample data for multiple analyses, as it saves both time and storage

**How to Set It Up**

Add the `ROCKETCHIP_DATA` variable to your configuration file (e.g. `.bashrc`, `.profile`, or another file appropriate for your file system) and set the path to your desired data storage location. For instance:

```
export ROCKETCHIP_DATA="/shared_drive/your_lab/raw_data/"
```

#### Option 2: Use the `--data` Argument

**Purpose**: This method allows you to specify the data directory directly when running Rocketchip from the command line.

**Benefits**

* Easier to implement, especially for individual/one-time analyses
* Allows for personalized organization 

**How to Set It Up**

You run it during the `rocketchip` command like so for the current working directory:

```
rocketchip --data .
```

Alternatively, you can specify another path:

```
rocketchip --data /shared_drive/your_lab/raw_data/
```

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
* Molecule (optional) - if you are interested in running Rocketchip for RIP-seq, then you should specify "RNA" here (e.g. `Molecule: RNA`), otherwise you can say "DNA" or leave it blank (it defaults to DNA)
	* Note that if you are using Rocketchip with RNA, it is recommended to use the STAR aligner

Note: For local read data, use the absolute paths as entries. For example, the read path can be `absolute/path/to/my/sample_id`. Do not put the `fq.gz` extension or the read direction (forward vs. reverse) in the path to your read names. Just end the path at the sample name. Based on if you tell Rocketchip whether the data is single- or paired-end, it will match the appropriate files automatically. For read data from the SRA, use the SRA ID for the sample instead of a path. This goes for both the samples and controls.

Here are various examples of project yaml files:
* [One sample with one replicate, no control](https://github.com/vhaghani26/rocketchip_tests/blob/main/run_times/human/human_SRR17409984.yaml)
* [Two samples with three replicates each, no control](https://github.com/vhaghani26/rocketchip_tests/blob/main/cut_and_tag_run/cut_and_tag.yaml)
* [Three samples with two replicates each, no control](https://github.com/vhaghani26/rocketchip_tests/blob/main/cut_and_tag_run/cut_and_run_se.yaml)
* [One sample with two replicates, one control with one replicate](https://github.com/vhaghani26/rocketchip_tests/blob/main/cut_and_tag_run/cut_and_run_pe.yaml)
* [One sample with two replicates, one control with two replicates](https://github.com/vhaghani26/rocketchip_tests/blob/main/replicability/replicability.yaml)

### Executing Rocketchip Locally

1. **Run Rocketchip**

Enter the directory containing the `project_file.yaml` that you have set up (you can rename this, just make sure to change the name in the command below). Assuming you have set `ROCKETCHIP_DATA`, all you need to do is run the following:

```
rocketchip project_file.yaml --output_file {output_file_name}
```

If you have not set `ROCKETCHIP_DATA`, you will need to specify it at the command line:

```
rocketchip project_file.yaml --output_file {output_file_name} --data {directory_to_store_the_data}
```

This will generate the Snakefile you have named `{output_file_name}` that we will run in the next step.

2. **Run Snakemake**

Now, you will run Snakemake. This follows Snakemake's command line usage, but at it's simplest, you can run:

```
snakemake -j 1 -s {output_file_name}
```

Increase `-j` to match the number of jobs you would like to parallelize.

### Executing Rocketchip via SLURM Submission

Running Rocketchip via SLURM provides several advantages over running it locally:

- **Resource Management**: SLURM allows for individual job submissions, helping to allocate resources effectively without overwhelming your local machine
- **Job Scheduling**: It manages job scheduling and monitoring, reducing the risk of disconnections and failures
- **Efficiency**: By optimizing resource allocation, SLURM ensures smoother and more efficient execution of Rocketchip
- **Reduced Failures**: If a job fails, SLURM can handle retries and manage job states, minimizing the need for manual intervention

Overall, using SLURM enhances the reliability and performance of running Rocketchip, especially when dealing with resource-intensive tasks and large sample numbers.

1. **Run Rocketchip**

This step is not resource intensive, so you may run it by doing the following. Enter the directory containing the `project_file.yaml` that you have set up (you can rename this, just make sure to change the name in the command below). Assuming you have set `ROCKETCHIP_DATA`, run:

```
rocketchip project_file.yaml --output_file {output_file_name}
```

If you have not set `ROCKETCHIP_DATA`, you will need to specify it at the command line:

```
rocketchip project_file.yaml --output_file {output_file_name} --data {directory_to_store_the_data}
```

This will generate the Snakefile you have named `{output_file_name}` that we will run in the next step.

2. **Download SLURM Configuration Files**

Navigate to the [SLURM Directory](https://github.com/vhaghani26/rocketchip/tree/main/slurm) and download `config.yaml` and `slurm-status.py`. Store them in a directory called `slurm/`. 

3. **Run Snakemake**

Snakemake manages the submission of jobs to SLURM, so wherever you run it, it will need to stay open. As such, I recommend running it in [screen](https://github.com/vhaghani26/python_focus_group/blob/main/Session_14/session_14.md) or something similar, like [tmux](https://github.com/tmux/tmux). Activate your conda environment. When you are ready, run Snakemake:

```
snakemake -s {output_file_name} --profile slurm/
```

This assumes you stored `slurm/` in the current working directory. Just make sure you put the proper path, whether absolute or relative, to the `slurm/` directory that contains `config.yaml` and `slurm-status.py`. `config.yaml` handles resource allocations and job submission details, while `slurm-status.py` is simply to check on the status of the jobs submitted to SLURM in order to track if jobs are running, completed, or failed.

#### Advanced Configuration

If you are interested in digging more deeply into optimizing your workflow for your samples, including changing resource allocations, configuring your SLURM profile submission differently, etc., refer to the [Snakemake SLURM Documentation](https://snakemake.readthedocs.io/en/v7.19.1/executing/cluster.html) or feel free to reach out to me if the documentation does not answer the questions you have.

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

To install it, you can carry out the following commands.

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
export PATH=$PATH:{your_directory}/cisgenome_project/bin
```

Now, either source your configuration file or restart your terminal.  To confirm proper installation, you can run:

```
seqpeak
```

This will display the options available for use, indicating that you are able to execute it.

## Questions, Comments, or Concerns?

Feel free to contact me at [vhaghani@ucdavis.edu](vhaghani@ucdavis.edu), and I will get back to you as soon as possible.
