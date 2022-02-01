# Rocketchip: A Comprehensive Bioinformatics Workflow for ChIP-Seq Data Analysis

vT.B.D.

**Rocketchip is currently being redone. The usage information listed below is no longer relevant/compatible with the current version of Rocketchip. The README.md will be redone shortly to include updated usage instructions and a tutorial for getting started.**

## What is Rocketchip?

Rocketchip is an automated bioinformatics workflow that downloads ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) and uses it to generate the files required for data visualization and peak delineation. It was created using the workflow manager, Snakemake. Both single- and paired-end data can be used, but they must be run separately. ChIP-seq data is downloaded directly from the SRA using the SRA Toolkit (sra-tools, v2.10.9). The raw sequence read file is then split into its respective paired-end read files by the SRA Toolkit (if the data is paired-end) and converted into FASTQ file format (sra-tools, v2.10.9). Reads are aligned to the genome of the user's choice from the UCSC Genome Browser using the Burrows-Wheeler Algorithm (BWA) Maximal Exact Match (MEM) software (bwa, v0.7.17). Samtools (samtools, v1.12) is used for file format conversion and deduplication of sequence data. Deeptools (deeptools, v3.5.1) is used to convert data to the bigwig file format, which can be used for visualization of ChIP-seq data in the UCSC Genome Browser or other visualization tools. Peaks are called using MACS2 (macs2, v2.2.7.1). FastQC (fastqc, v0.11.9) carries out a sequence quality control analysis pre-processing, on the raw sequence data, and post-processing, after sequence alignment takes place. Overall, this workflow carries out the major steps of ChIP-seq data analysis and generates output files to be used as figures and to be used in further analysis.

## Getting Started

In order to use Rocketchip, clone [this GitHub repository](https://github.com/vhaghani26/rocketchip.git) into the directory of your choice. It is highly recommended to try the provided tutorial to better understand specific usage of Rocketchip. 

## Tutorial

This is a small demo that shows you how to setup and run rocketchip analyses.

### Setup

1. Install conda if not already installed
2. Clone the rocketchip repository
3. Create the rocketchip conda enviornment

Download a conda installer such as mini-conda and install it with a command that looks something like the following:

```
sh Miniconda3-py39_4.10.3-Linux-x86_64.sh
```

Also install `mamba` which is better way to run conda installs

```
conda install mamba -n base -c conda-forge
```

Clone the rocketchip repository to wherever you typically keep repositories. For example this might be your home directory.

```
cd $HOME
git clone https://github.com/vhaghani26/rocketchip
```

Create the conda environment for rocketchip. This will take a few minutes to run.

```
cd rocketchip
mamba env create --file rocketchip.yaml
conda activate rocketchip
```

All of the commands here in **Setup** only need to be done once.

### Testing

Let's test rocketchip to make sure the installation is working. Create a new directory, possibly in your home directory for testing. We will do all of our testing in this directory. Once complete, you can remove it.

```
cd $HOME
mkdir test
cd test
```

Make sure your rocketchip conda environment is active. If not, activate it.

```
conda activate rocketchip
```

The rocketchip program requires three parameters `--genome`, `--sra`, and `--project`. We will also be using the `--data` and `--src` parameters, which you can omit later after testing (see below).

#### --genome

The first time you execute rocketchip with a particular genome, it will download and index the genome. This step happens only once. Subsequent analyses will use the local files. This test uses S. cereviseae (sacCer3) to minimize 
dowload and processing time.

#### --sra

The first time you execute rocketchip with a list of SRA identifiers, it will download them and extract the sequences. This step happens only once per SRA identifier. Subsequent analyses using previously called SRA identifiers will use the local files that have already been downloaded. The SRA files used in this test are small to minimize download and processing time.

* SRR12926698 (paired-end reads)
* SRR9257200 (single-end reads)

#### --project

Each analysis is stored in a project directory. In this test, we will be creating a couple projects: a paired-end layout wth narrow peak calling and a single-end layout with broad peak calling.

#### --data

Local genome and SRA files must be stored _somewhere_. In this test, we will put them in a directory called `cache`, inside our `test` directory. If you have followed the directions, your terminal is already in the `test` directory.

```
mkdir cache
```

#### Test 1: Paired-End, Narrow Peaks

Run the following command. This assumes your Rocketchip source and test directories are both located in the same directory (e.g. your home directory). It should take about 1 minute to download and process the genome and SRA files (assuming typical network speeds and CPUs).

```
../rocketchip/rocketchip --data cache --genome sacCer3 --src ../rocketchip --project demo1 --sra 
```

To start the analysis, change to the demo directory and run `snakemake`. This should take only a few minutes to run and uses minimal resources.

```
chdir demo
snakemake
```

Look in the XXX directory and examine the YYY file...

#### Test 2: Single-End, Broad Peaks

The following command is similar to the first test, but the SRA file comes from single-end sequencing rather than paired.

```
../rocketchip/rocketchip --data cache --genome sacCer3 --src ../rocketchip --project demo2 --sra SRR9257200 --broad
```

Look in the XXX directory and examine the YYY file...

## Post Demo Refinements

You can now remove the `test` directory.

To make subsequent analyses easier, you should do the following:

* define `ROCKETCHIP_DATA`
* define `ROCKETCHIP_SRC`
* add rocketchip to your `PATH`

Where do you want to put your local copies of genome and fastq files? This should be a shared location where multiple projects can reuse the same genome and fastq files so that you don't have to download them multiple times.

The rocketchip source directory can also be shared. Modify your `.profile`, `.bashrc`, `.zshrc`, or whatever your shell reads upon login with something like the following:

```
export ROCKETCHIP_DATA="/share/mylab/data/rocketchip
export ROCKETCHIP_SRC="/share/mylab/pkg/rocketchip
PATH="$PATH:$ROCKETCHIP_SRC"
```

## Rocketchip Outputs
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

## Repository Contents

### `snakefiles/` 

`snakefiles/` contains the different snakefiles corresponding to the type of analysis being run. Depending on the user's data, namely whether the analysis is for single-end reads, paired-end reads, narrow peak analysis, or broad peak analysis, Rocketchip will run the appropriate workflow. The user should not edit the files within `snakefiles/`.

### .gitignore
The `.gitignore` file contains the directories and consequently the data generated by each rule. Because the data is not pushed to git, it is strongly recommended to back up your data somewhere that it is protected.

## FAQs

### Can I use my own data instead of data from the SRA?

The tutorial above outlined the process for using publicly available ChIP-seq data hosted on the NCBI SRA. However (make notes about using personal data if possible)

### Can I run Rocketchip on the Cluster or an HPC?

Yes! In fact, with larger genomes, it is highly recommended to run the workflow using SLURM. If you are new to using conda or SLURM, it may be helpful to submit a test script that requests minimal resources to make sure conda works correctly. This is a SLURM template script you can use to make sure conda and SLURM are working for you before you submit a big resource- and time-intensive job. Make sure to change the information necessary, namely the information at the top preceeded by `#SBATCH` and the three conda-related lines.

```
#!/bin/bash
#
#SBATCH --mail-user=vhaghani@ucdavis.edu                       # User email to receive updates
#SBATCH --mail-type=ALL                                        # Get an email when the job begins, ends, or if it fails
#SBATCH -p production                                          # Partition, or queue, to assign to
#SBATCH -J test                                                # Name for job
#SBATCH -o test_slurm.j%j.out                                  # File to write STDOUT to
#SBATCH -e test_slurm.j%j.err                                  # File to write error output to
#SBATCH -N 1                                                   # Number of nodes/computers
#SBATCH -n 1                                                   # Number of cores
#SBATCH -c 1                                                   # Eight cores per task
#SBATCH -t 00:05:00                                            # Ask for no more than 5 minutes
#SBATCH --mem=1gb                                              # Ask for no more than 1 GB of memory
#SBATCH --chdir=/home/vhaghani/test_directory                  # Directory I want the job to run in

# Source .profile, .bashrc, or .zshrc so conda can be used (file depends on where you have your conda stuff)
source ~/.profile

# Initialize conda (path depends on where you have your conda stuff)
. /software/anaconda3/4.8.3/lssc0-linux/etc/profile.d/conda.sh

# Activate the Rocketchip environment
conda activate rocketchip

# Fail on weird errors
set -o nounset
set -o errexit
set -x

# Sample job
echo Hello World > helloworld.txt

# Print out various information about the job
env | grep SLURM                                               # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}                              # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
```

To submit the script, use the command `sbatch example.slurm`. Note: Run `dos2unix example.slurm` if you get an sbatch DOS line break error and resubmit.

Once you ensure that your job is running correctly, you can switch the resources and commands to be submitted so that the file looks more like this. Note that there is a big differences in the resources and time being requested.

```
#!/bin/bash
#
#SBATCH --mail-user=vhaghani@ucdavis.edu                       # User email to receive updates
#SBATCH --mail-type=ALL                                        # Get an email when the job begins, ends, or if it fails
#SBATCH -p production                                          # Partition, or queue, to assign to
#SBATCH -J rocketchip_project                                  # Name for job
#SBATCH -o rocketchip_slurm.j%j.out                            # File to write STDOUT to
#SBATCH -e rocketcihp_slurm.j%j.err                            # File to write error output to
#SBATCH -N 1                                                   # Number of nodes/computers
#SBATCH -n 16                                                  # Number of cores
#SBATCH -t 48:00:00                                            # Ask for no more than 48 hours
#SBATCH --mem=16gb                                             # Ask for no more than 16 GB of memory
#SBATCH --chdir=/home/vhaghani/rocketchip_project              # Directory I want the job to run in

# Source .profile, .bashrc, or .zshrc so conda can be used (file depends on where you have your conda stuff)
source ~/.profile

# Initialize conda (path depends on where you have your conda stuff)
. /software/anaconda3/4.8.3/lssc0-linux/etc/profile.d/conda.sh

# Activate the Rocketchip environment
conda activate rocketchip

# Fail on weird errors
set -o nounset
set -o errexit
set -x

# Run Rocketchip (note that this command is from the tutorial and will need to be changed to match your setup)
../rocketchip/rocketchip --data cache --genome sacCer3 --src ../rocketchip --project demo2 --sra SRR9257200 --broad

# Print out various information about the job
env | grep SLURM                                               # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}                              # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
```

Once the job is complete, you should have all the directories and outputs generated by the workflow. Take some time to explore and see what there is. Now you are ready to work with your own data files. Congratulations and good luck!


## Questions, Comments, or Concerns?

Feel free to contact me at [vhaghani@ucdavis.edu](vhaghani@ucdavis.edu), and I will get back to you as soon as possible.