# Rocketchip: A Comprehensive Bioinformatics Workflow for ChIP-Seq Data Analysis

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


#######################################################################################################


































### Conda Preparation
If you already have conda or an environment manager, continue to the "Software Installation" section. If not, the following information may be helpful to you.

Briefly, make sure conda is installed. Once conda is installed, it needs to be initialized. To initialize conda, use
```
conda init
```
To clean up your shell settings for a cleaner prompt, you can also use
```
echo "PS1='\w $ '" >> .bashrc
```

Close your terminal and then open a new one to make sure initialization is fully carried out.

Now you will set up your channels. Paste the following commands into your terminal.
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Software Installation
The only software necessary for running this workflow is Snakemake. In this tutorial, we will create an environment called `chipseq`. To create the environment and install Snakemake, paste the following into your terminal. Click `y` if/when prompted.
```
conda create -y --name chipseq -c conda-forge -c bioconda fastqc
conda activate chipseq
conda install -c bioconda -c conda-forge snakemake-minimal
```

This completes the software installation stage. Note that **every time** you are using snakemake, you will need to activate your environment using `conda activate chipseq`.

### Download the Repository
You will want to download the repository, as it contains the necessary software files and Snakefile needed for the analysis. In this tutorial, we will call the directory we are working in `chipseq_tutorial`. To do so, run
```
git clone https://github.com/vhaghani26/rocketchip.git chipseq_tutorial
```
This clones the repository and its contents into a directory called `chipseq_tutorial`. You should now see the contents outlined in the "Repository Contents" section outlined above.

### Edit the YAML File
Open the `paired_samples.yml` file and add each SRA ID you are interested in visualizing on its own line using the following format:
```
---
samples:
   - SRAID1
   - SRAID2
```
This file will be used as a configuration file to pipe in the samples used in the analysis. The syntax in this file is important, so be aware that errors relating to `wildcards` are likely caused by the entries in this file. Also note that the SRA ID is NOT the BioSample Accession Number, which starts with the "SAMN" prefix. The SRA ID begins with the "SRR" prefix.

In this tutorial, we will use the following SRA IDs (paired-end data):
```
---
samples:
   - SRR5785190
   - SRR6760417
```
Now save and close the YAML file.

### Running the Workflow Locally
To run the workflow locally, all you need to do now is paste this command into your terminal:
```
snakemake -j 4 -p --use-conda -s paired_end_snakefile
```
It is important to note that because we are working with "big data," this will take upwards of 24 hours to run from scratch. Although it can be done provided your internet connection allows, an alternative is to run it rule by rule or to use other methods, such as `screen` to ensure completion of the job. Furthermore, we are working with paired-end data in this tutorial. If you are working with single-end data (which can be verified by viewing the data entry on the SRA database), please use `-s single_end_snakefile` instead.

To run it rule by rule, carry out the following commands in this order:
```
snakemake -j 4 -p --use-conda -s paired_end_snakefile make_directories
snakemake -j 4 -p --use-conda -s paired_end_snakefile download_data
snakemake -j 4 -p --use-conda -s paired_end_snakefile sra_to_fastq
snakemake -j 4 -p --use-conda -s paired_end_snakefile fastqc_precheck
snakemake -j 4 -p --use-conda -s paired_end_snakefile download_genome
snakemake -j 4 -p --use-conda -s paired_end_snakefile process_genome
snakemake -j 4 -p --use-conda -s paired_end_snakefile set_alignment_reference
snakemake -j 4 -p --use-conda -s paired_end_snakefile align_reads
snakemake -j 4 -p --use-conda -s paired_end_snakefile sam_to_bam
snakemake -j 4 -p --use-conda -s paired_end_snakefile sam_fixmate
snakemake -j 4 -p --use-conda -s paired_end_snakefile sam_sort
snakemake -j 4 -p --use-conda -s paired_end_snakefile sam_markdup
snakemake -j 4 -p --use-conda -s paired_end_snakefile sam_index
snakemake -j 4 -p --use-conda -s paired_end_snakefile bam_to_bigwig
snakemake -j 4 -p --use-conda -s paired_end_snakefile call_peaks
snakemake -j 4 -p --use-conda -s paired_end_snakefile fastqc_postprocessing
```
Some rules, such as `set_alignment_reference` and `align_reads` take much longer to run than the rest of the workflow. The majority of the rules, however, will take no more than 20 minutes each to run using only the two samples provided in this tutorial. With each extra sample included, the time required increases. It may be helpful to keep the genome alignment reference files for future analysis since that's the most resource- and time-intensive step.

### Running the Workflow on the Cluster or on an HPC
It is most highly recommended to run the workflow through SLURM. Based on how conda was/is initialized and how/where your environment is set up, there will be things you should change in the SLURM script.

If you are new to using conda or SLURM, it may be helpful to submit a test script that requests minimal resources to make sure conda works correctly. This is a SLURM template script you can use to make sure conda and SLURM are working for you before you submit a big resource- and time-intensive job like this workflow. Make sure to change the information necessary, namely the information at the top preceeded by `#SBATCH` and the three conda-related lines.
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
#SBATCH --chdir=/home/vhaghani/chipseq_tutorial                # Directory I want the job to run in

# Source .profile or .bashrc so conda can be used (file depends on where you have your conda stuff)
source ~/.profile

# Initialize conda
. /software/anaconda3/4.8.3/lssc0-linux/etc/profile.d/conda.sh

# Activate your desired conda environment
conda activate chipseq

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
#SBATCH -J chipseq                                             # Name for job
#SBATCH -o chipseq_slurm.j%j.out                               # File to write STDOUT to
#SBATCH -e chipseq_slurm.j%j.err                               # File to write error output to
#SBATCH -N 1                                                   # Number of nodes/computers
#SBATCH -n 16                                                  # Number of cores
#SBATCH -t 48:00:00                                            # Ask for no more than 48 hours
#SBATCH --mem=16gb                                             # Ask for no more than 16 GB of memory
#SBATCH --chdir=/home/vhaghani/chipseq_tutorial                # Directory I want the job to run in

# Source .profile or .bashrc so conda can be used (file depends on where you have your conda stuff)
source ~/.profile

# Initialize conda
. /software/anaconda3/4.8.3/lssc0-linux/etc/profile.d/conda.sh

# Activate your desired conda environment
conda activate chipseq

# Fail on weird errors
set -o nounset
set -o errexit
set -x

# Run the snakemake!
snakemake -j 4 -p --use-conda -s [endedness]

# Print out various information about the job
env | grep SLURM                                               # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}                              # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch
```

Once the job is complete, you should have all the directories and outputs generated by the workflow. Take some time to explore and see what there is. Now you are ready to work with your own data files. Congratulations and good luck!

## Updates

Current version: 2.0.0

### Launch Phase (v1.0.0) - Deprecated

This version allowed for the analysis of only paired-end ChIP-seq data specific to mouse experiments.

### Cruise Phase (v2.0.0) - Available

Version 2.0.0, the currently available version of Rocketchip, allows for the use of single- and paired-end ChIP-seq data analysis.

### Encounter Phase (v3.0.0) - Under Construction
This update will include the ability to use different reference genomes for alignment.

Updated software description:
Rocketchip is an automated bioinformatics workflow that downloads and processes both single- and paired-end ChIP-seq data from the National Center for Biotechnology Information (NCBI) Sequence Read Archive (SRA) and uses it to generate the files required for data visualization and peak delineation. The ChIP-seq data is downloaded directly from the SRA using the SRA Toolkit (sra-tools, v2.10.9). The raw sequence read file is then either split into its respective paired-end read files by the SRA Toolkit and converted into FASTQ file format or maintained as a single file and converted to a FASTQ file. Reads are aligned to the reference genome assembly of the user's choice using the Burrows-Wheeler Algorithm (BWA) Maximal Exact Match (MEM) software (bwa, v0.7.17). Samtools (samtools, v1.12) is used for file format conversion and deduplication of sequence data. Deeptools (deeptools, v3.5.1) is used to convert data to the bigwig file format, which can be used for visualization of ChIP-seq data in the UCSC Genome Browser or other visualization tools. Peaks are called using MACS2 (macs2, v2.2.7.1). FastQC (fastqc, v0.11.9) carries out a sequence quality control analysis pre-processing, on the raw sequence data, and post-processing, after sequence alignment takes place. Overall, this analysis carries out the major steps of ChIP-seq data analysis and generates output files to be used as figures and to be used in further analysis.

### Extended Operations Phase

Future ideas for Rocketchip include the incorporation of track line generators for narrowPeak and bigwig files to visualize ChIP-seq peaks and read distribution, respectively, on the UCSC Genome Browser.

Feel free to contact me at [vhaghani@ucdavis.edu](vhaghani@ucdavis.edu) if you have any questions, comments, or concerns.


