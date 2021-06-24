#!/bin/bash
#
#SBATCH --mail-user=vhaghani@ucdavis.edu                       # User email to receive updates
#SBATCH --mail-type=ALL                                        # Get an email when the job begins, ends, or if it fails
#SBATCH -p production                                          # Partition, or queue, to assign to
#SBATCH -J mechip2                                             # Name for job
#SBATCH -o mechip2_slurm.j%j.out                               # File to write STDOUT to
#SBATCH -e mechip2_slurm.j%j.err                               # File to write error output to
#SBATCH --nodes=1                                              # Number of nodes/computers
#SBATCH --ntasks=32                                            # Number of cores
#SBATCH -c 8                                                   # Eight cores per task
#SBATCH -t 48:00:00                                            # Ask for no more than 48 hours
#SBATCH --mem=16gb                                             # Ask for no more than 16 GB of memory
#SBATCH --chdir=/share/korflab/home/viki/mechip2               # Directory I want the job to run in

# Source profile so conda can be used
source /share/korflab/home/viki/profile

# Initialize conda
. /software/anaconda3/4.8.3/lssc0-linux/etc/profile.d/conda.sh

# Activate your desired conda environment
conda activate /share/korflab/home/viki/.conda/mechip2

# Fail on weird errors
set -o nounset
set -o errexit
set -x

# Run the snakemake!
snakemake -j 4 -p --use-conda

# Print out various information about the job
env | grep SLURM                                               # Print out values of the current jobs SLURM environment variables
  
scontrol show job ${SLURM_JOB_ID}                              # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU' -P ${SLURM_JOB_ID}.batch

# Note: Run dos2unix {filename} if sbatch DOS line break error occurs