#!/usr/bin/env python3

'''
Note: I ran this script in the functional_test/ directory:
    python3 functional_test.py
'''

####################
## Import Modules ##
####################

import sys
from collections import defaultdict
import pandas as pd
import os
import textwrap
import subprocess 

###########################
## Set Working Variables ##
###########################

# User-specific variables
authors = 'Viktoria_Haghani'

# Combinatorial testing variables
controltypes = ["with_control", "no_control"]  
readtypes = ["paired"]
peaktypes = ["narrow"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]
num_tests = 1

#########################
## Delineate Functions ##
#########################

# Make project files 
def generate_project_files(authors, controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Set up directory structure for project files if needed
    if not os.path.exists(f'project_files/'):
        print(f'Directory project_files/ not found. Creating project_files/')
        os.system(f'mkdir project_files')
    if not os.path.exists(f'project_files/with_control/'):
        print(f'Directory project_files/with_control/ not found. Creating project_files/with_control/')
        os.system(f'mkdir project_files/with_control')    
    if not os.path.exists(f'project_files/no_control/'):
        print(f'Directory project_files/no_control/ not found. Creating project_files/no_control/')
        os.system(f'mkdir project_files/no_control')  
    
    # Start combinatorial project file generation
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if control == "with_control":
                                    proj_file_info = textwrap.dedent(f"""
                                    Author: {authors}
                                    Project: functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                                    Genome:
                                        Name: genome
                                        Location: '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'genome.fa')}'
                                    Reads:
                                        Samples:
                                            grp1: 
                                                - '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'exp_a')}'
                                                - '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'exp_b')}'
                                        Controls:
                                            ctl1: 
                                                - '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'input')}'
                                    Readtype: {readtype}
                                    Peaktype: {peaktype}
                                    Aligner: {aligner}
                                    Deduplicator: {deduplicator}
                                    Peakcaller: {peakcaller}
                                    Threads: 1
                                    """)
                                elif control == "no_control":
                                    if peakcaller == "cisgenome" or peakcaller == "pepr": continue
                                    proj_file_info = textwrap.dedent(f"""
                                    Author: {authors}
                                    Project: functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                                    Genome:
                                        Name: genome
                                        Location: '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'genome.fa')}'
                                    Reads:
                                        Samples:
                                            grp1: 
                                            - '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'exp_a')}'
                                            - '{os.path.join(os.getcwd(), 'seq_data', f'{readtype}_{peaktype}', f'test_{i}', 'exp_b')}'
                                        Controls:
                                    Readtype: {readtype}
                                    Peaktype: {peaktype}
                                    Aligner: {aligner}
                                    Deduplicator: {deduplicator}
                                    Peakcaller: {peakcaller}
                                    Threads: 1
                                    """)
                                print(f'Generating project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml...')
                                os.system(f'touch project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')
                                with open(f'project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml', 'w') as f:
                                    f.write(f'{proj_file_info}')
                                os.system(f'sed -i \'1d\' project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml')

# Make Snakefiles
def generate_snakefiles(controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Set up directory structure for Snakefiles if needed
    if not os.path.exists(f'snakefiles/'):
        print(f'Directory snakefiles/ not found. Creating snakefiles/')
        os.system(f'mkdir snakefiles')
    if not os.path.exists(f'snakefiles/with_control/'):
        print(f'Directory snakefiles/with_control/ not found. Creating snakefiles/with_control/')
        os.system(f'mkdir snakefiles/with_control')    
    if not os.path.exists(f'snakefiles/no_control/'):
        print(f'Directory snakefiles/no_control/ not found. Creating snakefiles/no_control/')
        os.system(f'mkdir snakefiles/no_control') 

    # Start combinatorial Snakefile generation
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                else:
                                    # Make the directory structure if it does not already exist
                                    snakefile_dir = f'snakefiles/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
                                    if not os.path.exists(snakefile_dir):
                                        os.makedirs(snakefile_dir)
                                    
                                    # Create the snakefiles using Rocketchip
                                    print(f"Generating {snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}")
                                    os.system(f'rocketchip project_files/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}.yaml --data seq_data/{readtype}_{peaktype}/test_{i} --output_file {snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    
                                    # MACS3 has an issue with small read files requiring the --nomodel flag, so I will manually add it for the single-end data that are having problems with peak-calling
                                    if readtype == "single" and peakcaller == "macs3":
                                        print(f'Adding --nomodel flag in MACS3 for {snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                        file_to_open = f'{snakefile_dir}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}'
                                        # Read in the file
                                        with open(file_to_open, 'r') as file:
                                            filedata = file.read()
                                        # Replace the target string
                                        filedata = filedata.replace('macs3 callpeak ', 'macs3 callpeak --nomodel ')
                                        # Write the file out again
                                        with open(file_to_open, 'w') as file:
                                            file.write(filedata)

# Run Snakefiles
def run_snakefiles(controltypes, readtypes, peaktypes, aligners, peakcallers, deduplicators, num_tests):
    # Run Snakemake for all combinations
    for control in controltypes:
        for readtype in readtypes:
            for peaktype in peaktypes:
                for aligner in aligners:
                    for peakcaller in peakcallers:
                        for deduplicator in deduplicators:
                            for i in range(1, num_tests + 1):
                                if (control == "no_control") and (peakcaller == "cisgenome" or peakcaller == "pepr"):
                                    continue
                                else:
                                    # Change into snakefile directory
                                    os.chdir(f'snakefiles/{control}/functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    os.system('pwd')
                                    # Run snakefile
                                    os.system(f'snakemake -j 4 -s functional_test_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}_{control}')
                                    # Go back to original directory
                                    os.chdir(f'../../../')

    
####################
## Run Everything ##
####################

# Generate project_files
generate_project_files(authors = authors, controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests)

# Generate Snakefiles
generate_snakefiles(controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests)

# Run Snakefiles
run_snakefiles(controltypes = controltypes, readtypes = readtypes, peaktypes = peaktypes, aligners = aligners, peakcallers = peakcallers, deduplicators = deduplicators, num_tests = num_tests)