#!/usr/bin/env python3

import os
import yaml

paired_end = []
single_end = []
links = []

# Load in sample data from configuration file
with open("samples.yaml", "r") as yamlfile:
    config = yaml.load(yamlfile, Loader = yaml.FullLoader)

# Generate links corresponding to each sample
for sample in config["samples"]:
    links.append(f'https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={sample}')
    
# Download sample HTML files
for link, sample in zip(links, config["samples"]):
    os.system(f'curl {link} > 01_raw_data/{sample}.html')
# Determine if sample data is single or paired end
    with open(f"01_raw_data/{sample}.html", "r") as fp:
        readfile = fp.read()
        if '<td>PAIRED</td>' in readfile:
            paired_end.append(f'{sample}')
        else:
            single_end.append(f'{sample}')

# Check if directories exist and create them if they don't
# The directories are made here to avoid the error "mv: cannot move ...: Not a directory" 
directories = ["01_raw_data/html_files/",
               "01_raw_data/html_files/single/",
               "01_raw_data/html_files/paired/",
               "01_raw_data/sra_files/",
               "01_raw_data/sra_files/single/",
               "01_raw_data/sra_files/paired/",
               "01_raw_data/sra_files/",
               "01_raw_data/sra_files/single/",
               "01_raw_data/sra_files/paired/"
               ]
for directory in directories:
    isExist = os.path.exists(directory)
    if isExist == True:
        print(f'The directory {directory} already exists and will not be overridden')
    else:
        os.system(f'mkdir {directory}')
        print(f'The directory {directory} has been created')
    
# Add file name tags to distinguish paired and single end reads     
for sample in config["samples"]:
    if sample in paired_end:
        os.system(f'mv 01_raw_data/{sample}.html 01_raw_data/html_files/paired/{sample}_paired.html')
        print(f'{sample} is a paired-end read')
    else:
        os.system(f'mv 01_raw_data/{sample}.html 01_raw_data/html_files/single/{sample}_single.html')
        print(f'{sample} is a single-end read')
    os.system(f'touch 01_raw_data/{sample}_placeholder.txt')