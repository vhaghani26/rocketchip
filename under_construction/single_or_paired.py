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

# Add file name tags to distinguish paired and single end reads            
for sample in config["samples"]:
    if sample in paired_end:
        os.system(f'mv 01_raw_data/{sample}.html 01_raw_data/{sample}_paired.html')
    else:
        os.system(f'mv 01_raw_data/{sample}.html 01_raw_data/{sample}_single.html')
    os.system(f'touch 01_raw_data/{sample}_placeholder.txt')