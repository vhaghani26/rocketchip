#!/usr/bin/env python3

import os
import yaml

paired_end = []
single_end = []
links = []

with open("samples.yaml", "r") as yamlfile:
    config = yaml.load(yamlfile)

for sample in config["samples"]:
    links.append(f'https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={sample}')
for link, sample in zip(links, config["samples"]):
    os.system(f'curl {link} > 01_raw_data/{sample}.html')
    with open(f"01_raw_data/{sample}.html", "r") as fp:
        readfile = fp.read()
        if '<td>PAIRED</td>' in readfile:
            paired_end.append(f'{sample}')
        else:
            single_end.append(f'{sample}')