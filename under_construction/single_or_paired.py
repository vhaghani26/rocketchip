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
        os.system(f'touch {sample}_alignment_snakefile')
        print(f'{sample} is a paired-end read')
    else:
        os.system(f'mv 01_raw_data/{sample}.html 01_raw_data/{sample}_single.html')
        os.system(f'touch {sample}_alignment_snakefile')
        print(f'{sample} is a single-end read')
    with open(f'{sample}_alignment_snakefile', 'w') as fp:
        if sample in paired_end:
            paired_lines = [
                'rule sra_to_fastq_paired:\n',
                '    message: "Converting SRA file to FASTQ file format"\n',
                '    conda: "00_conda_software/chip_sra.yml"\n',
                f'    input: "01_raw_data/{sample}/{sample}.sra"\n',
                '    output:\n',
                f'        "01_raw_data/{sample}_1.fastq.gz",\n',
                f'        "01_raw_data/{sample}_2.fastq.gz"\n',
                f'    log: "00_logs/{sample}_sra_to_fastq.log"\n',
                '    shell: "fastq-dump {input} --split-files --gzip --outdir 01_raw_data/ 2> {log}"\n',
                '\n',
                'rule fastqc_precheck:\n',
                '    message: "Running quality control on samples pre-processing"\n',
                '    conda: "00_conda_software/chip_fastqc.yml"\n',
                '    input:\n',
                f'        r1 = "01_raw_data/{sample}_1.fastq.gz",\n',
                f'        r2 = "01_raw_data/{sample}_2.fastq.gz",\n',
                f'        dependency = "01_raw_data/{sample}_paired.html"\n',
                '    output:\n',
                f'        "02_fastqc_analysis/{sample}_1_fastqc.html",\n',
                f'        "02_fastqc_analysis/{sample}_1_fastqc.zip",\n',
                f'        "02_fastqc_analysis/{sample}_2_fastqc.html",\n',
                f'        "02_fastqc_analysis/{sample}_2_fastqc.zip"\n',
                '    log:\n',
                f'        r1 = "00_logs/{sample}_fastqc_precheck_r1.log",\n',
                f'        r2 = "00_logs/{sample}_fastqc_precheck_r2.log"\n',
                '    shell: """\n',
                '    fastqc {input.r1} --outdir 02_fastqc_analysis/ 2> {log.r1}\n',
                '    fastqc {input.r2} --outdir 02_fastqc_analysis/ 2> {log.r2}\n',
                '    """\n',
                '\n',
                'rule align_reads:\n',
                '    message: "Aligning paired end reads to GRCm39/mm39 reference genome"\n',
                '    conda: "00_conda_software/chip_bwa.yml"\n',
                '    input:\n',
                f'        r1 = "01_raw_data/{sample}_1.fastq.gz",\n',
                f'        r2 = "01_raw_data/{sample}_2.fastq.gz",\n',
                '        genome = multiext("01_raw_data/mm39", ".amb", ".ann", ".bwt", ".pac", ".sa")\n',
                f'    output: "03_sam_files/{sample}.sam"\n',
                f'    log: "00_logs/{sample}_align_reads_err.log"\n',
                '    shell: "bwa mem 01_raw_data/mm39 {input.r1} {input.r2} > {output} 2> {log}"\n',
            ]
            fp.writelines(paired_lines)
        else:
            print(f'{sample} is a single-end read')