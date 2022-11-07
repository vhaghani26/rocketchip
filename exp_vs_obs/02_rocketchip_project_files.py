#!/usr/bin/env python3

import os

readtypes = ["paired", "single"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "no_deduplication", "sambamba", "picard"]

for readtype in readtypes:
    for peaktype in peaktypes:
        for aligner in aligners:
            for peakcaller in peakcallers:
                for deduplicator in deduplicators:
                    for i in range(1,7):
                        os.system(f'rocketchip project_files/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}.yaml --output_file snakefiles/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}')

# rocketchip exp_vs_obs/project_files/exp_vs_obs_paired_narrow_bwa_mem_macs3_samtools_test_1.yaml --data exp_vs_obs/seq_data/paired_narrow/test_1 --src . --output_file exp_vs_obs/snakefiles/exp_vs_obs_paired_narrow_bwa_mem_macs3_samtools_test_1