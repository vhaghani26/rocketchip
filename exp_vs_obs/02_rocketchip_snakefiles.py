#!/usr/bin/env python3

'''
Note: I ran this script in the main rocketchip directory using:
    python3 exp_vs_obs/02_rocketchip_snakefiles.py
'''

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
                        if not os.path.exists(f'exp_vs_obs/snakefiles/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}'):
                            os.system(f'mkdir exp_vs_obs/snakefiles/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}')
                        os.system(f'python3 rocketchip exp_vs_obs/project_files/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}.yaml --data exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i} --src . --output_file exp_vs_obs/snakefiles/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}')