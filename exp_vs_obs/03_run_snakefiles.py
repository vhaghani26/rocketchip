#!/usr/bin/env python3

import os

'''
Note: I ran this script in the main rocketchip directory using:
    python3 exp_vs_obs/03_run_snakefiles.py
'''

readtypes = ["paired"]#, "single"]
peaktypes = ["broad"] #["narrow", "broad"]
aligners = ["bowtie2"] #["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["cisgenome"] #["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["no_deduplication"] #["samtools", "no_deduplication", "sambamba", "picard"]

for readtype in readtypes:
    for peaktype in peaktypes:
        for aligner in aligners:
            for peakcaller in peakcallers:
                for deduplicator in deduplicators:
                    for i in range(1,7):
                        # Dry run first
                        os.system(f'snakemake -j 4 -s exp_vs_obs/snakefiles/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}')