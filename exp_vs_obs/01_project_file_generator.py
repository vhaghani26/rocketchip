#!/usr/bin/env python3

'''
Note: I ran this script in the exp_vs_obs/ directory:
    python3 01_project_file_generator.py
'''

import textwrap
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
                        proj_file_info = textwrap.dedent(f"""
                        Author: Viktoria_Haghani_and_Aditi_Goyal_and_Alan_Zhang
                        Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                        Genome:
                          Name: genome
                          Location: '/mnt/c/Users/Viki Haghani/Documents/Github/rocketchip/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/genome.fa'
                        Reads:
                          Samples:
                            grp1: 
                              - '/mnt/c/Users/Viki Haghani/Documents/Github/rocketchip/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_a'
                              - '/mnt/c/Users/Viki Haghani/Documents/Github/rocketchip/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/exp_b'
                          Controls:
                            ctl1: 
                              - '/mnt/c/Users/Viki Haghani/Documents/Github/rocketchip/exp_vs_obs/seq_data/{readtype}_{peaktype}/test_{i}/input'
                        Readtype: {readtype}
                        Peaktype: {peaktype}
                        Aligner: {aligner}
                        Deduplicator: {deduplicator}
                        Peakcaller: {peakcaller}
                        Threads: 1
                        """)
                        os.system(f'touch project_files/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}.yaml')
                        with open(f'project_files/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}.yaml', 'w') as f:
                            f.write(f'{proj_file_info}')
                        os.system(f'sed -i \'1d\' project_files/exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}_test_{i}.yaml')
