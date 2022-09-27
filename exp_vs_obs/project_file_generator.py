#!/usr/bin/env python3

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
                        if readtype == "paired":
                            read_data_path_1 = f"/seq_data/{readtype}_{peaktype}/test_{i}/exp_a_1.fastq.gz"
                            read_data_path_2 = f"/seq_data/{readtype}_{peaktype}/test_{i}/exp_a_2.fastq.gz"
                            read_data_path_3 = f"/seq_data/{readtype}_{peaktype}/test_{i}/exp_b_1.fastq.gz"
                            read_data_path_4 = f"/seq_data/{readtype}_{peaktype}/test_{i}/exp_b_2.fastq.gz"
                            proj_file_info = textwrap.dedent(f"""
                            Author: Viktoria Haghani & Aditi Goyal & Alan Zhang
                            Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                            Genome:
                              Name: genome
                              Location: /seq_data/{readtype}_{peaktype}/test_{i}/genome.fa
                            Reads:
                              Samples:
                                grp1: 
                                  - {read_data_path_1}
                                  - {read_data_path_2}
                                grp2:
                                  - {read_data_path_3}
                                  - {read_data_path_4}
                              Controls:
                                ctl1: 
                                  - /seq_data/{readtype}_{peaktype}/test_{i}/input_1.fastq.gz
                                  - /seq_data/{readtype}_{peaktype}/test_{i}/input_2.fastq.gz
                            Readtype: {readtype}
                            Peaktype: {peaktype}
                            Aligner: {aligner}
                            Deduplicator: {deduplicator}
                            Peakcaller: {peakcaller}
                            Threads: 1
                            """)
                        elif readtype == "single":
                            read_data_path_1 = f"/seq_data/{readtype}_{peaktype}/test_{i}/exp_a.fastq.gz"
                            read_data_path_2 = f"/seq_data/{readtype}_{peaktype}/test_{i}/exp_b.fastq.gz"
                            proj_file_info = textwrap.dedent(f"""
                            Author: Viktoria Haghani & Aditi Goyal & Alan Zhang
                            Project: exp_vs_obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                            Genome:
                              Name: genome
                              Location: /seq_data/{readtype}_{peaktype}/test_{i}/genome.fa
                            Reads:
                              Samples:
                                grp1: 
                                  - {read_data_path_1}
                                grp2:
                                  - {read_data_path_2}
                              Controls:
                                ctl1: 
                                  - /seq_data/{readtype}_{peaktype}/test_{i}/input.fastq.gz
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
