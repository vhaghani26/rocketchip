#!/usr/bin/env python3

readtypes = ["paired", "single"]
peaktypes = ["narrow", "broad"]
aligners = ["bwa_mem", "bowtie2", "STAR"]
peakcallers = ["macs3", "cisgenome", "genrich", "pepr"]
deduplicators = ["samtools", "umi_tools", "sambamba", "picard"]

for readtype in readtypes:
    for peaktype in peaktypes:
        for aligner in aligners:
            for peakcaller in peakcallers:
                for deduplicator in deduplicators:
                    print(readtype, peaktype, aligner, peakcaller, deduplicator)




"""
print(f'Author: Viktoria Haghani & Aditi Goyal/n 
        Project: Expected vs Observed Peaks/n
        Genome:/n
          Name: /n
          Location: /n
        Reads:/n
          Samples:/n
            grp1: /n
            grp2:/n
          Controls:/n
            ctl1: /n
        Readtype: /n
        Peaktype: /n
        Aligner: /n
        Deduplicator: /n
        Peakcaller: /n
        Threads:')
"""