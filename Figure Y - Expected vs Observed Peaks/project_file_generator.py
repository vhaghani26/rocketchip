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
                    if readtype == "paired" and peaktype == "narrow":
                        proj_file_info = f"""
                        Author: Viktoria Haghani & Aditi Goyal
                        Project: Exp_vs_Obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                        Genome:
                          Name: 
                          Location: 
                        Reads:
                          Samples:
                            grp1: 
                              - 
                              - 
                            grp2:
                              - 
                          Controls:
                            ctl1: 
                              - 
                        Readtype: {readtype}
                        Peaktype: {peaktype}
                        Aligner: {aligner}
                        Deduplicator: {deduplicator}
                        Peakcaller: {peakcaller}
                        Threads: 1
                        """
                    elif readtype == "paired" and peaktype == "broad":
                        proj_file_info = f"""
                        Author: Viktoria Haghani & Aditi Goyal
                        Project: Exp_vs_Obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                        Genome:
                          Name: 
                          Location: 
                        Reads:
                          Samples:
                            grp1: 
                              - 
                              - 
                            grp2:
                              - 
                          Controls:
                            ctl1: 
                              - 
                        Readtype: {readtype}
                        Peaktype: {peaktype}
                        Aligner: {aligner}
                        Deduplicator: {deduplicator}
                        Peakcaller: {peakcaller}
                        Threads: 1
                        """
                    elif readtype == "single" and peaktype == "narrow":
                        proj_file_info = f"""
                        Author: Viktoria Haghani & Aditi Goyal
                        Project: Exp_vs_Obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                        Genome:
                          Name: 
                          Location: 
                        Reads:
                          Samples:
                            grp1: 
                              - 
                              - 
                            grp2:
                              - 
                          Controls:
                            ctl1: 
                              - 
                        Readtype: {readtype}
                        Peaktype: {peaktype}
                        Aligner: {aligner}
                        Deduplicator: {deduplicator}
                        Peakcaller: {peakcaller}
                        Threads: 1
                        """
                    elif readtype == "single" and peaktype == "broad":
                        proj_file_info = f"""
                        Author: Viktoria Haghani & Aditi Goyal
                        Project: Exp_vs_Obs_{readtype}_{peaktype}_{aligner}_{peakcaller}_{deduplicator}
                        Genome:
                          Name: 
                          Location: 
                        Reads:
                          Samples:
                            grp1: 
                              - 
                              - 
                            grp2:
                              - 
                          Controls:
                            ctl1: 
                              - 
                        Readtype: {readtype}
                        Peaktype: {peaktype}
                        Aligner: {aligner}
                        Deduplicator: {deduplicator}
                        Peakcaller: {peakcaller}
                        Threads: 1
                        """
                    print()
                    print(proj_file_info)