#!/usr/bin/env python3

import ssl
import urllib.request
import argparse
import yaml
import sys
import os
import re

# Define the argument parser
parser = argparse.ArgumentParser(description='Make Snakefiles')
parser.add_argument('configfile', type=str, metavar='<path>',
    help='Path to configuration file. See README for details')
parser.add_argument('--data', type=str, metavar='<str>', required=False, default=os.getcwd(),
    help='Location to store data files (i.e. genome and raw read files); default: current working directory')
parser.add_argument('--output_file', type=str, metavar='<str>', required=False,
    help='Output snakefile name (default: STDOUT)')

# Fix SSL Certificate bug (https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error, https://stackoverflow.com/questions/35569042/ssl-certificate-verify-failed-with-python3)
context = ssl._create_unverified_context()


###############
## Functions ##
###############

# Check for any invalid titles in config file
def read_config(config):
    # Titles
    titles = ['Author', 'Project', 'Genome', 'Reads', 'Readtype', 'Peaktype',
    'Aligner', 'Deduplicator', 'Peakcaller', 'Threads']
    for title in config:
        if title not in titles:
            sys.exit(f'Error: Invalid title in configfile: {title}')
        if config[title] is None:
            sys.exit(f'Error: No input detected in config file in section \'{title}\'') 
    
    # Genome
    genome_info = ['Name', 'Location']
    for title in config['Genome']:
        if title not in genome_info:
            sys.exit(f'Error: Invalid title in configfile: {title}')
        if config['Genome'][title] is None:
            sys.exit(f'Error: No input detected in config file in section \'Genome: {title}\'')
    
    # Reads
    reads_info  = ['Samples', 'Controls']
    for title in config['Reads']:
        if title not in reads_info:
            sys.exit(f'Error: Invalid title in config file: {title}')
    
    # User input
    if config['Reads']['Samples'] is None:
        sys.exit(f'Error: Invalid sample input: no sample group detected')
    else:
        for group in config['Reads']['Samples']:
            spl_group = config['Reads']['Samples'][group]
            if spl_group is None or len(spl_group) < 1:
                sys.exit(f'Error: Invalid sample input: no samples detected in group {group}')
    
    if config['Reads']['Controls'] is not None:
        for group in config['Reads']['Controls']:
            ctl_group = config['Reads']['Controls'][group]
            if ctl_group is None or len(ctl_group) < 1:
                sys.exit(f'Error: Invalid control sample input: no samples detected in group {group}')
            
    options = {
        'Readtype': ['single', 'paired'],
        'Peaktype': ['narrow', 'broad'],
        'Aligner': ['bwa_mem', 'bowtie2', 'STAR'],
        'Deduplicator': ['samtools', 'picard', 'sambamba', 'no_deduplication'],
        'Peakcaller': ['macs3', 'genrich', 'pepr', 'cisgenome']
        }
    for option in options:
        user_input = config[option]
        if user_input not in options[option]:
            sys.exit(f'Error: Invalid {option} in configfile: {user_input}\nInput one of the following {option}s: {options[option]} (case sensitive)')

    
    return (
        config['Genome']['Name'].replace(' ', '_'),
        config['Genome']['Location'],
        config['Readtype'],
        config['Peaktype'],
        config['Aligner'],
        config['Deduplicator'],
        config['Peakcaller'],
        config['Reads']['Controls'] is not None,
        config['Threads']
    )


def dict2snake(reads, dict_type, read_src):
    o = []
    if dict_type == 'names':
        s = 'Samples'
        o.append(f'{s} = {{\n')
        if read_src == 'sra':
            for group in reads[s]: o.append(f'\t\'{group}\': {reads[s][group]},\n')
        else:
            for group in reads[s]:
                read_names = []
                for read in reads[s][group]:
                    read_name = os.path.basename(read)
                    read_name = re.sub(' ', '\ ', read_name)
                    read_names.append(read_name)
                o.append(f'\t\'{group}\': {read_names},\n')
        o.append(f'}}\n')
        
        c = 'Controls'
        if reads['Controls'] is not None:
            o.append(f'{c} = {{\n')
            if read_src == 'sra':
                for group in reads[c]: o.append(f'\t\'{group}\': {reads[c][group]},\n')
            else:
                for group in reads[c]:
                    read_names = []
                    for read in reads[c][group]:
                        read_name = os.path.basename(read)
                        read_name = re.sub(' ', '\ ', read_name)
                        read_names.append(read_name)
                    o.append(f'\t\'{group}\': {read_names},\n')
            o.append(f'}}\n')
            
    elif dict_type == 'paths':
        o.append(f'reads_path = {{\n')
        for sample in reads:
            if type(reads[sample]) == list:
                for item, i in zip(reads[sample], range(len(reads[sample]))):
                    if ' ' in item:
                        reads[sample][i] = item
                        samps = "', '".join(reads[sample])
                if ' ' in reads[sample][0]:
                    o.append(f'\t\'{sample}\': [\'{samps}\'],\n')
                else:
                    o.append(f'\t\'{sample}\': {reads[sample]},\n')
            else:
                o.append(f'\t\'{sample}\': \'{reads[sample]}\',\n')
                
        o.append(f'}}\n')
    
    return ''.join(o)

# Main function
def main():
    arg = parser.parse_args()

    #####################################
    ## Environment and file  locations ##
    #####################################

    # Data environment
    if arg.data: os.environ['ROCKETCHIP_DATA'] = arg.data
    if 'ROCKETCHIP_DATA' not in os.environ:
        sys.exit('Error: ROCKETCHIP_DATA not in environment; set ROCKETCHIP_DATA or use --data')
    if not os.path.isdir(os.environ['ROCKETCHIP_DATA']):
        sys.exit('Error: ROCKETCHIP_DATA is not a directory')
    DATA = os.path.abspath(os.environ['ROCKETCHIP_DATA'])
    DATA = re.sub(' ', '\ ', DATA)

    # Source environment
    SRC = os.path.dirname(os.path.abspath(__file__)) 
    SRC = re.sub(' ', '\ ', SRC)

    # Current and config file locations
    CURRENT_DIR     = os.getcwd()
    CURRENT_DIR = re.sub(' ', '\ ', CURRENT_DIR)

    CONFIGFILE_PATH = os.path.abspath(arg.configfile)
    CONFIGFILE_PATH = re.sub(' ', '\ ', CONFIGFILE_PATH)

    CONFIGFILE_DIR  = os.path.dirname(CONFIGFILE_PATH)
    CONFIGFILE_DIR = re.sub(' ', '\ ', CONFIGFILE_DIR)

    RULES_DIR       = f'{SRC}/rules'

    #####################
    ## Read Configfile ##
    #####################

    config = None
    with open(arg.configfile, 'r') as yamlfile: config = yaml.safe_load(yamlfile)
    (GENOME_NAME, GENOME_LOCATION, READTYPE, PEAKTYPE,
    ALIGNER, DEDUPLICATOR, PEAKCALLER, CONTROL, THREADS) = read_config(config)

    #######################
    ## Genome Management ##
    #######################
    if os.path.isfile(f'{CONFIGFILE_DIR}/{GENOME_LOCATION}'):
        GENOME_LOCATION = os.path.abspath(f'{CONFIGFILE_DIR}/{GENOME_LOCATION}')
        GENOME_LOCATION = re.sub(' ', '\ ', GENOME_LOCATION)
    if not os.path.isabs(GENOME_LOCATION) and not GENOME_LOCATION.startswith(('http://', 'https://')):
        sys.exit('Error: Invalid path to local genome or url to remote genome in configfile in section Genome: Location')

    ######################
    ## Reads Management ##
    ######################
    lcr  = {}
    sra = []
    READ_SRC = None
    for sample_type in config['Reads']:
        if config['Reads'][sample_type] is None: continue
        for group in config['Reads'][sample_type]:
            for read in config['Reads'][sample_type][group]:
                read_dir  = os.path.dirname(read)
                read_name = os.path.basename(read)
                read_name = re.sub(' ', '\ ', read_name)
                # Handle SRA entries
                if read_name.startswith(('SRR', 'ERR', 'DRR')):
                    if read_name in sra: continue
                    sys.stderr.write(f'SRA entry detected... assessing read type (paired/single) for {read_name}\n')
                    page = urllib.request.urlopen(f'https://www.ncbi.nlm.nih.gov/sra/?term={read_name}', context = context)
                    html = page.read().decode('UTF-8')
                    page.close()
                    match = re.search('Layout: <span>(.{6})</span>', html)
                    readtype = None
                    if match:
                        if re.search('SINGLE', match.group()): readtype = 'single'
                        elif re.search('PAIRED', match.group()): readtype = 'paired'
                    if readtype is None:
                        sys.exit(f'Could not assess paired/end for {read_name}, aborting')
                    if readtype != READTYPE:
                        sys.exit(f'Inconsistent read type...\n{read_name} has type: {readtype} while config file input is Readtype: {READTYPE}')
                    sys.stderr.write(f'{read_name} has read type: {readtype}\n')
                    sra.append(read_name)
                # Handle local data
                elif os.path.isdir(f'{CONFIGFILE_DIR}/{read_dir}') or os.path.isabs(read_dir):
                    if READTYPE == 'single':
                        if os.path.isabs(read_dir):
                            read_path = f'{read_dir}/{read_name}.fastq.gz'
                        else:
                            read_path = f'{CONFIGFILE_DIR}/{read_dir}/{read_name}.fastq.gz'
                        tmp_read_path = read_path.split(r'\ ')
                        tmp_read_path = ' '.join(tmp_read_path)                    
                        if os.path.isfile(tmp_read_path):
                            if read_name not in lcr: lcr[read_name] = read_path
                        else:
                            sys.exit(f'Error: Missing {read_name}.fastq.gz in {read_dir}; Double check input for Readtype or Reads')
                    elif READTYPE == 'paired':
                        if os.path.isabs(read_dir):
                            read_1_path = f'{read_dir}/{read_name}_1.fastq.gz'
                            read_2_path = f'{read_dir}/{read_name}_2.fastq.gz'
                        else:
                            read_1_path = f'{CONFIGFILE_DIR}/{read_dir}/{read_name}_1.fastq.gz'
                            read_2_path = f'{CONFIGFILE_DIR}/{read_dir}/{read_name}_2.fastq.gz'
                        tmp_read_1_path = read_1_path.split(r'\ ')
                        tmp_read_1_path = ' '.join(tmp_read_1_path)
                        tmp_read_2_path = read_2_path.split(r'\ ')
                        tmp_read_2_path = ' '.join(tmp_read_2_path)
                        two_reads = os.path.exists(tmp_read_1_path) and os.path.exists(tmp_read_2_path)
                        if two_reads:
                            if read_name not in lcr:
                                lcr[read_name] = [read_1_path, read_2_path]
                        else:
                            sys.exit(f'Error: Missing one or both reads from {read_name} in {read_dir}; Double check input for Readtype or Reads')
                    else:
                        sys.exit(f'Error: Invalid read type: {READTYPE}')
                else:
                    sys.exit(f'Error: Invalid read: {read_name}; neither SRA identifier nor detected as local read')

    if len(lcr) > 0 and len(sra) > 0:
        sys.exit(f'Error: Mix of local fastq reads and remote SRR identifiers detected')
    if len(lcr) <= 0 and len(sra) <= 0:
        sys.exit(f'Error: 0 reads detected')
    if len(lcr) > 0: READ_SRC = 'lcr'
    else: READ_SRC = 'sra'

    #########################################
    ## Check for peakcaller specific erros ##
    #########################################

    # cisgenome
    if PEAKCALLER == 'cisgenome':
        # accept no replicates when there is 0 control
        if config['Reads']['Controls'] is None:
            sys.exit(f'Error: cisgenome requires at least 1 control sample; no control sample detected')
                    
    # pepr
    if PEAKCALLER == 'pepr':
        for spl_group in config['Reads']['Samples']:
            num_samples = len(config['Reads']['Samples'][spl_group])
            if num_samples < 2:
                sys.exit(f'Error: pepr must read in at least 2 sample replicates per experiment group and at least 1 control sample; {num_samples} sample(s) provided in experiment group \'{spl_group}\'')
        if config['Reads']['Controls'] is None:
            sys.exit(f'Error: pepr must read in at least 2 expriment replicates per experiment group and at least 1 control sample; 0 control provided')

    ######################
    ## Make Snake Rules ##
    ######################

    snakerules = []

    # Project file content
    project_file = '\'\'\'\n' + yaml.dump(config, sort_keys=False) + '\'\'\'\n'
    snakerules.append(project_file)

    # Wildcard constraints
    wildcard_constraints = f'wildcard_constraints:\n\tsample=\'[a-zA-Z0-9_]+\',\n\tcontrol=\'[a-zA-Z0-9_]+\'\n'
    snakerules.append(wildcard_constraints)

    # Reads
    reads = dict2snake(config['Reads'], dict_type = 'names', read_src = READ_SRC)
    snakerules.append(reads)
    if READ_SRC == 'lcr':
        lcr_path =  dict2snake(lcr, dict_type = 'paths', read_src = READ_SRC)
        snakerules.append(lcr_path)

    snakerules.append(f'src = \'{SRC}\'\n')

    # Read parser
    READ_PARSER = None
    if CONTROL:
        READ_PARSER = f'{RULES_DIR}/parse_reads/parse_reads_with_control.txt'
    else:
        READ_PARSER = f'{RULES_DIR}/parse_reads/parse_reads_no_control.txt'
    tmp_READ_PARSER = READ_PARSER.split(r'\ ')
    tmp_READ_PARSER = ' '.join(tmp_READ_PARSER)
    with open(tmp_READ_PARSER) as fh: snakerules.append(fh.read())

    # Rule all: fastq files, peak files, bigwig file
    RULE_ALL_FASTQ = None
    RULE_ALL_PEAKCALLER = None
    RULE_ALL_BIGWIG = f'{RULES_DIR}/rule_all/bigwig/bigwig.txt'
    if CONTROL:
        RULE_ALL_FASTQ = f'{RULES_DIR}/rule_all/fastq/fastq_with_control_{READTYPE}.txt'
        RULE_ALL_PEAKCALLER = f'{RULES_DIR}/rule_all/peaks/{PEAKCALLER}_with_control.txt'
    else:
        RULE_ALL_FASTQ = f'{RULES_DIR}/rule_all/fastq/fastq_no_control_{READTYPE}.txt'
        RULE_ALL_PEAKCALLER = f'{RULES_DIR}/rule_all/peaks/{PEAKCALLER}_no_control.txt'
    for rule in (RULE_ALL_FASTQ, RULE_ALL_PEAKCALLER, RULE_ALL_BIGWIG):
        tmp_rule = rule.split(r'\ ')
        tmp_rule = ' '.join(tmp_rule)
        with open(tmp_rule) as fh: snakerules.append(fh.read())

    # Make directories
    MAKE_DIR = f'{RULES_DIR}/make_dir/make_dir.txt'
    tmp_MAKE_DIR = MAKE_DIR.split(r'\ ')
    tmp_MAKE_DIR = ' '.join(tmp_MAKE_DIR)
    with open(tmp_MAKE_DIR) as fh:
        rule = fh.read()
        assert('PEAKCALLER' in rule)
        rule = rule.replace('PEAKCALLER', PEAKCALLER)
        snakerules.append(rule)

    # Download/copy genome:
    GET_GENOME = None
    if os.path.isabs(GENOME_LOCATION): GET_GENOME = f'{RULES_DIR}/get_genome/local_genome.txt'
    elif GENOME_LOCATION.endswith(f'{GENOME_NAME}/bigZips/{GENOME_NAME}.fa.gz'): GET_GENOME = f'{RULES_DIR}/get_genome/remote_fa_gz_genome.txt'
    elif GENOME_LOCATION.endswith(f'{GENOME_NAME}/bigZips/chromFa.tar.gz'): GET_GENOME = f'{RULES_DIR}/get_genome/remote_tar_gz_genome.txt'
    elif GENOME_LOCATION.endswith(f'hg38.chromFa.tar.gz'): GET_GENOME = f'{RULES_DIR}/get_genome/remote_hg38_genome.txt'
    else: GET_GENOME = f'{RULES_DIR}/get_genome/remote_tar_gz_genome.txt'
    tmp_GET_GENOME = GET_GENOME.split(r'\ ')
    tmp_GET_GENOME = ' '.join(tmp_GET_GENOME)
    with open(tmp_GET_GENOME) as fh:
        rule = fh.read()
        assert('GENOME_LOCATION' in rule and 'GENOME' in rule and 'DATA' in rule)
        rule = rule.replace('DATA', DATA)
        rule = rule.replace('GENOME_LOCATION', GENOME_LOCATION)
        rule = rule.replace('GENOME', GENOME_NAME)
        snakerules.append(rule)
        
    # Index genome:
    INDEX_GENOME = f'{RULES_DIR}/index_genome/{ALIGNER}.txt'
    tmp_INDEX_GENOME = INDEX_GENOME.split(r'\ ')
    tmp_INDEX_GENOME = ' '.join(tmp_INDEX_GENOME)
    with open(tmp_INDEX_GENOME) as fh:
        rule = fh.read()
        assert('GENOME' in rule and 'DATA' in rule)
        rule = rule.replace('DATA', DATA)
        rule = rule.replace('GENOME', GENOME_NAME)
        snakerules.append(rule)

    # Download/link reads:
    if len(sra) > 0:
        GET_READS = f'{RULES_DIR}/get_reads/remote_reads_{READTYPE}.txt'
        tmp_GET_READS = GET_READS.split(r'\ ')
        tmp_GET_READS = ' '.join(tmp_GET_READS)
        with open(tmp_GET_READS) as fh:
            rule = fh.read()
            assert('DATA' in rule)
            rule = rule.replace('DATA', DATA)
            snakerules.append(rule)
    else:
        GET_READS = f'{RULES_DIR}/get_reads/local_reads_{READTYPE}.txt'
        tmp_GET_READS = GET_READS.split(r'\ ')
        tmp_GET_READS = ' '.join(tmp_GET_READS)
        with open(tmp_GET_READS) as fh: snakerules.append(fh.read())


    # Fastq preprocessing
    FASTQ_PREPROCESS = f'{RULES_DIR}/fastqc/preprocessing_{READTYPE}.txt'
    tmp_FASTQ_PREPROCESS = FASTQ_PREPROCESS.split(r'\ ')
    tmp_FASTQ_PREPROCESS = ' '.join(tmp_FASTQ_PREPROCESS)
    with open(tmp_FASTQ_PREPROCESS) as fh: snakerules.append(fh.read())

    # Align reads
    ALIGN_READS = f'{RULES_DIR}/align_reads/{ALIGNER}_{READTYPE}.txt'
    tmp_ALIGN_READS = ALIGN_READS.split(r'\ ')
    tmp_ALIGN_READS = ' '.join(tmp_ALIGN_READS)
    with open(tmp_ALIGN_READS) as fh:
        rule = fh.read()
        assert('GENOME' in rule and 'THREADS' in rule)
        rule = rule.replace('GENOME', GENOME_NAME)
        rule = rule.replace('THREADS', str(THREADS))
        snakerules.append(rule)
        
    # Sort reads
    SORT_READS = f'{RULES_DIR}/sort_index_reads/samtools_sort.txt'
    tmp_SORT_READS = SORT_READS.split(r'\ ')
    tmp_SORT_READS = ' '.join(tmp_SORT_READS)
    with open(tmp_SORT_READS) as fh:
        rule = fh.read()
        assert('THREADS' in rule)
        rule = rule.replace('THREADS', str(THREADS))
        snakerules.append(rule)

    # Markdup
    MARKDUP = None
    if DEDUPLICATOR == 'no_deduplication': MARKDUP = f'{RULES_DIR}/mark_dup/fake_dedup.txt'
    else: MARKDUP = f'{RULES_DIR}/mark_dup/{DEDUPLICATOR}_dedup.txt'
    tmp_MARKDUP = MARKDUP.split(r'\ ')
    tmp_MARKDUP = ' '.join(tmp_MARKDUP)
    with open(tmp_MARKDUP) as fh:
        rule = fh.read()
        if DEDUPLICATOR in ('sambamba', 'samtools'):
            assert('THREADS' in rule)
            rule = rule.replace('THREADS', str(THREADS))
        snakerules.append(rule)

    # Index reads
    INDEX_READS = f'{RULES_DIR}/sort_index_reads/samtools_index.txt'
    tmp_INDEX_READS = INDEX_READS.split(r'\ ')
    tmp_INDEX_READS = ' '.join(tmp_INDEX_READS)
    with open(tmp_INDEX_READS) as fh:
        rule = fh.read()
        assert('THREADS' in rule)
        rule = rule.replace('THREADS', str(THREADS))
        snakerules.append(rule)

    # Bam to bigwig
    BAM_TO_BIGWIG = f'{RULES_DIR}/bam_to_bigwig/bam_to_bigwig.txt'
    tmp_BAM_TO_BIGWIG = BAM_TO_BIGWIG.split(r'\ ')
    tmp_BAM_TO_BIGWIG = ' '.join(tmp_BAM_TO_BIGWIG)
    with open(tmp_BAM_TO_BIGWIG) as fh: snakerules.append(fh.read())

    # Call peaks
    CALL_PEAKS = None
    if CONTROL:
        CALL_PEAKS = f'{RULES_DIR}/call_peaks/{PEAKCALLER}/{PEAKCALLER}_with_control.txt'
    else:
        CALL_PEAKS = f'{RULES_DIR}/call_peaks/{PEAKCALLER}/{PEAKCALLER}_no_control.txt'
    tmp_CALL_PEAKS = CALL_PEAKS.split(r'\ ')
    tmp_CALL_PEAKS = ' '.join(tmp_CALL_PEAKS)    
    with open(tmp_CALL_PEAKS) as fh:
        rule = fh.read()
        # cisgenome specific
        if PEAKCALLER == 'cisgenome':
            assert('THREADS' in rule)
            rule = rule.replace('THREADS', str(THREADS))
            if not CONTROL:
                assert('GENOME' in rule)
                rule = rule.replace('GENOME', GENOME_NAME)
            rule = rule.replace('./tools', f'{os.environ["ROCKETCHIP_SRC"]}/tools')
        # genrich specific
        if PEAKCALLER == 'genrich':
            assert('THREADS' in rule and 'GENRICH_READ_TYPE' in rule)
            rule = rule.replace('THREADS', str(THREADS))
            if READTYPE == 'single': rule = rule.replace('GENRICH_READ_TYPE', '-y')
            else: rule = rule.replace('GENRICH_READ_TYPE', '')
        # macs3 specific
        if PEAKCALLER == 'macs3':
            assert('MACS3_READ_TYPE' in rule and 'MACS3_PEAK_TYPE' in rule)
            if READTYPE == 'single': rule = rule.replace('MACS3_READ_TYPE', 'BAM')
            else: rule = rule.replace('MACS3_READ_TYPE', 'BAMPE')
            if PEAKTYPE == 'narrow': rule = rule.replace('MACS3_PEAK_TYPE', '')
            else: rule = rule.replace('MACS3_PEAK_TYPE','--broad')
        # pepr specific
        if PEAKCALLER == 'pepr':
            assert('PEPR_READ_TYPE' in rule and 'PEPR_PEAK_TYPE' in rule and 'THREADS' in rule)
            rule = rule.replace('THREADS', str(THREADS))
            if READTYPE == 'single': rule = rule.replace('PEPR_READ_TYPE','bam')
            else: rule = rule.replace('PEPR_READ_TYPE','bampe')
            if PEAKTYPE == 'narrow': rule = rule.replace('PEPR_PEAK_TYPE','sharp')
            else: rule = rule.replace('PEPR_PEAK_TYPE','broad')
        
        snakerules.append(rule)

    # Fastq postprocessing
    FASTQ_PREPROCESS = f'{RULES_DIR}/fastqc/postprocessing_{READTYPE}.txt'
    tmp_FASTQ_PREPROCESS = FASTQ_PREPROCESS.split(r'\ ')
    tmp_FASTQ_PREPROCESS = ' '.join(tmp_FASTQ_PREPROCESS)  
    with open(tmp_FASTQ_PREPROCESS) as fh: snakerules.append(fh.read())

    if arg.output_file:
        with open(arg.output_file, 'w') as fh: fh.write('\n'.join(snakerules))
    else:
        print('\n'.join(snakerules))

# Entry point
if __name__ == "__main__":
    main()