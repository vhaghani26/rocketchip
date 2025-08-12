#!/usr/bin/env python3

import ssl
import urllib.request
import argparse
import yaml
import sys
import os
import re
import glob

try:
    from importlib.metadata import version # For Python > 3.8
except ImportError:
    from importlib_metadata import version  # For Python < 3.8

# Get the version
def get_version():
    try:
        return version("rocketchip")
    except Exception:
        return "unknown"
        
# Define the argument parser
parser = argparse.ArgumentParser(description='Make Snakefiles')
parser.add_argument('configfile', type=str, metavar='<path>',
    help='Path to configuration file. See README for details')
parser.add_argument('--data', type=str, metavar='<str>', required=False, default=os.getcwd(),
    help='Location to store data files (i.e. genome and raw read files); default: current working directory')
parser.add_argument('--output_file', type=str, metavar='<str>', required=False,
    help='Output snakefile name (default: STDOUT)')
parser.add_argument('--version', action='version', version=f'%(prog)s {get_version()}',
                   help='Print the version and exit')

# Fix SSL Certificate bug (https://stackoverflow.com/questions/27835619/urllib-and-ssl-certificate-verify-failed-error, https://stackoverflow.com/questions/35569042/ssl-certificate-verify-failed-with-python3)
context = ssl._create_unverified_context()


###############
## Functions ##
###############
   
# Check for any invalid titles in config file
def read_config(config):
    # Titles
    titles = ['Author', 'Project', 'Genome', 'Reads', 'Readtype', 'Peaktype',
    'Aligner', 'Deduplicator', 'Peakcaller', 'Threads', 'Molecule']
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

    # Molecule
    molecule = config.get('Molecule', 'DNA').upper()
    if molecule not in ['DNA', 'RNA']:
        molecule = 'DNA'
    
    return (
        config['Genome']['Name'].replace(' ', '_'),
        config['Genome']['Location'],
        config['Readtype'],
        config['Peaktype'],
        config['Aligner'],
        config['Deduplicator'],
        config['Peakcaller'],
        config['Reads']['Controls'] is not None,
        config['Threads'],
        molecule
    )

def get_read_pair(sample_id, read_dir, readtype):
    patterns = [
        # Common paired end patterns
        (f'{sample_id}_1.fastq.gz', f'{sample_id}_2.fastq.gz'),
        (f'{sample_id}_R1.fastq.gz', f'{sample_id}_R2.fastq.gz'),
        (f'{sample_id}_1.fq.gz', f'{sample_id}_2.fq.gz'),
        (f'{sample_id}_R1.fq.gz', f'{sample_id}_R2.fq.gz'),
        (f'{sample_id}_R1_001.fastq.gz', f'{sample_id}_R2_001.fastq.gz'),
        (f'{sample_id}_R1_001.fq.gz', f'{sample_id}_R2_001.fq.gz'),
        # Flexible paired end patterns
        (f'{sample_id}*_1.fastq.gz', f'{sample_id}*_2.fastq.gz'),
        (f'{sample_id}*_R1.fastq.gz', f'{sample_id}*_R2.fastq.gz'),
        (f'{sample_id}*_1.fq.gz', f'{sample_id}*_2.fq.gz'),
        (f'{sample_id}*_R1.fq.gz', f'{sample_id}*_R2.fq.gz'),
        (f'{sample_id}*_R1_001.fastq.gz', f'{sample_id}*_R2_001.fastq.gz'),
        (f'{sample_id}*_R1_001.fq.gz', f'{sample_id}*_R2_001.fq.gz'),
        # Common single end patterns
        (f'{sample_id}.fastq.gz', None),
        (f'{sample_id}.fq.gz', None),
        # Flexible single end patterns
        (f'{sample_id}*.fastq.gz', None),
        (f'{sample_id}*.fq.gz', None),
    ]
    
    for r1_pattern, r2_pattern in patterns:
        r1_matches = glob.glob(os.path.join(read_dir, r1_pattern))
        if r1_matches:
            if readtype == 'paired' and r2_pattern:
                r2_matches = glob.glob(os.path.join(read_dir, r2_pattern))
                if r2_matches:
                    return (r1_matches[0], r2_matches[0])
            else:
                return (r1_matches[0], None)
    return (None, None)

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

####################
## Main Execution ##
####################

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
    
    lcr = {}
    sra = []
    READ_SRC = None
    sample_reads = {}

    for sample_type in config['Reads']:
        if config['Reads'][sample_type] is None:
            continue
            
        for group in config['Reads'][sample_type]:
            sample_reads[group] = {}  # Initialize group
            for read in config['Reads'][sample_type][group]:
                read_id = read  # Use consistent variable name
                read_name = os.path.basename(read_id)  # For SRA handling
                
                # Handle SRA entries
                if read_name.startswith(('SRR', 'ERR', 'DRR')):
                    if read_name in sra: continue
                    sys.stderr.write(f'SRA entry detected... assessing read type (paired/single) for {read_name}\n')
                    page = urllib.request.urlopen(f'https://www.ncbi.nlm.nih.gov/sra/?term={read_name}', context=context)
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
                    continue
                    
                # Handle local files
                read_dir = os.path.dirname(read_id) if os.path.isabs(read_id) else os.path.join(CONFIGFILE_DIR, os.path.dirname(read_id))
                sample_id = os.path.basename(read_id)
                
                # First try exact matching
                r1, r2 = get_read_pair(sample_id, read_dir, READTYPE)
                
                # If no matches, try more flexible matching
                if r1 is None:
                    all_files = [f for f in os.listdir(read_dir) 
                               if f.endswith(('.fastq.gz', '.fq.gz'))]
                    matching_files = [f for f in all_files 
                                    if f.startswith(sample_id)]
                    
                    if READTYPE == 'paired':
                        r1_candidates = [f for f in matching_files 
                                       if '_1.' in f or '_R1.' in f]
                        r2_candidates = [f for f in matching_files 
                                       if '_2.' in f or '_R2.' in f]
                        if r1_candidates and r2_candidates:
                            r1 = os.path.join(read_dir, r1_candidates[0])
                            r2 = os.path.join(read_dir, r2_candidates[0])
                    else:
                        r1_candidates = [f for f in matching_files 
                                        if not ('_2.' in f or '_R2.' in f)]
                        if r1_candidates:
                            r1 = os.path.join(read_dir, r1_candidates[0])
                
                # Verify files are found
                if r1 is None:
                    sys.exit(f'Error: Could not find read files for sample {sample_id} in {read_dir}')
                    
                # Store for verification and processing
                sample_reads[group][sample_id] = {
                    'forward': r1,
                    'reverse': r2 if READTYPE == 'paired' else None
                }
                
                # Store in lcr dict
                if READTYPE == 'paired':
                    lcr[sample_id] = [r1, r2]
                else:
                    lcr[sample_id] = r1

    ##########################
    ## Verify File Matching ##
    ##########################

    print("\nSample File Assignment Verification:\n")
    print(f"Read Type: {READTYPE}\n")

    for group, samples in sample_reads.items():
        print(f"Group: {group}")
        for sample_id, reads in samples.items():
            print(f"  Sample ID: {sample_id}")
            print(f"    Forward Reads: {reads['forward']}")
            if reads['reverse']:
                print(f"    Reverse Reads: {reads['reverse']}")
        print()

    confirm = input("Are these file assignments correct? (y/n) ").lower()
    if confirm not in ('y', 'yes'):
        print("Please edit the names of your input files (dirty fix) or config file (recommended fix), then try again.")
        sys.exit(1)

    if len(lcr) > 0 and len(sra) > 0:
        sys.exit(f'Error: Mix of local fastq reads and remote SRR identifiers detected')
    if len(lcr) <= 0 and len(sra) <= 0:
        sys.exit(f'Error: 0 reads detected')
    if len(lcr) > 0: 
        READ_SRC = 'lcr'
    else: 
        READ_SRC = 'sra'

    #########################################
    ## Check for peakcaller specific erros ##
    #########################################

    # Cisgenome
    if PEAKCALLER == 'cisgenome':
        # Accept no replicates when there is 0 control
        if config['Reads']['Controls'] is None:
            sys.exit(f'Error: cisgenome requires at least 1 control sample; no control sample detected')
                    
    # PePr
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
    if ALIGNER == 'STAR':
        molecule_type = 'RNA' if molecule == 'RNA' else 'DNA'
        ALIGN_READS = f'{RULES_DIR}/align_reads/STAR_{molecule_type}_{READTYPE}.txt'
    else:
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