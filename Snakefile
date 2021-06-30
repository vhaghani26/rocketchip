configfile: "samples.yaml"

wildcard_constraints:
    sample='[a-zA-Z0-9]+'

print(f'Starting ChIP-seq data analysis workflow for samples: {config["samples"]}')
   
rule all:
    input:     
        expand("02_fastqc_analysis/{sample}_1_fastqc.html", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_1_fastqc.zip", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.html", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.zip", sample=config["samples"]),
        expand("05_bigwig_files/{sample}.bw", sample=config["samples"])

rule make_directories:
    message: "Making directories for data organization"
    output:
        directory("00_logs/"),
        directory("01_raw_data/"),
        directory("02_fastqc_analysis/"),
        directory("03_sam_files/"),
        directory("04_bam_files/"),
        directory("05_bigwig_files/"),
        directory("06_macs2_peaks/")
    shell: """
        mkdir 00_logs
        mkdir 01_raw_data
        mkdir 02_fastqc_analysis
        mkdir 03_sam_files
        mkdir 04_bam_files
        mkdir 05_bigwig_files
        mkdir 06_macs2_peaks
    """

rule download_data_wc:
    input: expand("01_raw_data/{sample}/{sample}.sra", sample=config["samples"])

rule download_data:
    message: "Downloading raw data files"
    conda: "00_conda_software/chip_sra.yml"
    output: "01_raw_data/{sample}/{sample}.sra"
    log: "00_logs/{sample}_download_data.log"
    shell: """
    prefetch {wildcards.sample} > {log}
    mv {wildcards.sample}/ 01_raw_data/
    """

rule split_paired_reads_wc:
    input:
        expand("01_raw_data/{sample}_1.fastq.gz", sample=config["samples"]),
        expand("01_raw_data/{sample}_2.fastq.gz", sample=config["samples"])
   
rule split_paired_reads:
    message: "Splitting paired end reads into separate files"
    conda: "00_conda_software/chip_sra.yml"
    input: "01_raw_data/{sample}/{sample}.sra"
    output:
        "01_raw_data/{sample}_1.fastq.gz",
        "01_raw_data/{sample}_2.fastq.gz"
    log: "00_logs/{sample}_split_paired_reads.log"
    shell: "fastq-dump {input} --split-files --gzip --outdir 01_raw_data/ 2> {log}"
    
rule fastqc_precheck:
    message: "Running quality control on samples pre-processing"
    conda: "00_conda_software/chip_fastqc.yml"
    input:
        r1 = "01_raw_data/{sample}_1.fastq.gz",
        r2 = "01_raw_data/{sample}_2.fastq.gz"
    output:
        "02_fastqc_analysis/{sample}_1_fastqc.html",
        "02_fastqc_analysis/{sample}_1_fastqc.zip",
        "02_fastqc_analysis/{sample}_2_fastqc.html",
        "02_fastqc_analysis/{sample}_2_fastqc.zip"
    log: 
        r1 = "00_logs/{sample}_fastqc_precheck_r1.log",
        r2 = "00_logs/{sample}_fastqc_precheck_r2.log"
    shell: """
    fastqc {input.r1} --outdir 02_fastqc_analysis/ 2> {log.r1}
    fastqc {input.r2} --outdir 02_fastqc_analysis/ 2> {log.r2}
    """  

rule download_genome:
    message: "Downloading GRCm39/mm39 mouse genome from the UCSC Genome Browser"
    output: "01_raw_data/mm39.chromFa.tar.gz"
    log: "00_logs/download_genome.log"
    shell: "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz -O {output} 2> {log}"

rule process_genome:
    message: "Decompressing genome. Concatenating individual chromosome files to create full assembly. Removing chromosome sequence files"
    input: "01_raw_data/mm39.chromFa.tar.gz"
    output: "01_raw_data/mm39.fa"
    log: "00_logs/process_genome.log"
    shell: """
    tar zvfx {input} --directory 01_raw_data/
    cat 01_raw_data/*.fa > {output} 2> {log}
    rm 01_raw_data/chr*.fa
    """
    
rule set_alignment_reference:
    message: "Setting GRCm39/mm39 mouse genome assembly as reference genome for alignment" 
    conda: "00_conda_software/chip_bwa.yml"
    input: "01_raw_data/mm39.fa"
    output: multiext("01_raw_data/mm39", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log: "00_logs/set_alignment_reference.log"
    shell: """
    bwa index -p mm39 -a bwtsw {input} 2> {log}
    mv mm39* 01_raw_data/
    """ 
    
rule align_reads_wc:
    input: 
        expand("03_sam_files/{sample}.sam", sample=config["samples"])
        
rule align_reads:
    message: "Aligning paired end reads to GRCm39/mm39 reference genome"
    conda: "00_conda_software/chip_bwa.yml"
    input:
        r1 = "01_raw_data/{sample}_1.fastq.gz",
        r2 = "01_raw_data/{sample}_2.fastq.gz",
        genome = multiext("01_raw_data/mm39", ".amb", ".ann", ".bwt", ".pac", ".sa")
    output: "03_sam_files/{sample}.sam"
    log: "00_logs/{sample}_align_reads_err.log"
    shell: "bwa mem 01_raw_data/mm39 {input.r1} {input.r2} > {output} 2> {log}"
    
rule sam_to_bam:
    message: "Converting SAM to BAM file format"
    conda: "00_conda_software/chip_samtools.yml"
    input: "03_sam_files/{sample}.sam"
    output: "04_bam_files/{sample}.bam"
    log: "00_logs/{sample}_sam_to_bam.log"
    shell: "samtools view -b {input} > {output} 2> {log}"

rule sam_fixmate:
    message: "Removing secondary and unmapped reads. Adding tags to reads for deduplication"
    conda: "00_conda_software/chip_samtools.yml"
    input: "04_bam_files/{sample}.bam"
    output: "04_bam_files/{sample}.namesorted.fixmate.bam"
    log: "00_logs/{sample}_sam_fixmate.log"
    shell: "samtools fixmate -rcm -O bam {input} {output} 2> {log}"

rule sam_sort:
    message: "Sorting reads by chromosome coordinates"
    conda: "00_conda_software/chip_samtools.yml"
    input: "04_bam_files/{sample}.namesorted.fixmate.bam"
    output: "04_bam_files/{sample}.coorsorted.fixmate.bam"
    log: "00_logs/{sample}_sam_sort.log"
    shell: "samtools sort {input} -o {output} 2> {log}"

rule sam_markdup:
    message: "Marking and removing duplicates"
    conda: "00_conda_software/chip_samtools.yml"
    input: "04_bam_files/{sample}.coorsorted.fixmate.bam"
    output: "04_bam_files/{sample}.coorsorted.dedup.bam"
    log: "00_logs/{sample}_sam_markdup.log"
    shell: "samtools markdup -r --mode s {input} {output} 2> {log}"

rule sam_index:
    message: "Indexing deduplicated BAM file"
    conda: "00_conda_software/chip_samtools.yml"
    input: "04_bam_files/{sample}.coorsorted.dedup.bam"
    output: "04_bam_files/{sample}.coorsorted.dedup.bam.bai", 
    log: "00_logs/{sample}_sam_index.log"
    shell: "samtools index {input} 2> {log}"

rule bam_to_bigwig_wc:
    input:
        expand("05_bigwig_files/{sample}.bw", sample=config["samples"])
        
rule bam_to_bigwig:
    message: "Converting BAM file format to bigwig file format for visualization"
    conda: "00_conda_software/chip_deeptools.yml"
    input:
        index = "04_bam_files/{sample}.coorsorted.dedup.bam.bai",
        bam = "04_bam_files/{sample}.coorsorted.dedup.bam"
    output: "05_bigwig_files/{sample}.bw"
    log: "00_logs/{sample}_bam_to_bigwig.log"
    shell: "bamCoverage -b {input} -o {output} 2> {log}"

rule call_peaks_wc:
    input:
        "06_macs2_peaks/{sample}_peaks.narrowPeak",
        "06_macs2_peaks/{sample}_peaks.xls",
        "06_macs2_peaks/{sample}_summits.bed",
        "06_macs2_peaks/{sample}_model.R",
        "06_macs2_peaks/{sample}_control_lambda.bdg",
        "06_macs2_peaks/{sample}_treat_pileup.bdg"')

# Need to run through with one sample to make sure outputs work
rule call_peaks:
    message: "Calling ChIP-seq peaks"
    conda: "00_conda_software/chip_macs2.yml"
    input: "04_bam_files/{sample}.coorsorted.dedup.bam"
    output: multiext('"06_macs2_peaks/{sample}", "_peaks.narrowPeak", "_peaks.xls", "_summits.bed", "_model.R", "_control_lambda.bdg", "_treat_pileup.bdg"')
    log: "00_logs/{sample}_macs2_peaks.log"
    shell: "macs2 callpeak -t {input} -f BAM -n {wildcards.sample} --outdir 06_macs2_peaks/ 2> {log}"