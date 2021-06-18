configfile: "samples.yaml"

print(f'Starting ChIP-seq data analysis workflow for samples: {config["samples"]}')
   
rule all:
    input: 
        "01_raw_data/mm39.amb",
        "01_raw_data/mm39.ann",
        "01_raw_data/mm39.bwt",
        "01_raw_data/mm39.pac",
        "01_raw_data/mm39.sa", 
        expand("01_raw_data/{sample}/{sample}.sra", sample=config["samples"]),
        expand("01_raw_data/{sample}_1.fastq.gz", sample=config["samples"]),
        expand("01_raw_data/{sample}_2.fastq.gz", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_1_fastqc.html", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_1_fastqc.zip", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.html", sample=config["samples"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.zip", sample=config["samples"]),
        expand("04_bam_files/{sample}.coorsorted.dedup.bam.bai", sample=config["samples"]),
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

rule download_data:
    message: "Downloading raw data files"
    conda: "chip_seq_environment.yml"
    params:
        lambda wildcards: config[wildcards.samples]
    output: "01_raw_data/{sample}/{sample}.sra"
    log: "00_logs/{sample}_download_data.log"
    shell: """
    echo 'prefetch {params.id} 2> {log}'
    echo 'mv {params.id}/ 01_raw_data/'
    touch {output}
    """
    
rule split_paired_reads:
    message: "Splitting paired end reads into separate files"
    conda: "chip_seq_environment.yml"
    input: "01_raw_data/{sample}/{sample}.sra"
    output:
        "01_raw_data/{sample}_1.fastq.gz",
        "01_raw_data/{sample}_2.fastq.gz"
    log: "00_logs/{sample}_split_paired_reads.log"
    shell: "echo 'fastq-dump {input} --split-files --gzip --outdir 01_raw_data/ 2> 00_logs/{log}'"
    
rule fastqc_precheck_r1:
    message: "Running quality control on samples pre-processing"
    conda: "chip_seq_environment.yml"
    input: "01_raw_data/{sample}_1.fastq.gz"
    output:
        "02_fastqc_analysis/{sample}_1_fastqc.html",
        "02_fastqc_analysis/{sample}_1_fastqc.zip"
    log: "00_logs/{sample}_fastqc_precheck_r1.log"
    shell: "fastqc {input} --outdir 02_fastqc_analysis/ 2> 00_logs/{log}"

rule fastqc_precheck_r2:
    message: "Running quality control on samples pre-processing"
    conda: "chip_seq_environment.yml"
    input: "01_raw_data/{sample}_2.fastq.gz"
    output:
        "02_fastqc_analysis/{sample}_2_fastqc.html",
        "02_fastqc_analysis/{sample}_2_fastqc.zip"
    log: "00_logs/{sample}_fastqc_precheck_r2.log"
    shell: "fastqc {input} --outdir 02_fastqc_analysis/ 2> 00_logs/{log}"   

rule download_genome:
    message: "Downloading GRCm39/mm39 mouse genome from the UCSC Genome Browser"
    output: "01_raw_data/mm39.chromFa.tar.gz"
    log: "00_logs/download_genome.log"
    shell: "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz -O {output} 2> 00_logs/{log}"

rule process_genome:
    message: "Decompressing genome. Concatenating individual chromosome files to create full assembly. Removing chromosome sequence files"
    input: "01_raw_data/mm39.chromFa.tar.gz"
    output: "01_raw_data/mm39.fa"
    log: "00_logs/process_genome.log"
    shell: """
    tar zvfx {input} --directory 01_raw_data/
    cat 01_raw_data/*.fa > {output} 2> 00_logs/{log}
    rm 01_raw_data/chr*.fa
    """
    
rule set_alignment_reference:
    message: "Setting GRCm39/mm39 mouse genome assembly as reference genome for alignment" 
    conda: "chip_seq_environment.yml"
    input: "01_raw_data/mm39.fa"
    output: multiext("01_raw_data/mm39", ".amb", ".ann", ".bwt", ".pac", ".sa")
    log: "00_logs/set_alignment_reference.log"
    shell: """
    bwa index -p mm39 -a bwtsw {input} 2> {log}
    mv mm39* 01_raw_data/
    """ 

rule align_reads:
    message: "Aligning paired end reads to GRCm39/mm39 reference genome"
    conda: "chip_seq_environment.yml"
    input:
        r1 = "01_raw_data/{sample}_1.fastq.gz",
        r2 = "01_raw_data/{sample}_2.fastq.gz"
    output: "03_sam_files/{sample}.sam"
    log: "00_logs/{sample}_align_reads_err.log"
    threads: 8
    shell: "bwa mem 01_raw_data/mm39 {input.r1} {input.r2} > {output} 2> {log}"
    
rule sam_to_bam:
    message: "Converting SAM to BAM file format"
    conda: "chip_seq_environment.yml"
    input: "03_sam_files/{sample}.sam"
    output: "04_bam_files/{sample}.bam"
    log: "00_logs/{sample}_sam_to_bam.log"
    shell: "samtools view -b {input} > {output} 2> 00_logs/{log}"

rule sam_fixmate:
    message: "Removing secondary and unmapped reads. Adding tags to reads for deduplication"
    conda: "chip_seq_environment.yml"
    input: "04_bam_files/{sample}.bam"
    output: "04_bam_files/{sample}.namesorted.fixmate.bam"
    log: "00_logs/{sample}_sam_fixmate.log"
    shell: "samtools fixmate -rcm -O bam {input} {output} 2> 00_logs/{log}"

rule sam_sort:
    message: "Sorting reads by chromosome coordinates"
    conda: "chip_seq_environment.yml"
    input: "04_bam_files/{sample}.namesorted.fixmate.bam"
    output: "04_bam_files/{sample}.coorsorted.fixmate.bam"
    log: "00_logs/{sample}_sam_sort.log"
    shell: "samtools sort {input} -o {output} 2> 00_logs/{log}"

rule sam_markdup:
    message: "Marking and removing duplicates"
    conda: "chip_seq_environment.yml"
    input: "04_bam_files/{sample}.coorsorted.fixmate.bam"
    output: "04_bam_files/{sample}.coorsorted.dedup.bam"
    log: "00_logs/{sample}_sam_markdup.log"
    shell: "samtools markdup -r --mode s {input} {output} 2> 00_logs/{log}"

rule sam_index:
    message: "Indexing deduplicated BAM file"
    conda: "chip_seq_environment.yml"
    input: "04_bam_files/{sample}.coorsorted.dedup.bam"
    output: "04_bam_files/{sample}.coorsorted.dedup.bam.bai", 
    log: "00_logs/{sample}_sam_index.log"
    shell: "samtools index {input} 2> 00_logs/{log}"

rule bam_to_bigwig:
    message: "Converting BAM file format to bigwig file format for visualization"
    conda: "chip_seq_environment.yml"
    input: "04_bam_files/{sample}.coorsorted.dedup.bam"
    output: "05_bigwig_files/{sample}.bw"
    log: "00_logs/{sample}_bam_to_bigwig.log"
    shell: "bamCoverage -b {input} -o {output} 2> 00_logs/{log}"

# Need to add {wildcard.samples} or have {sample} string name for -n option
# Need to run through with one sample to make sure outputs work
#rule call_peaks:
#    message: "Calling ChIP-seq peaks"
#    conda: "chip_seq_environment.yml"
#    input: "04_bam_files/{sample}.coorsorted.dedup.bam"
#    output: multiext('"06_macs2_peaks/{sample}", "_peaks.narrowPeak", "_peaks.xls", "_summits.bed", "_model.R", "_control_lambda.bdg", "_treat_pileup.bdg"')
#    params:
#        samp_name = "{sample}"
#    log: "00_logs/{sample}_macs2_peaks.log"
#    shell: "macs2 callpeak -t {input} -f BAM -n {params.samp_name} --outdir 06_macs2_peaks/ 2> 00_logs/{log}"