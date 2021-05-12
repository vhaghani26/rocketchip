import os

SAMPLES = []

with open("samples.txt", "r") as a_file:
    for line in a_file:
        if not line.lstrip().startswith('#'):
            SAMPLES.append(line[:-1])
            
rule all:
    input: 
        "01_raw_data/mm39.amb",
        "01_raw_data/mm39.ann",
        "01_raw_data/mm39.bwt",
        "01_raw_data/mm39.pac",
        "01_raw_data/mm39.sa", 
        expand("02_fastqc_analysis/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("02_fastqc_analysis/{sample}_1_fastqc.zip", sample=SAMPLES),
        expand("02_fastqc_analysis/{sample}_2_fastqc.html", sample=SAMPLES),
        expand("02_fastqc_analysis/{sample}_2_fastqc.zip", sample=SAMPLES),
        expand("04_bam_files/{sample}.coorsorted.dedup.bam.bai", sample=SAMPLES),
        expand("05_bigwig_files/{sample}.bw", sample=SAMPLES)

rule make_directories:
    message: "Making directories for data organization"
    output:
        directory("01_raw_data/"),
        directory("02_fastqc_analysis/"),
        directory("03_sam_files/"),
        directory("04_bam_files/"),
        directory("05_bigwig_files/"),

    shell: """
        mkdir 01_raw_data
        mkdir 02_fastqc_analysis
        mkdir 03_sam_files
        mkdir 04_bam_files
        mkdir 05_bigwig_files
    """

rule download_data:
    message: "Downloading raw data files"
    conda: "chip_seq_environment.yml"
    output: expand("01_raw_data/{sample}/{sample}.sra", sample=SAMPLES)
    shell: """
        for i in $( grep -v "^#" samples.txt ); do
            prefetch $i
            mv $i/ 01_raw_data/
        done
    """

rule split_paired_reads:
    message: "Splitting paired end reads into separate files"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}/{sample}.sra", sample=SAMPLES)
    output:
        expand("01_raw_data/{sample}_1.fastq.gz", sample=SAMPLES),
        expand("01_raw_data/{sample}_2.fastq.gz", sample=SAMPLES)
    shell: "fastq-dump {input} --split-files --gzip --outdir 01_raw_data/"
    
rule fastqc_precheck_r1:
    message: "Running quality control on samples pre-processing"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}_1.fastq.gz", sample=SAMPLES),
    output:
        expand("02_fastqc_analysis/{sample}_1_fastqc.html", sample=SAMPLES),
        expand("02_fastqc_analysis/{sample}_1_fastqc.zip", sample=SAMPLES),
    shell: "fastqc {input} --outdir 02_fastqc_analysis/"

rule fastqc_precheck_r2:
    message: "Running quality control on samples pre-processing"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}_2.fastq.gz", sample=SAMPLES),
    output:
        expand("02_fastqc_analysis/{sample}_2_fastqc.html", sample=SAMPLES),
        expand("02_fastqc_analysis/{sample}_2_fastqc.zip", sample=SAMPLES),
    shell: "fastqc {input} --outdir 02_fastqc_analysis/"   

rule download_genome:
    message: "Downloading GRCm39/mm39 mouse genome from the UCSC Genome Browser"
    output: "01_raw_data/mm39.chromFa.tar.gz"
    shell: "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz -O {output}"

rule process_genome:
    message: "Decompressing genome. Concatenating individual chromosome files to create full assembly. Removing chromosome sequence files"
    input: "01_raw_data/mm39.chromFa.tar.gz"
    output: "01_raw_data/mm39.fa"
    shell: """
    tar zvfx {input} --directory 01_raw_data/
    cat 01_raw_data/*.fa > {output}
    rm 01_raw_data/chr*.fa
    """
    
rule set_alignment_reference:
    message: "Setting GRCm39/mm39 mouse genome assembly as reference genome for alignment" 
    conda: "chip_seq_environment.yml"
    input: "01_raw_data/mm39.fa"
    output:
        "01_raw_data/mm39.amb",
        "01_raw_data/mm39.ann",
        "01_raw_data/mm39.bwt",
        "01_raw_data/mm39.pac",
        "01_raw_data/mm39.sa"
    shell: """
    bwa index -p mm39 -a bwtsw {input}
    mv mm39* 01_raw_data/
    """ 

rule align_reads:
    message: "Aligned paired end reads to GRCm39/mm39 reference genome"
    conda: "chip_seq_environment.yml"
    input:
        r1 = expand("01_raw_data/{sample}_1.fastq.gz", sample=SAMPLES),
        r2 = expand("01_raw_data/{sample}_2.fastq.gz", sample=SAMPLES)
    output: expand("03_sam_files/{sample}.sam", sample=SAMPLES)
    shell: "bwa mem 01_raw_data/mm39 {input.r1} {input.r2} > {output}"
    
rule sam_to_bam:
    message: "Converting SAM to BAM file format"
    conda: "chip_seq_environment.yml"
    input: expand("03_sam_files/{sample}.sam", sample=SAMPLES)
    output: expand("04_bam_files/{sample}.bam", sample=SAMPLES)
    shell: "samtools view -b {input} > {output}"

rule sam_fixmate:
    message: "Removing secondary and unmapped reads. Adding tags to reads for deduplication"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.bam", sample = SAMPLES)
    output: expand("04_bam_files/{sample}.namesorted.fixmate.bam", sample=SAMPLES)
    shell: "samtools fixmate -rcm -O bam {input} {output}"

rule sam_sort:
    message: "Sorting reads by chromosome coordinates"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.namesorted.fixmate.bam", sample=SAMPLES)
    output: expand("04_bam_files/{sample}.coorsorted.fixmate.bam", sample=SAMPLES)
    shell: "samtools sort {input} -o {output}"

rule sam_markdup:
    message: "Marking and removing duplicates"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.coorsorted.fixmate.bam", sample=SAMPLES)
    output: expand("04_bam_files/{sample}.coorsorted.dedup.bam", sample=SAMPLES)
    shell: "samtools markdup -r --mode s {input} {output}"

rule sam_index:
    message: "Indexing deduplicated BAM file"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.coorsorted.dedup.bam", sample=SAMPLES)
    output: expand("04_bam_files/{sample}.coorsorted.dedup.bam.bai", sample=SAMPLES), 
    shell: """
    samtools index {input}
    mv *.bai 04_bam_files/
    """

rule bam_to_bigwig:
    message: "Converting BAM file format to bigwig file format for visualization"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.coorsorted.dedup.bam", sample=SAMPLES)
    output: expand("05_bigwig_files/{sample}.bw", sample=SAMPLES)
    shell: "bamCoverage -b {input} -o {output}"
