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
        expand("03_bam_files/{sample}.coorsorted.dedup.bam.bai", sample=SAMPLES),
        expand("04_bigwig_files/{sample}.bw", sample=SAMPLES)

rule organize_data:
    message: "Making directories for data organization"
    output:
        directory("01_raw_data/"),
        directory("02_sam_files/"),
        directory("03_bam_files/"),
        directory("04_bigwig_files/")
    shell: """
        mkdir 01_raw_data
        mkdir 02_sam_files
        mkdir 03_bam_files
        mkdir 04_bigwig_files
    """

rule download_data:
    message: "Downloading raw data files"
    conda: "chip_seq_environment.yml"
    output: directory(expand("01_raw_data/{sample}/", sample=SAMPLES))
    shell: """
        for i in $( grep -v "^#" samples.txt ); do
            prefetch $i
            mv $i/ 01_raw_data/$i/
        done
    """

rule split_paired_reads:
    message: "Splitting paired end reads into separate files"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}/{sample}.sra", sample=SAMPLES)
    output:
        expand("{sample}_1.fastq.gz", sample=SAMPLES),
        expand("{sample}_2.fastq.gz", sample=SAMPLES)
    shell: "fastq-dump {input} --split-files --gzip"

rule download_genome:
    message: "Downloading GRCm39/mm39 mouse genome from the UCSC Genome Browser"
    output: "mm39.chromFa.tar.gz"
    shell: "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz -O {output}"
    
rule decompress_genome:
    message: "Decompressing genome"
    input: "mm39.chromFa.tar.gz"
    shell: "tar zvfx {input}"

rule concatenate_chromosomes:
    message: "Concatenating individual chromosome files to create full assembly"
    output: "mm39.fa"
    shell: "cat *.fa > {output}" 

rule delete_chromosome_files:
    message: "Removing chromosome sequence files"
    shell: "rm chr*.fa"  
    
rule set_alignment_reference:
    message: "Setting GRCm39/mm39 mouse genome assembly as reference genome for alignment" 
    conda: "chip_seq_environment.yml"
    input: "mm39.fa"
    output:
        "mm39.amb",
        "mm39.ann",
        "mm39.bwt",
        "mm39.pac",
        "mm39.sa"
    shell: "bwa index -p mm39 -a bwtsw {input}" 

rule align_reads:
    message: "Aligned paired end reads to GRCm39/mm39 reference genome"
    conda: "chip_seq_environment.yml"
    input:
        r1 = expand("{sample}_1.fastq.gz", sample=SAMPLES),
        r2 = expand("02_sam_files/{sample}_2.fastq.gz", sample=SAMPLES)
    output: expand("{sample}.sam", sample=SAMPLES)
    shell: "bwa mem mm39 {input.r1} {input.r2} > {output}"
    
rule sam_to_bam:
    message: "Converting SAM to BAM file format"
    conda: "chip_seq_environment.yml"
    input: expand("{sample}.sam", sample=SAMPLES)
    output: expand("{sample}.bam", sample=SAMPLES)
    shell: "samtools view -b {input} > {output}"

rule sam_fixmate:
    message: "Removing secondary and unmapped reads. Adding tags to reads for deduplication"
    conda: "chip_seq_environment.yml"
    input: expand("{sample}.bam", sample = SAMPLES)
    output: expand("{sample}.namesorted.fixmate.bam", sample=SAMPLES)
    shell: "samtools fixmate -rcm -O bam {input} {output}"

rule sam_sort:
    message: "Sorting reads by chromosome coordinates"
    conda: "chip_seq_environment.yml"
    input: expand("{sample}.namesorted.fixmate.bam", sample=SAMPLES)
    output: expand("{sample}.coorsorted.fixmate.bam", sample=SAMPLES)
    shell: "samtools sort {input} -o {output}"

rule sam_markdup:
    message: "Marking and removing duplicates"
    conda: "chip_seq_environment.yml"
    input: expand("{sample}.coorsorted.fixmate.bam", sample=SAMPLES)
    output: expand("{sample}.coorsorted.dedup.bam", sample=SAMPLES)
    shell: "samtools markdup -r --mode s {input} {output}"

rule sam_index:
    message: "Indexing deduplicated BAM file"
    conda: "chip_seq_environment.yml"
    input: expand("{sample}.coorsorted.dedup.bam", sample=SAMPLES)
    output: expand("{sample}.coorsorted.dedup.bam.bai", sample=SAMPLES), 
    shell: "samtools index {input}"

rule bam_to_bigwig:
    message: "Converting BAM file format to bigwig file format for visualization"
    conda: "chip_seq_environment.yml"
    input: expand("{sample}.coorsorted.dedup.bam", sample=SAMPLES)
    output: expand("{sample}.bw", sample=SAMPLES)
    shell: "bamCoverage -b {input} -o {output}"

rule organize_data_1:
    message: "Organizing data and output files"
    run:
        with open("samples.txt", "r") as a_file:
            for line in a_file:
                if not line.lstrip().startswith('#'):
                    # Moving .sra files out of SRA folders
                    os.system(f"mv {line[:-1]}/{line[:-1]}.sra 01_raw_data/{line[:-1]}.sra")
                    # Deleting empty SRA folder
                    os.system(f"rm -rf {line[:-1]}/")
                    # Moving reference genome files to 01_raw_data directory
                    os.system("mv mm39.amb 01_raw_data/mm39.amb")
                    os.system("mv mm39.bwt 01_raw_data/mm39.bwt")
                    os.system("mv mm39.fa 01_raw_data/mm39.fa")
                    os.system("mv mm39.sa 01_raw_data/mm39.sa")
                    os.system("mv mm39.ann 01_raw_data/mm39.ann")
                    os.system("mv mm39.pac 01_raw_data/mm39.pac")
                    # Deleting chromosome file folder for reference genome
                    os.system("rm mm39.chromFa.tar.gz")
                    # Moving fastq files to 01_raw_data directory
                    os.system(f"mv {line[:-1]}_1.fastq.gz 01_raw_data/{line[:-1]}_1.fastq.gz")
                    os.system(f"mv {line[:-1]}_2.fastq.gz 01_raw_data/{line[:-1]}_2.fastq.gz")
                    # Moving SAM files to 02_sam_files directory
                    os.system(f"mv {line[:-1]}.sam 02_sam_files/{line[:-1]}.sam")
                    # Moving BAM files to 03_bam_files directory
                    os.system(f"mv {line[:-1]}.bam 03_bam_files/{line[:-1]}.bam")
                    os.system(f"mv {line[:-1]}.namesorted.fixmate.bam 03_bam_files/{line[:-1]}.namesorted.fixmate.bam")
                    os.system(f"mv {line[:-1]}.coorsorted.fixmate.bam 03_bam_files/{line[:-1]}.coorsorted.fixmate.bam")
                    os.system(f"mv {line[:-1]}.coorsorted.dedup.bam 03_bam_files/{line[:-1]}.coorsorted.dedup.bam")
                    os.system(f"mv {line[:-1]}.coorsorted.dedup.bam.bai 03_bam_files/{line[:-1]}.coorsorted.dedup.bam.bai")
                    # Moving bigwig files to 04_bigwig_files directory
                    os.system(f"mv {line[:-1]}.bw 04_bigwig_files/{line[:-1]}.bw")