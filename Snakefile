configfile: "samples.yaml"

print("Starting ChIP-seq data analysis workflow")
            
rule all:
    input: 
        "01_raw_data/mm39.amb",
        "01_raw_data/mm39.ann",
        "01_raw_data/mm39.bwt",
        "01_raw_data/mm39.pac",
        "01_raw_data/mm39.sa", 
        expand("02_fastqc_analysis/{sample}_1_fastqc.html", sample=config["sample"]),
        expand("02_fastqc_analysis/{sample}_1_fastqc.zip", sample=config["sample"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.html", sample=config["sample"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.zip", sample=config["sample"]),
        expand("04_bam_files/{sample}.coorsorted.dedup.bam.bai", sample=config["sample"]),
        expand("05_bigwig_files/{sample}.bw", sample=config["sample"])

rule make_directories:
    message: "Making directories for data organization"
    output:
        directory("00_logs/"),
        directory("01_raw_data/"),
        directory("02_fastqc_analysis/"),
        directory("03_sam_files/"),
        directory("04_bam_files/"),
        directory("05_bigwig_files/")
    shell: """
        mkdir 00_logs
        mkdir 01_raw_data
        mkdir 02_fastqc_analysis
        mkdir 03_sam_files
        mkdir 04_bam_files
        mkdir 05_bigwig_files
    """

rule download_data:
    message: "Downloading raw data files"
    conda: "chip_seq_environment.yml"
    #output: expand("01_raw_data/{sample}/{sample}.sra", sample=config["sample"])
    shell: "prefetch {config[sample]}"

rule split_paired_reads:
    message: "Splitting paired end reads into separate files"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}/{sample}.sra", sample=config["sample"])
    #output:
        #expand("01_raw_data/{sample}_1.fastq.gz", sample=config["sample"]),
        #expand("01_raw_data/{sample}_2.fastq.gz", sample=config["sample"])
    shell: "fastq-dump {input} --split-files --gzip --outdir 01_raw_data/"
    
rule fastqc_precheck_r1:
    message: "Running quality control on samples pre-processing"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}_1.fastq.gz", sample=config["sample"]),
    output:
        expand("02_fastqc_analysis/{sample}_1_fastqc.html", sample=config["sample"]),
        expand("02_fastqc_analysis/{sample}_1_fastqc.zip", sample=config["sample"]),
    shell: "fastqc {input} --outdir 02_fastqc_analysis/"

rule fastqc_precheck_r2:
    message: "Running quality control on samples pre-processing"
    conda: "chip_seq_environment.yml"
    input: expand("01_raw_data/{sample}_2.fastq.gz", sample=config["sample"]),
    output:
        expand("02_fastqc_analysis/{sample}_2_fastqc.html", sample=config["sample"]),
        expand("02_fastqc_analysis/{sample}_2_fastqc.zip", sample=config["sample"]),
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
    output: multiext("01_raw_data/mm39", ".amb", ".ann", ".bwt", ".pac", ".sa")
    shell: """
    bwa index -p mm39 -a bwtsw {input}
    mv mm39* 01_raw_data/
    """ 

rule align_reads:
    message: "Aligning paired end reads to GRCm39/mm39 reference genome"
    conda: "chip_seq_environment.yml"
    input:
        r1 = expand("01_raw_data/{sample}_1.fastq.gz", sample=config["sample"]),
        r2 = expand("01_raw_data/{sample}_2.fastq.gz", sample=config["sample"])
    output: expand("03_sam_files/{sample}.sam", sample=config["sample"])
    log: expand("00_logs/{sample}_align_reads_err.log", sample=config["sample"])
    shell: "bwa mem 01_raw_data/mm39 {input.r1} {input.r2} > {output} 2> {log}"
    
rule sam_to_bam:
    message: "Converting SAM to BAM file format"
    conda: "chip_seq_environment.yml"
    input: expand("03_sam_files/{sample}.sam", sample=config["sample"])
    output: expand("04_bam_files/{sample}.bam", sample=config["sample"])
    shell: "samtools view -b {input} > {output}"

rule sam_fixmate:
    message: "Removing secondary and unmapped reads. Adding tags to reads for deduplication"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.bam", sample=config["sample"])
    output: expand("04_bam_files/{sample}.namesorted.fixmate.bam", sample=config["sample"])
    shell: "samtools fixmate -rcm -O bam {input} {output}"

rule sam_sort:
    message: "Sorting reads by chromosome coordinates"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.namesorted.fixmate.bam", sample=config["sample"])
    output: expand("04_bam_files/{sample}.coorsorted.fixmate.bam", sample=config["sample"])
    shell: "samtools sort {input} -o {output}"

rule sam_markdup:
    message: "Marking and removing duplicates"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.coorsorted.fixmate.bam", sample=config["sample"])
    output: expand("04_bam_files/{sample}.coorsorted.dedup.bam", sample=config["sample"])
    shell: "samtools markdup -r --mode s {input} {output}"

rule sam_index:
    message: "Indexing deduplicated BAM file"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.coorsorted.dedup.bam", sample=config["sample"])
    output: expand("04_bam_files/{sample}.coorsorted.dedup.bam.bai", sample=config["sample"]), 
    shell: """
    samtools index {input}
    """

rule bam_to_bigwig:
    message: "Converting BAM file format to bigwig file format for visualization"
    conda: "chip_seq_environment.yml"
    input: expand("04_bam_files/{sample}.coorsorted.dedup.bam", sample=config["sample"])
    output: expand("05_bigwig_files/{sample}.bw", sample=config["sample"])
    shell: "bamCoverage -b {input} -o {output}"
