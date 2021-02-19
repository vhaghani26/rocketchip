# Deduplication using samtools
samtools fixmate -rcm -O bam /share/korflab/project/background/02.2_bwa_mem_bam/{file}.bam /share/korflab/project/background/03_dedup/{file}.namesorted.fixmate.bam 2> /share/korflab/project/background/00_std_errs/{file}_fixmate.log
samtools sort -o /share/korflab/project/background/03_dedup/{file}.coorsorted.fixmate.bam /share/korflab/project/background/03_dedup/{file}.namesorted.fixmate.bam 2> /share/korflab/project/background/00_std_errs/{file}_samsort.log
samtools markdup -r --mode s /share/korflab/project/background/03_dedup/{file}.coorsorted.fixmate.bam /share/korflab/project/background/03_dedup/{file}.coorsorted.dedup.bam 2> /share/korflab/project/background/00_std_errs/{file}_markdup.log
samtools index /share/korflab/project/background/03_dedup/{file}.coorsorted.dedup.bam 2> /share/korflab/project/background/00_std_errs/{file}_index.log

# Converting deduplicated BAM file to bigWig file for visualization in UCSC GB
bamCoverage -b /share/korflab/project/background/03_dedup/{file}.coorsorted.dedup.bam -o /share/korflab/project/background/04_visual_ucscgb/{file}.bw 2> /share/korflab/project/background/00_std_errs/{file}_bamcoverage.log

# Create track lines for data
# Go to https://genome.ucsc.edu/s/vhaghani/IgG_controls_chipseq to view tracks
# Open IgG_controls_chipseq_tracklines.txt to see the tracklines I used for visualizing the data 






######################################################################################################

rule download_data:
    conda: ############### SRA_toolkit
    output: "{sample}.sra"
    shell: "prefetch {SRA code}" ######

rule split_paired_reads:
    conda: ###### SRA_toolkit
    input: "{sample}.sra"
    output:
    shell: ##### "fastq-dump ~/Code/Rotation_2_Project/SRA_Toolkit/sra/{SRA Code}.sra --split-files --outdir ../files"

rule gzip_data:
    input: 
    output:
    shell: #### "gzip files/SRR*"

rule download_genome:
    output: "mm39.chromFa.tar.gz"
    shell: "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz"
    
rule decompress_genome:
    input: "mm39.chromFa.tar.gz"
    output: #######
    shell: "tar zvfx {input}"

rule concatenate_chromosomes:
    input: ######
    output: "mm39.fa"
    shell: "cat *.fa > {output}" #######

rule delete_chromosome_files:
    input:
    output:
    shell: "rm chr*.fa   " ###### 
    
rule set_alignment_reference:
    conda: ############ BWA
    input: "mm39.fa"
    output:
    shell: "bwa index -p mm39 -a bwtsw {input}" 

rule align_reads:
    conda: ###### BWA
    input:
        r1 = "{sample}_1.fastq.gz",
        r2 = "{sample}_2.fastq.gz"
    output: "{sample}.sam"
    shell: "bwa mem mm39 {sample}_1.fastq.gz {sample}_2.fastq.gz > {sample}.sam"
    
rule sam_to_bam:
    conda: ############samtools
    input: "{sample}.sam"
    output: "{sample}.bam"
    shell: "samtools view -b {input} > {output}"

    
