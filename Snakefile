import os

#rule all:
#    input: 
#        expand("{sample}.sra", sample=SAMPLES)

rule download_data:
    message: "Downloading raw data files"
    output: protected("01_raw_data/{sample}.sra")
    run:
        with open("samples.txt", "r") as a_file:
            for line in a_file:
                if not line.lstrip().startswith('#'):
                    os.system(f"echo {line} > {output}")

#rule split_paired_reads:
#    input: "{sample}.sra"
#    output:
#    shell: "fastq-dump {sample}.sra --split-files --outdir ../files"

#rule gzip_data:
#    input: 
#    output:
#    shell: "gzip files/SRR*"

rule download_genome:
    message: "Downloading GRCm39/mm39 mouse genome from the UCSC Genome Browser"
    output: "01_raw_data/mm39.chromFa.tar.gz"
    shell: "wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.chromFa.tar.gz -O {output}"
    
rule decompress_genome:
    input: "01_raw_data/mm39.chromFa.tar.gz"
    shell: "tar zvfx {input}"

#rule concatenate_chromosomes:
#    input: 
#    output: protected("mm39.fa")
#    shell: "cat *.fa > {output}" 

#rule delete_chromosome_files:
#    input:
#    output:
#    shell: "rm chr*.fa   "  
    
#rule set_alignment_reference:
#    input: "mm39.fa"
#    output:
#    shell: "bwa index -p mm39 -a bwtsw {input}" 

#rule align_reads:
#    input:
#        r1 = "{sample}_1.fastq.gz",
#        r2 = "{sample}_2.fastq.gz"
#    output: "{sample}.sam"
#    shell: "bwa mem mm39 {sample}_1.fastq.gz {sample}_2.fastq.gz > {sample}.sam"
    
#rule sam_to_bam:
#    input: "{sample}.sam"
#    output: "{sample}.bam"
#    shell: "samtools view -b {input} > {output}"

#rule sam_fixmate:
#    input: "{sample}.bam"
#    output: "{sample}.namesorted.fixmate.bam"
#    shell: "samtools fixmate -rcm -O bam {input} {output}"

#rule sam_sort:
#    input: "{sample}.namesorted.fixmate.bam"
#    output: "{sample}.coorsorted.fixmate.bam "
#    shell: "samtools sort -o {output} {input}"

#rule sam_markdup:
#    input: "{sample}.coorsorted.fixmate.bam"
#    output: "{sample}.coorsorted.dedup.bam"
#    shell: "samtools markdup -r --mode s {input} {output}"

#rule sam_index:
#    input: "{sample}.coorsorted.dedup.bam"
#    output: "{sample}.indexed.dedup.bam"
#    shell: "samtools index {input}"

#rule bam_to_bigwig:
#    input: "{sample}.indexed.dedup.bam"
#    output: "{sample}.bw"
#    shell: "bamCoverage -b {input} -o {output}"
