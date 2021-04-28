import os

SAMPLES = []

with open("samples.txt", "r") as a_file:
    for line in a_file:
        if not line.lstrip().startswith('#'):
            SAMPLES.append(line[:-1])

#rule all:
#    input: 
#        expand("{sample}.bw", sample=SAMPLES)

# Having trouble with output formatting, but data gets downloaded
rule download_data:
    message: "Downloading raw data files"
    #output: expand("{sample}/", sample=SAMPLES)
    run:
        for sample in SAMPLES:
                os.system(f"prefetch {sample}")

rule split_paired_reads:
    message: "Splitting paired end reads into separate files"
    input: expand("{sample}/{sample}.sra", sample=SAMPLES)
    output:
        protected(expand("{sample}_1.fastq.gz", sample=SAMPLES)),
        protected(expand("{sample}_2.fastq.gz", sample=SAMPLES))
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
    output: protected("mm39.fa")
    shell: "cat *.fa > {output}" 

rule delete_chromosome_files:
    message: "Removing chromosome sequence files"
    shell: "rm chr*.fa"  
    
rule set_alignment_reference:
    message: "Setting GRCm39/mm39 mouse genome assembly as reference genome for alignment" 
    input: "mm39.fa"
    output:
        protected("mm39.amb"),
        protected("mm39.ann"),
        protected("mm39.bwt"),
        protected("mm39.pac"),
        protected("mm39.sa")
    shell: "bwa index -p mm39 -a bwtsw {input}" 

rule align_reads:
    message: "Aligned paired end reads to GRCm39/mm39 reference genome"
    input:
        r1 = expand("{sample}_1.fastq.gz", sample=SAMPLES),
        r2 = expand("{sample}_2.fastq.gz", sample=SAMPLES)
    output: expand("{sample}.sam", sample=SAMPLES)
    shell: "bwa mem mm39 {input.r1} {input.r2} > {output}"
    
rule sam_to_bam:
    message: "Converting SAM to BAM file format"
    input: expand("{sample}.sam", sample=SAMPLES)
    output: expand("{sample}.bam", sample=SAMPLES)
    shell: "samtools view -b {input} > {output}"

rule sam_fixmate:
    message: "Removing secondary and unmapped reads. Adding tags to reads for deduplication"
    input: expand("{sample}.bam", sample = SAMPLES)
    output: expand("{sample}.namesorted.fixmate.bam", sample=SAMPLES)
    shell: "samtools fixmate -rcm -O bam {input} {output}"

rule sam_sort:
    message: "Sorting reads by chromosome coordinates"
    input: expand("{sample}.namesorted.fixmate.bam", sample=SAMPLES)
    output: expand("{sample}.coorsorted.fixmate.bam ", sample=SAMPLES)
    shell: "samtools sort -o {output} {input}"

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

rule organize_data:
    message: "Organizing data and output files"
    run:
        # Make 01_raw_data directory
        os.system("mkdir 01_raw_data")
        with open("samples.txt", "r") as a_file:
            for line in a_file:
                if not line.lstrip().startswith('#'):
                    # Moving .sra files out of SRA folders
                    os.system(f"mv {line[:-1]}/{line[:-1]}.sra 01_raw_data/{line[:-1]}.sra")
                    # Deleting empty SRA folder
                    os.system(f"rm -rf {line[:-1]}/")
