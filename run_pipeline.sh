#!/bin/bash

sudo chown $USER:$USER /mnt/d/Thesis/demo/ERR331/ERR9134331.fastq

# Quality control with FastQC
fastqc /mnt/d/Thesis/demo/ERR331/ERR9134331.fastq

# Step 1: Trimming
trimmomatic SE -threads 4 \
  /mnt/d/Thesis/demo/ERR331/ERR9134331.fastq \
  /mnt/d/Thesis/demo/ERR331/ERR9134331_trimmed.fastq \
  TRAILING:10 -phred33

fastqc /mnt/d/Thesis/demo/ERR331/ERR9134331_trimmed.fastq -o /mnt/d/Thesis/demo/ERR331

# Step 2: Alignment
hisat2 -q --rna-strandness R \
  -x /mnt/d/Thesis/HISAT2/grch38/grch38/genome \
  -U /mnt/d/Thesis/demo/ERR331/ERR9134331_trimmed.fastq \
  | samtools sort -o /mnt/d/Thesis/demo/ERR331/demo_aligned.bam

# Step 3: Feature counting
featureCounts -S 2 \
  -a /mnt/d/Thesis/hg38/Homo_sapiens.GRCh38.106.gtf \
  -o /mnt/d/Thesis/demo/ERR31/ERR9134331_featurecounts.txt \
  /mnt/d/Thesis/demo/ERR331/demo_aligned.bam
