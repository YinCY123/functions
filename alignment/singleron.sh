#!/bin/bash
genome="/home/yincy/BioHome/index/STAR_GRCh38"
fastq_dir=fastq
output_dir="singleron_out/"
white="/home/yincy/BioHome/github/celescope-mobiu/celescope/data/chemistry/mobiu-1"


STAR --runThreadN 100 \
     --genomeDir $genome \
     --readFilesIn ${fastq_dir}/SRR11666954_R1.fastq.gz ${fastq_dir}/SRR11666954_R2.fastq.gz \
     --readFilesCommand zcat \
     --soloType CB_UMI_Complex \
     --soloCBposition 0_0_0_8 0_15_0_23 0_30_0_38 \
     --soloUMIposition 0_40_0_51 \
     --soloBarcodeReadLength 52 \
     --soloCBmatchWLtype 1MM \
     --soloFeatures Gene GeneFull_Ex50pAS \
     --soloCellFilter EmptyDrops_CR \
     --outFilterMatchNmin 50 \
     --soloCBwhitelist ${white}/bc1.txt ${white}/bc2.txt ${white}/bc3.txt \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
     --outFileNamePrefix $output_dir
