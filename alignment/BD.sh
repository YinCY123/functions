#!/bin/bash

genomeDir="/home/yincy/BioHome/index/STAR_GRCh38"
seqDir="fastq"
cls_dir="/home/yincy/BioHome/github/functions/alignment/STARsolo_DB"
gtf="/home/yincy/BioHome/reference/human/gencode.v49.primary_assembly.annotation.gtf"


STAR \
        --outFileNamePrefix "${seqDir}/SRR11666955" \
        --genomeDir "${genomeDir}" \
        --runThreadN 100 \
        --readFilesIn fastq/SRR11666955_R1.fastq.gz fastq/SRR11666955_R2.fastq.gz \
        --readFilesCommand zcat \
        --sjdbGTFfile "$gtf" \
        --soloType CB_UMI_Complex \
        --soloCBmatchWLtype 1MM \
        --soloCBwhitelist "${cls_dir}/BD_CLS1.txt" "${cls_dir}/BD_CLS2.txt" "${cls_dir}/BD_CLS3.txt" \
        --soloUMIlen 8 \
        --soloCBposition 0_0_0_8 0_21_0_29 0_43_0_51 \
        --soloUMIposition 0_52_0_59 \
        --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
