#!/bin/bash

feather="/home/yincy/BioHome/scenic/human/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather"
tbl="/home/yincy/BioHome/scenic/motif2tf_annotation/motifs-v9-nr.hgnc-m0.001-o0.0.tbl"
TF="/home/yincy/BioHome/datasets/TF/Homo_sapiens_TF.txt"

loom="/home/yincy/BioHome/projects/cns_250806/results/data/scenic/macrophage_MDSC_monocyte_mtx.loom"
output_dir="results/data/scenic"

# 01 grn
echo "run step 1: grn"
pyscenic grn \
        $loom \
        $TF \
        --num_workers 50 \
        --output ${output_dir}/grn_out.tsv \
        --method grnboost2

# 02 cistarget

echo "run step 2 ctx"
pyscenic ctx \
        ${output_dir}/grn_out.tsv $feather \
        --annotations_fname $tbl \
        --expression_mtx_fname $loom \
        --mode dask_multiprocessing \
        --output ${output_dir}/cistarget_out.csv \
        --num_workers 50 \
        --mask_dropouts


# 03 AUCell
echo "run step 3 AUCell"
pyscenic aucell \
      $loom \
      ${output_dir}/cistarget_out.csv \
      --output ${output_dir}/aucell_out.csv \
      --num_workers 50 
