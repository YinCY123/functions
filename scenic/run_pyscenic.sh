#!/bin/bash

run_pyscenic(){
    # usage/help
    echo " "
    local _usage="Usage: run_pyscenic -l <input.loom> [-t <tf_list.txt>] [-f <rankings.feather> ...] [-a <motif2tf.tbl>] [-o <output_dir>] [-w <num_workers>]\n\n\
    Description:\n\
    Run the three SCENIC steps (GRN, cisTarget, AUCell) on a Loom matrix.\n\n\
    Required:\n\
    -l  Path to input Loom expression matrix.\n\n\
    Optional (with defaults):\n\
    -t  Transcription factor list (txt, one TF per line).\n\
        Default: /home/yincy/BioHome/datasets/TF/Homo_sapiens_TF.txt\n\
    -f  cisTarget rankings feather file. Pass multiple -f flags for several ranking databases.\n\
        Default: /home/yincy/BioHome/scenic/human/hg38/mv_v10_clust/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather\n\
    -a  Motif-to-TF annotation table (.tbl).\n\
        Default: /home/yincy/BioHome/scenic/motif2tf_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl\n\
    -o  Output directory (created if missing). Default: ./scenic_out\n\
    -n  Number of workers for PySCENIC. Default: 30. Do not using to many threads, this will cause error.\n\
    -h  Show this help and exit.\n\n\
    Example:\n\
    run_pyscenic -l data/input.loom -o results/scenic \n"

    # initialize variables (defaults)
    local loom
    local output_dir=./scenic_out
    local transcription_factor=/home/yincy/BioHome/datasets/TF/Homo_sapiens_TF.txt
    local feather_default=/home/yincy/BioHome/scenic/human/hg38/mv_v10_clust/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
    local -a feather_files=("$feather_default")
    local tbl=/home/yincy/BioHome/scenic/motif2tf_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
    local num_workers=30 # do not use to many threads, this with cause error

    # parse args
    while (( "$#" )); do
        case "$1" in
            -l|--loom)
                loom="$2"; shift 2 ;;
            -t|--transcription_factor)
                transcription_factor="$2"; shift 2 ;;
            -f|--feather)
                local IFS=','
                local -a _feathers
                read -r -a _feathers <<< "$2"
                if [[ ${#feather_files[@]} -eq 1 && "${feather_files[0]}" == "$feather_default" ]]; then
                    feather_files=()
                fi
                feather_files+=("${_feathers[@]}")
                shift 2 ;;
            -a|--tbl)
                tbl="$2"; shift 2 ;;
            -o|--output_dir)
                output_dir="$2"; shift 2 ;;
            --num_workers|-n)
                num_workers="$2"; shift 2 ;;
            -h|--help)
                echo -e "$_usage"; return 0 ;;
            --)
                shift; break ;;
            -*)
                echo "Unknown option: $1" >&2; echo -e "$_usage"; return 1 ;;
            *)
                break ;;
        esac
    done

    # validate loom file exist
    if [[ -z "$loom" ]]; then
        echo "Error: -l <input.loom> is required." >&2
        echo -e "$_usage"
        return 1
    fi

    # ensure output dir exists
    mkdir -p "$output_dir"


    # 01 grn
    echo "run step 1: grn"
    pyscenic grn \
            "$loom" \
            "$transcription_factor" \
            --num_workers "$num_workers" \
            --output ${output_dir}/01_grn_out.tsv \
            --method grnboost2

    # 02 cistarget

    echo "run step 2 ctx"
    if [[ ${#feather_files[@]} -eq 0 ]]; then
        echo "Error: at least one -f <rankings.feather> must be provided." >&2
        return 1
    fi
    pyscenic ctx \
            ${output_dir}/01_grn_out.tsv "${feather_files[@]}" \
            --annotations_fname "$tbl" \
            --expression_mtx_fname "$loom" \
            --mode dask_multiprocessing \
            --output ${output_dir}/02_cistarget_out.csv \
            --num_workers "$num_workers" \
            --mask_dropouts


    # 03 AUCell
    echo "run step 3 AUCell"
    pyscenic aucell \
        "$loom" \
        ${output_dir}/02_cistarget_out.csv \
        --output ${output_dir}/03_aucell_out.csv \
        --num_workers "$num_workers"

    # 04 Binarize AUCell matrix
    echo "run step 4 binarize AUCell"
    python << EOF
import pandas as pd
from pyscenic.binarization import binarize

# Load AUC matrix
auc_df = pd.read_csv('${output_dir}/03_aucell_out.csv', index_col=0)

# Binarize using automatic threshold detection based on bimodal distribution
binarized_df, thresholds = binarize(auc_df)

# Save binarized matrix
binarized_df.to_csv('${output_dir}/04_aucell_binarized_out.csv')

# Save thresholds used for each regulon
thresholds.to_csv('${output_dir}/04_aucell_thresholds.csv', header=['threshold'])

print(f'Binarization completed. Used thresholds saved to ${output_dir}/04_aucell_thresholds.csv')
EOF

    echo "SCENIC analysis is done!"
    return 0
}
