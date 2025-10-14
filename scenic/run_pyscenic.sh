#!/bin/bash

run_pyscenic(){
    # usage/help
    local _usage="Usage: run_pyscenic -l <input.loom> [-t <tf_list.txt>] [-f <rankings.feather>] [-a <motif2tf.tbl>] [-o <output_dir>] [-w <num_workers>]\n\n\
    Description:\n\
    Run the three SCENIC steps (GRN, cisTarget, AUCell) on a Loom matrix.\n\n\
    Required:\n\
    -l  Path to input Loom expression matrix.\n\n\
    Optional (with defaults):\n\
    -t  Transcription factor list (txt, one TF per line).\n\
        Default: /home/yincy/BioHome/datasets/TF/Homo_sapiens_TF.txt\n\
    -f  cisTarget rankings feather file.\n\
        Default: /home/yincy/BioHome/scenic/human/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather\n\
    -a  Motif-to-TF annotation table (.tbl).\n\
        Default: /home/yincy/BioHome/scenic/motif2tf_annotation/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl\n\
    -o  Output directory (created if missing). Default: ./scenic_out\n\
    -w  Number of workers for PySCENIC. Default: 50\n\
    -h  Show this help and exit.\n\n\
    Example:\n\
    run_pyscenic -l data/input.loom -o results/scenic \n"

    # initialize variables (defaults)
    local loom
    local output_dir=./scenic_out
    local transcription_factor=/home/yincy/BioHome/datasets/TF/Homo_sapiens_TF.txt
    local feather=/home/yincy/BioHome/scenic/human/hg38/mc9nr/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
    local tbl=/home/yincy/BioHome/scenic/motif2tf_annotation/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
    local num_workers=50

    # parse args
    while (( "$#" )); do
        case "$1" in
            -l|--loom)
                loom="$2"; shift 2 ;;
            -t|--transcription_factor)
                transcription_factor="$2"; shift 2 ;;
            -f|--feather)
                feather="$2"; shift 2 ;;
            -a|--annotations|-b|--tbl)
                tbl="$2"; shift 2 ;;
            -o|--output_dir)
                output_dir="$2"; shift 2 ;;
            -w|--num_workers|-n)
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
            --output ${output_dir}/grn_out.tsv \
            --method grnboost2

    # 02 cistarget

    echo "run step 2 ctx"
    pyscenic ctx \
            ${output_dir}/grn_out.tsv "$feather" \
            --annotations_fname "$tbl" \
            --expression_mtx_fname "$loom" \
            --mode dask_multiprocessing \
            --output ${output_dir}/cistarget_out.csv \
            --num_workers "$num_workers" \
            --mask_dropouts


    # 03 AUCell
    echo "run step 3 AUCell"
    pyscenic aucell \
        "$loom" \
        ${output_dir}/cistarget_out.csv \
        --output ${output_dir}/aucell_out.csv \
        --num_workers "$num_workers" 

    echo "SCENIC analysis is done!"
    return 0
}
