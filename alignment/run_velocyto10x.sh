#!/bin/bash

run_velocyto10x(){
    if [[ $1 == "-h"  || $1 == "--help" ]]; then
        echo ""
        echo "--dir                 cellranger output directory"
        echo "-m                    mask file"
        echo "gtf                   gtf annotation file"
        echo "-l                    the logic to use for the filtering (default: Default)"
        echo "-@                    number of threads to use to sort the bam by cellID file using samtools (deafult: 128)"
        echo "--samtools-memory     the number of MB used for every thread by samtools to sort the bam file"
        echo "-t                    the dtype of the loom file layers, if more then 6000 molecular/reads per gene per cell are exprected set uint32 to avoid truncation (deafult: uint32)"
        echo "-d                    debugging purpose only, it will dump a molecular mapping report to hdf5"
        echo ""
        return 0
    fi

    local dir
    local m
    local l=Default
    local samtools_threads=128
    local samtools_memory=50000
    local t=uint32
    local d=0


    while [[ $# -gt 0 ]]; do
        case $1 in
            --dir)
                dir="$2"
                shift 2
                ;;
            -m|--mask)
                m="$2"
                shift 2
                ;;
            --gtf)
                gtf="$2"
                shift 2
                ;;
            -l|--logic)
                l="$2"
                shift 2
                ;;
            --samtools_threads)
                samtools_threads="$2"
                shift 2
                ;;
            --samtools_memory)
                samtools_memory="$2"
                shift 2
                ;;
            -t|--dtype)
                t="$2"
                shift 2
                ;;
            -d|--dump)
                d="$2"
                shift 2
                ;;
            *)
                echo "unknown parameter: $1"
                return 1
                ;;
        esac
    done

    # Check if required --dir parameter is provided
    if [[ -z "$dir" ]]; then
        echo "Error: --dir parameter is required"
        echo "Usage: run_velocyto10x --dir <cellranger_output_directory> [options]"
        return 1
    fi


    # Check if directory exists and is not empty
    if [[ ! -d "$dir" ]]; then
        echo "Error: Directory '$dir' does not exist"
        return 1
    fi
    
    local samples=( "$dir"/* )
    
    # Check if directory contains any files/directories
    if [[ ! -e "${samples[0]}" ]]; then
        echo "Error: Directory '$dir' is empty"
        return 1
    fi

    for sample in "${samples[@]}"; do
        echo "Processing sample: $sample"
        velocyto run10x \
            "$sample" \
            "$m" \
            -l "$l" \
            --samtools-threads "$samtools_threads" \
            --samtools-memory "$samtools_memory" \
            -t "$t" \
            -d "$d" \
            -vvv
        
        echo "Finished processing $sample"
    done
}