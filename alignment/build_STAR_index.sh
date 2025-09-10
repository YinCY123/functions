#!/bin/bash

build_STAR_index() {
    # Help message
    if [[ "$1" == "--help" || "$1" == "-h" ]]; then
        echo "Usage: build_STAR_index [options]"
        echo "Mandatory options:"
        echo "  --genomeDir, -g     Path to output STAR genome index directory"
        echo "  --fasta, -f         Path to genome FASTA file"
        echo "  --gtf, -a           Path to annotation GTF file"
        echo "Optional parameters:"
        echo "  --sjdbOverhang, -s  Length of the genomic sequence around annotated junctions (default: 100)"
        echo "  --nthreads, -n      Number of threads to use (default: 8)"
        echo ""
        return 0
    fi

    # Default values
    local genomeDir=""
    local fasta=""
    local gtf=""
    local sjdbOverhang=100
    local nthreads=8

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --genomeDir|-g)
                genomeDir="$2"
                shift 2
                ;;
            --fasta|-f)
                fasta="$2"
                shift 2
                ;;
            --gtf|-a)
                gtf="$2"
                shift 2
                ;;
            --sjdbOverhang|-s)
                sjdbOverhang="$2"
                shift 2
                ;;
            --nthreads|-n)
                nthreads="$2"
                shift 2
                ;;
            *)
                echo "Unknown parameter: $1"
                return 1
                ;;
        esac
    done

    # Check mandatory arguments
    if [[ -z "$genomeDir" || -z "$fasta" || -z "$gtf" ]]; then
        echo "Error: --genomeDir, --fasta, and --gtf are required."
        return 1
    fi

    # Check input files
    if [[ ! -f "$fasta" ]]; then
        echo "Error: FASTA file $fasta not found."
        return 1
    fi
    if [[ ! -f "$gtf" ]]; then
        echo "Error: GTF file $gtf not found."
        return 1
    fi

    # Create genomeDir if it doesn't exist
    if [[ ! -d "$genomeDir" ]]; then
        mkdir -p "$genomeDir"
    fi

    # Run STAR genomeGenerate
    STAR \
        --runMode genomeGenerate \
        --genomeDir "$genomeDir" \
        --genomeFastaFiles "$fasta" \
        --sjdbGTFfile "$gtf" \
        --sjdbOverhang "$sjdbOverhang" \
        --runThreadN "$nthreads"

    if [[ $? -eq 0 ]]; then
        echo "STAR genome index built successfully in $genomeDir"
    else
        echo "Error: STAR genome index build failed."
        return 1
    fi
} 