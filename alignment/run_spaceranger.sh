#!/bin/bash

run_spaceranger() {
    # Help message
    if [[ "$1" == "--help" || "$1" == "-h" ]]; then
        echo "Usage: run_spaceranger [options]"
        echo "Mandatory options:"
        echo "  --fastq, -q         Path to directory containing FASTQ files"
        echo "  --transcriptome, -t Path to Spaceranger reference transcriptome"
        echo "  --output, -o        Path to output directory"
        echo "Optional parameters:"
        echo "  --image, -i         Path to tissue image file (optional)"
        echo "  --create-bam, -b    Whether to create BAM file (true/false, default: false)"
        echo "  --cores, -c         Number of cores to use (default: 8)"
        echo "  --mem, -m           Memory in GB to use (default: 32)"
        echo ""
        return 0
    fi

    # Default values
    local fastq=""
    local transcriptome=""
    local output=""
    local image=""
    local slide=""
    local unknown_slide="visium-1"
    local create_bam="false"
    local cores=8
    local mem=32

    # Parse arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --fastq|-q)
                fastq="$2"
                shift 2
                ;;
            --transcriptome|-t)
                transcriptome="$2"
                shift 2
                ;;
            --output|-o)
                output="$2"
                shift 2
                ;;
            --image|-i)
                image="$2"
                shift 2
                ;;
            --unknown_slide)
                unknown_slide="$2"
                shift 2
                ;;
            --create-bam|-b)
                create_bam="$2"
                shift 2
                ;;
            --cores|-c)
                cores="$2"
                shift 2
                ;;
            --mem|-m)
                mem="$2"
                shift 2
                ;;
            *)
                echo "Unknown parameter: $1"
                return 1
                ;;
        esac
    done

    # Check mandatory arguments
    if [[ -z "$fastq" || -z "$transcriptome" || -z "$output" ]]; then
        echo "Error: --fastq, --transcriptome, and --output are required."
        return 1
    fi
    if [[ ! -d "$fastq" ]]; then
        echo "Error: FASTQ directory $fastq not found."
        return 1
    fi
    if [[ ! -d "$transcriptome" ]]; then
        echo "Error: Transcriptome directory $transcriptome not found."
        return 1
    fi
    if [[ ! -d "$output" ]]; then
        mkdir -p "$output"
    fi

    # Find all unique sample names based on CellRanger/SpaceRanger naming convention
    local found_any=0
    local sample_names=()
    while IFS= read -r file; do
        # Extract sample name as everything before _S1_L00
        base=$(basename "$file")
        sample_name=${base%%_S1_L00*}
        # Only add if not already in the array
        if [[ ! " ${sample_names[@]} " =~ " ${sample_name} " ]]; then
            sample_names+=("$sample_name")
        fi
        found_any=1
    done < <(find "$fastq" -maxdepth 1 -type f -name "*_S1_L00*_R1_001.fastq.gz" | sort)

    if [[ $found_any -eq 0 ]]; then
        echo "No FASTQ files matching *_S1_L00*_R1_001.fastq.gz found in $fastq."
        return 1
    fi
    
    # slide
    if [[ "$slide" == "" ]]; then
        local slideParam="--unknown-slide=$unknown_slide"
    fi

    # Loop through each unique sample name
    for sample_name in "${sample_names[@]}"; do
        echo "Processing sample: $sample_name"
        cmd=(spaceranger count \
            --id="$sample_name" \
            --transcriptome="$transcriptome" \
            --fastqs="$fastq" \
            --sample="$sample_name" \
            --localcores="$cores" \
            --localmem="$mem" \
            "$slideParam" \
            --create-bam="$create_bam" \
            --output-dir="$output")
        if [[ -n "$image" ]]; then
            cmd+=(--image="$image")
        fi
        "${cmd[@]}"
        if [[ $? -eq 0 ]]; then
            echo "Spaceranger analysis completed successfully for $sample_name."
        else
            echo "Error: Spaceranger analysis failed for $sample_name."
        fi
    done
}
