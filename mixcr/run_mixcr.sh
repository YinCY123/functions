#!/bin/bash

run_mixcr() {
    # Parse key=value arguments
    for arg in "$@"; do
        case $arg in
            preset=*) preset="${arg#*=}" ;;
            tag_pattern=*) tag_pattern="${arg#*=}" ;;
            input_dir=*) input_dir="${arg#*=}" ;;
            output_dir=*) output_dir="${arg#*=}" ;;
            species=*) species="${arg#*=}" ;;
        esac
    done

    # Show help if requested
    if [[ "$1" == "-h" || "$1" == "--help" || "$show_help" == "true" ]]; then
        cat <<EOF
Usage: run_mixcr input_dir=<path> output_dir=<path> [species=<code>] [help=true]

Parameters:
  preset      MiXCR preset to use (e.g., 10x-sc-xcr-vdj, rna-seq, takara-mouse-rna-tcr-smarter, takara-mouse-rna-tcr-umi-smarseq, or takara-sc-human-rna-tcr-smarter), To list all built-in presets run "mixcr listPresets".
  tag_pattern   UMI and Cell barcode pattern for R1
  input_dir   Directory containing FASTQ files (named sample_R1.fastq.gz and sample_R2.fastq.gz)
  output_dir  Directory to store MiXCR outputs
  species     Species code (default: hsa for human)
  help        Show this help message

Examples:
  run_mixcr input_dir=./data output_dir=./results/tcr species=hsa
  run_mixcr input_dir=./data output_dir=./results/bcr

Notes:
  - Replace <10x-scrna-tcr-preset> and <10x-scrna-ig-preset> with the correct MiXCR presets for your 10x chemistry.
  - The function loops over all *_R1.fastq.gz files in input_dir and expects matching *_R2.fastq.gz files.
EOF
        return 0
    fi

    # Defaults
    species=${species:-hsa}
    preset=${preset:-"generic-single-cell-gex-with-umi"}
    # generic-single-cell-gex-with-umi, generic-single-cell-gex, 10x-sc-5gex, split-seq-3gex, rna-seq, bd-sc-xcr-rhapsody-full-length
    
    # Tag pattern for Drop-seq: cell barcode + UMI, then transcript sequence
    # Format: ^(CELL:N{X})(UMI:N{Y})(R1:*)
    # MiXCR requires molecule tag names to start with "umi" or "mi" (case-insensitive)
    # Default: 16bp cell barcode + 10bp UMI (adjust based on your Drop-seq protocol)
    # tag_pattern=${tag_pattern:-'(CELL1:N{9})N{12}(CELL2:N{9})N{13}(CELL3:N{9})(UMI:N{8})\^(R2:*)'}
    tag_pattern=${tag_pattern:-"(CELL:N{16})(UMI:N{12})\^(R2:*)"} #  GSE290927 (3' v3), GSE286827 (3' v3.1)
    # tag_pattern=${tag_pattern:-"(CELL:N{16})(UMI:N{10})\^(R2:*)"} # GSE273138 (5' v1.1, 16C10U, R2 barcode),  GSE252490(3' v2), GSE206332 (3' v2)


    if [[ -z "$input_dir" || -z "$output_dir" ]]; then
        echo "Error: input_dir and output_dir must be specified"
        return 1
    fi

    mkdir -p "$output_dir"

    # avoid literal globbing when no files match
    shopt -s nullglob
    for R1 in "$input_dir"/*_R1.fastq.gz; do
        SAMPLE=$(basename "$R1" _R1.fastq.gz)
        R2="$input_dir/${SAMPLE}_R2.fastq.gz"

        # Skip if paired R2 is missing
        if [[ ! -f "$R2" ]]; then
            echo "Warning: missing R2 for sample $SAMPLE, skipping."
            continue
        fi

        # Run MiXCR analysis 
        mixcr -Xmx100g analyze "$preset" \
                --tag-pattern "$tag_pattern" \
                --species "$species" \
                "$R1" "$R2" \
                "$output_dir/$SAMPLE" 
                # --tag-pattern "$tag_pattern" 

        # Determine which .clns file was produced (try .clns then .contigs.clns)
        if [[ -f "$output_dir/$SAMPLE.clns" ]]; then
            CLNS_INPUT="$output_dir/$SAMPLE.clns"
        elif [[ -f "$output_dir/${SAMPLE}.contigs.clns" ]]; then
            CLNS_INPUT="$output_dir/${SAMPLE}.contigs.clns"
        else
            echo "Error: no .clns file found for $SAMPLE in $output_dir, skipping export."
            continue
        fi

        # Export clones with per-cell barcodes using the detected CLNS file
        OUTPUT_FILE="$output_dir/$SAMPLE.clones.tsv"
        
        # Remove output file if it exists to avoid conflicts
        [[ -f "$OUTPUT_FILE" ]] && rm -f "$OUTPUT_FILE"
        
        # Run exportClones command
        # Note: MiXCR uses -fieldName flags, not --fields
        # Fields are chosen to be compatible with loadContigs from scRepertoire
        # Required columns for scRepertoire::loadContigs(format="MiXCR"):
        #   - tagValueCELL (from -tags Cell, may need renaming if exported as tags.Cell)
        #   - topChains (from -topChains)
        #   - readCount (from -readCount)
        # Note: Flag names are case-sensitive! Use -cloneId (lowercase 'd'), not -cloneID
        mixcr -Xmx100g exportClones \
            -cloneId \
            -readCount \
            -topChains \
            -tags Cell \
            "$CLNS_INPUT" \
            "$OUTPUT_FILE" 2>&1
        
        # Check if output file was successfully created
        if [[ ! -f "$OUTPUT_FILE" ]]; then
            echo "Error: exportClones failed for $SAMPLE - output file not created: $OUTPUT_FILE"
            continue
        elif [[ ! -s "$OUTPUT_FILE" ]]; then
            echo "Warning: exportClones produced empty file for $SAMPLE: $OUTPUT_FILE"
            continue
        else
            echo "Successfully exported clones for $SAMPLE to $OUTPUT_FILE"
        fi
    done
    shopt -u nullglob
}