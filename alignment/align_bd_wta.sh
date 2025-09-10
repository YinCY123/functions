#!/bin/bash

# Function to align BD WTA files (single or multiple)
align_bd_wta() {
    # Default values
    local input=""
    local output=""
    local reference_genome=""
    local threads=4
    local max_parallel_jobs=4
    local verbose=false
    local dry_run=false

    # Function to display help information
    show_help() {
        echo "BD WTA Alignment Tool"
        echo "====================="
        echo
        echo "This tool aligns BD WTA (Whole Transcriptome Analysis) FASTQ files using STAR aligner."
        echo
        echo "Usage:"
        echo "  align_bd_wta [OPTIONS]"
        echo
        echo "Options:"
        echo "  -i, --input <path>           Input FASTQ file or directory containing FASTQ files"
        echo "  -o, --output <path>          Output BAM file or directory for BAM files"
        echo "  -r, --reference <path>       Path to STAR reference genome directory"
        echo "  -t, --threads <number>       Number of threads per sample (default: 4)"
        echo "  -j, --jobs <number>          Maximum number of parallel jobs (default: 4)"
        echo "  -v, --verbose               Enable verbose output"
        echo "  -d, --dry-run               Show what would be done without actually doing it"
        echo "  -h, --help                  Show this help message"
        echo
        echo "Examples:"
        echo "  1. Process a single FASTQ file:"
        echo "     align_bd_wta -i sample.fastq.gz -o output.bam -r /path/to/reference -t 8"
        echo
        echo "  2. Process a directory of FASTQ files:"
        echo "     align_bd_wta -i /path/to/fastq_files -o /path/to/output -r /path/to/reference -t 8 -j 4"
        echo
        echo "Input File Requirements:"
        echo "  - Supported file extensions: .fastq, .fq, .fastq.gz, .fq.gz"
        echo "  - For paired-end reads, files should be named with _R1 and _R2 suffixes"
        echo
        echo "Output:"
        echo "  - Creates sorted BAM files with the naming pattern: <sample_name>Aligned.sortedByCoord.out.bam"
        echo "  - Generates BAM index files (.bai)"
        echo
        echo "Requirements:"
        echo "  - STAR aligner must be installed and available in PATH"
        echo "  - samtools must be installed and available in PATH"
        echo "  - Sufficient disk space for BAM files and temporary files"
        echo
        echo "For more information, please refer to the STAR documentation:"
        echo "https://github.com/alexdobin/STAR"
    }

    # Parse command line options
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input)
                input="$2"
                shift 2
                ;;
            -o|--output)
                output="$2"
                shift 2
                ;;
            -r|--reference)
                reference_genome="$2"
                shift 2
                ;;
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -j|--jobs)
                max_parallel_jobs="$2"
                shift 2
                ;;
            -v|--verbose)
                verbose=true
                shift
                ;;
            -d|--dry-run)
                dry_run=true
                shift
                ;;
            -h|--help)
                show_help
                return 0
                ;;
            *)
                echo "Error: Unknown option: $1"
                echo "Use 'align_bd_wta --help' for usage information"
                return 1
                ;;
        esac
    done

    # Validate required parameters
    if [ -z "$input" ] || [ -z "$output" ] || [ -z "$reference_genome" ]; then
        echo "Error: Missing required parameters"
        echo "Required parameters: -i/--input, -o/--output, -r/--reference"
        echo "Use 'align_bd_wta --help' for usage information"
        return 1
    fi

    # Validate numeric parameters
    if ! [[ "$threads" =~ ^[0-9]+$ ]] || [ "$threads" -lt 1 ]; then
        echo "Error: Threads must be a positive integer"
        return 1
    fi

    if ! [[ "$max_parallel_jobs" =~ ^[0-9]+$ ]] || [ "$max_parallel_jobs" -lt 1 ]; then
        echo "Error: Maximum parallel jobs must be a positive integer"
        return 1
    fi

    # Validate reference genome
    if [ ! -f "$reference_genome" ]; then
        echo "Error: Reference genome file not found: $reference_genome"
        return 1
    fi

    # Function to process a single file
    process_single_file() {
        local input_fastq="$1"
        local output_bam="$2"
        local ref_genome="$3"
        local num_threads="$4"

        if [ "$dry_run" = true ]; then
            echo "[DRY RUN] Would process: $input_fastq -> $output_bam"
            return 0
        fi

        # Create output directory if it doesn't exist
        local output_dir=$(dirname "$output_bam")
        mkdir -p "$output_dir"

        if [ "$verbose" = true ]; then
            echo "Starting BD WTA alignment for $(basename "$input_fastq")..."
            echo "Input: $input_fastq"
            echo "Output: $output_bam"
            echo "Reference: $ref_genome"
            echo "Threads: $num_threads"
        fi

        # Run STAR alignment
        STAR \
            --runMode alignReads \
            --runThreadN "$num_threads" \
            --genomeDir "$ref_genome" \
            --readFilesIn "$input_fastq" \
            --outFileNamePrefix "${output_bam%.bam}" \
            --outSAMtype BAM SortedByCoordinate \
            --outBAMsortingThreadN "$num_threads" \
            --outSAMunmapped Within \
            --outSAMattributes Standard \
            --twopassMode Basic \
            --chimOutType WithinBAM \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --alignIntronMax 100000 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.04 \
            --outFilterType BySJout \
            --outFilterScoreMinOverLread 0.33 \
            --outFilterMatchNminOverLread 0.33 \
            --limitSjdbInsertNmax 1000000 \
            --limitBAMsortRAM 10000000000

        # Check if alignment was successful
        if [ $? -eq 0 ]; then
            if [ "$verbose" = true ]; then
                echo "Alignment completed successfully for $(basename "$input_fastq")"
            fi
            # Index the BAM file
            samtools index "${output_bam%.bam}Aligned.sortedByCoord.out.bam"
            if [ $? -eq 0 ]; then
                if [ "$verbose" = true ]; then
                    echo "BAM file indexed successfully for $(basename "$input_fastq")"
                fi
                return 0
            else
                echo "Error: Failed to index BAM file for $(basename "$input_fastq")"
                return 1
            fi
        else
            echo "Error: Alignment failed for $(basename "$input_fastq")"
            return 1
        fi
    }

    # Check if input is a directory or a file
    if [ -d "$input" ]; then
        # Process multiple files
        if [ "$verbose" = true ]; then
            echo "Processing directory of FASTQ files..."
        fi
        
        # Create output directory
        mkdir -p "$output"

        # Find all FASTQ files in the input directory
        local fastq_files=($(find "$input" -name "*.fastq" -o -name "*.fq" -o -name "*.fastq.gz" -o -name "*.fq.gz"))
        
        if [ ${#fastq_files[@]} -eq 0 ]; then
            echo "Error: No FASTQ files found in $input"
            return 1
        fi

        if [ "$verbose" = true ]; then
            echo "Found ${#fastq_files[@]} FASTQ files to process"
            echo "Processing with $threads threads per sample"
            echo "Maximum parallel jobs: $max_parallel_jobs"
        fi

        # Process files in parallel
        local i=0
        for fastq in "${fastq_files[@]}"; do
            # Get base name without extension
            local base_name=$(basename "$fastq" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
            local output_bam="$output/${base_name}.bam"
            
            # Run alignment in background
            process_single_file "$fastq" "$output_bam" "$reference_genome" "$threads" &
            
            # Increment counter
            ((i++))
            
            # If we've reached max parallel jobs, wait for one to finish
            if [ $i -eq $max_parallel_jobs ]; then
                wait -n
                ((i--))
            fi
        done

        # Wait for remaining jobs to finish
        wait

        if [ "$verbose" = true ]; then
            echo "All alignments completed"
        fi
        
        # Check if all BAM files were created successfully
        local success=true
        for fastq in "${fastq_files[@]}"; do
            local base_name=$(basename "$fastq" | sed -E 's/\.(fastq|fq)(\.gz)?$//')
            local output_bam="$output/${base_name}Aligned.sortedByCoord.out.bam"
            if [ ! -f "$output_bam" ]; then
                echo "Error: Failed to create BAM file for $base_name"
                success=false
            fi
        done

        if [ "$success" = true ]; then
            if [ "$verbose" = true ]; then
                echo "All samples processed successfully"
            fi
            return 0
        else
            echo "Some samples failed to process"
            return 1
        fi

    elif [ -f "$input" ]; then
        # Process single file
        process_single_file "$input" "$output" "$reference_genome" "$threads"
        return $?
    else
        echo "Error: Input path does not exist: $input"
        return 1
    fi
}

# Export the function so it can be used in other scripts
export -f align_bd_wta
