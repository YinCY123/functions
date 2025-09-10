trim_bd_wta() {
    # Default values
    threads=8
    verbose=false
    cellbarcode=16
    umi=10
    input_dir=""
    output_dir=""

    # Help message
    show_help() {
        echo "Usage: trim_bd_wta_r2 [options]"
        echo "Options:"
        echo "  -i, --input DIR     Input directory containing sample_R2.fastq.gz files"
        echo "  -o, --output DIR    Output directory for trimmed files (default: input_dir_trimmed)"
        echo "  -t, --threads N     Number of threads to use (default: 8)"
        echo "  -v, --verbose       Show verbose output"
        echo "  -h, --help          Show this help message"
        echo ""
        echo "Example:"
        echo "  trim_bd_wta_r2 -i /path/to/input -o /path/to/output -t 16"
    }

    # Parse options
    while [[ $# -gt 0 ]]; do
        case $1 in
            -i|--input_dir)
                input_dir="$2"
                shift 2
                ;;
            -o|--output_dir)
                output_dir="$2"
                shift 2
                ;;
            -t|--threads)
                threads="$2"
                shift 2
                ;;
            -v|--verbose)
                verbose=true
                shift
                ;;
            -h|--help)
                show_help
                return 0
                ;;
            *)
                echo "Unknown option: $1"
                show_help
                return 1
                ;;
        esac
    done

    # Check if input directory is provided
    if [ -z "$input_dir" ]; then
        echo "Error: Input directory is required"
        show_help
        return 1
    fi

    # Set default output directory if not provided
    if [ -z "$output_dir" ]; then
        output_dir="${input_dir}_trimmed"
    fi

    # Check if input directory exists
    if [ ! -d "$input_dir" ]; then
        echo "Error: Input directory $input_dir does not exist"
        return 1
    fi

    # Create output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Show settings if verbose
    if [ "$verbose" = true ]; then
        echo "Settings:"
        echo "  Input directory: $input_dir"
        echo "  Output directory: $output_dir"
        echo "  Threads: $threads"
        echo "----------------------------------------"
    fi

    # Find all matching files in the input directory
    echo "Processing files in $input_dir..."
    
    # Process each matching file
    for input_file in "$input_dir"/*.fastq.gz; do
        # Skip if no files found
        [ -f "$input_file" ] || continue

        # Get filename without path
        filename=$(basename "$input_file" .fastq.gz)
        
        # Create output filename
        output_file="$output_dir/${filename%.*}_trimmed.fastq.gz"

        if [ "$verbose" = true ]; then
            echo "Processing $filename..."
            echo "Output will be saved to $output_file"
        fi

        # Run cutadapt with specified threads
        cutadapt -j "$threads" -l $((cellbarcode + umi)) -o "$output_file" "$input_file"

        # Check if cutadapt was successful
        if [ $? -eq 0 ]; then
            if [ "$verbose" = true ]; then
                echo "Successfully trimmed $filename"
                echo "Output saved to $output_file"
                echo "----------------------------------------"
            fi
        else
            echo "Error: Trimming failed for $filename"
            echo "----------------------------------------"
        fi
    done
}