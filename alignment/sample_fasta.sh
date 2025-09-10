#!/bin/bash

# Function to sample FASTA files
sample_fasta() {
    # Check if required arguments are provided
    if [ "$#" -lt 4 ]; then
        echo "Usage: sample_fasta <input_fa1> <input_fa2> <output_fa1> <output_fa2> [num_reads|percentage]"
        echo "Arguments:"
        echo "  input_fa1    : Input FASTA file 1"
        echo "  input_fa2    : Input FASTA file 2"
        echo "  output_fa1   : Output FASTA file 1"
        echo "  output_fa2   : Output FASTA file 2"
        echo "  num_reads    : Number of reads to sample (default: 1000)"
        echo "  percentage   : Percentage of reads to sample (e.g., 50 for 50%)"
        echo ""
        echo "Note: If the last argument ends with '%', it will be treated as a percentage"
        return 1
    fi

    # Check if pv is installed
    if ! command -v pv &> /dev/null; then
        echo "Error: 'pv' command is required but not installed."
        echo "Please install it using:"
        echo "  Ubuntu/Debian: sudo apt-get install pv"
        echo "  CentOS/RHEL: sudo yum install pv"
        echo "  macOS: brew install pv"
        return 1
    fi

    # Set variables
    local input_fa1="$1"
    local input_fa2="$2"
    local output_fa1="$3"
    local output_fa2="$4"
    local sample_param="${5:-1000}"  # Default to 1000 reads if not specified

    # Check if input files exist
    if [ ! -f "$input_fa1" ] || [ ! -f "$input_fa2" ]; then
        echo "Error: Input files do not exist"
        return 1
    fi

    # Create temporary directory
    local temp_dir=$(mktemp -d)
    trap 'rm -rf "$temp_dir"' EXIT

    # Extract headers from both files
    echo "Extracting headers..."
    grep "^>" "$input_fa1" | cut -d' ' -f1 > "$temp_dir/headers1.txt"
    grep "^>" "$input_fa2" | cut -d' ' -f1 > "$temp_dir/headers2.txt"

    # Find common headers
    echo "Finding common headers..."
    comm -12 "$temp_dir/headers1.txt" "$temp_dir/headers2.txt" > "$temp_dir/common_headers.txt"

    # Calculate number of reads to sample
    local total_reads=$(wc -l < "$temp_dir/common_headers.txt")
    local num_reads

    # Check if the parameter is a percentage
    if [[ "$sample_param" == *"%" ]]; then
        # Remove the % sign and convert to decimal
        local percentage=${sample_param%\%}
        # Calculate number of reads based on percentage
        num_reads=$(echo "scale=0; $total_reads * $percentage / 100" | bc)
        echo "Sampling $percentage% of reads ($num_reads out of $total_reads)..."
    else
        num_reads=$sample_param
        echo "Sampling $num_reads reads out of $total_reads..."
    fi

    # Randomly sample headers
    shuf -n "$num_reads" "$temp_dir/common_headers.txt" > "$temp_dir/sampled_headers.txt"

    # Create a temporary file with the headers to extract
    echo "Creating extraction patterns..."
    sed 's/^/^/' "$temp_dir/sampled_headers.txt" > "$temp_dir/patterns.txt"

    # Function to extract reads with progress bar
    extract_reads() {
        local input_file="$1"
        local output_file="$2"
        local pattern_file="$3"
        local total_patterns=$(wc -l < "$pattern_file")
        
        echo "Extracting reads from $(basename "$input_file")..."
        # Clear output file if it exists
        > "$output_file"
        
        # Process each pattern and show progress
        while IFS= read -r pattern; do
            awk -v RS=">" -v pattern="$pattern" '
                $0 ~ pattern {
                    print ">" $0
                }
            ' "$input_file" >> "$output_file"
            echo -n "." >&2
        done < "$pattern_file" | pv -l -s "$total_patterns" -N "Progress" > /dev/null
        echo >&2
    }

    # Extract reads from both files with progress bars
    extract_reads "$input_fa1" "$output_fa1" "$temp_dir/patterns.txt"
    extract_reads "$input_fa2" "$output_fa2" "$temp_dir/patterns.txt"

    # Count the number of reads in output files
    local count1=$(grep -c "^>" "$output_fa1")
    local count2=$(grep -c "^>" "$output_fa2")

    echo "Sampling complete!"
    echo "Number of reads in output files:"
    echo "  $output_fa1: $count1 reads"
    echo "  $output_fa2: $count2 reads"

    # Clean up
    rm -rf "$temp_dir"
}

# Example usage:
# sample_fasta input_1.fa input_2.fa output_1.fa output_2.fa 1000    # Sample 1000 reads
# sample_fasta input_1.fa input_2.fa output_1.fa output_2.fa 50%     # Sample 50% of reads