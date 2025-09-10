extract_fastq_barcode() {

    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        echo "Usage: extract_fastq_barcode <input_folder> <output_folder> <file_pattern> <n_bases>"
        echo
        echo "Extracts the first <n_bases> bases from each sequence and quality line in FASTQ files."
        echo
        echo "Arguments:"
        echo "  <input_folder>   Folder containing input FASTQ files"
        echo "  <output_folder>  Folder to save processed FASTQ files"
        echo "  <file_pattern>   Pattern to match FASTQ files (without extension)"
        echo "  <n_bases>        Number of bases to extract from each sequence/quality line"
        return 0
    fi

    local input=$1
    local output=$2
    local file=$3
    local n_bases=$4

    for fastq_file in "$input"/*${file}.fastq.gz; do
        output_file="$output/$(basename "${fastq_file%.fastq.gz}").fastq.gz"
        zcat "$fastq_file" | \
        awk -v n="$n_bases" '{
            if(NR%4==2) $0=substr($0,1,n);
            else if(NR%4==0) $0=substr($0,1,n);
            print
        }' | gzip > "$output_file"
        echo "Processed $fastq_file -> $output_file"
    done
}
