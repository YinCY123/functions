extract_unmapped_reads() {
    if [ $# -ne 2 ]; then
        echo "Usage: extract_unmapped_reads <input_directory> <output_directory>"
        return 1
    fi

    input_dir=$1
    output_dir=$2

    mkdir -p "$output_dir"

    for bam_file in "$input_dir"/*.bam; do
        if [ -f "$bam_file" ]; then
            base_name=$(basename "$bam_file" .bam)
            echo "processing $base_name..."
            # Extract unmapped reads
            samtools view -f 4 "$bam_file" | awk -v prefix="$output_dir/$base_name" '
            {
                barcode = ""; 

                # Scan optional fields for barcode (CB or CR)
                for (i=12; i<=NF; i++) { 
                    if ($i ~ /^CB:Z:/ || $i ~ /^CR:Z:/) { 
                        split($i, arr, ":")
                        barcode = arr[3]
                        break
                    }
                }

                if (barcode != "") {
                    print "@"$1"\n"barcode"\n+\n"barcode > prefix "_R1.fastq"  # Barcode (R1)
                } else {
                    print "Warning: No barcode found for read " $1 > "/dev/stderr"
                }

                print "@"$1"\n"$10"\n+\n"$11 > prefix "_R2.fastq"  # cDNA (R2)
            }'

            echo "Processed: $bam_file -> $output_dir/${base_name}_R1.fastq and $output_dir/${base_name}_R2.fastq"
        fi
    done

    echo "Batch processing completed!"
}
