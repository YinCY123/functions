# Function to run Cell Ranger for sequence alignment with samples derived from FASTQ directory
run_cellranger() {
    local fastq_dir=""
    local reference=""
    local output_dir=""
    local chemistry=auto
    local r1_length=0
    local r2_length=0
    local include_introns=true
    local jobmode=local
    local create_bam=true
    local localcores=10
    local localvmem=85
    local localmem=85
    local maxjobs=1

    # Check for help option before processing other parameters
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        echo ""
        echo " Mendatory parameters: "
        echo "  -q, --fastq_dir <dir>           Directory containing FASTQ files."
        echo "  -g, --reference <path>          Path to the reference genome."
        echo "  -o, --output_dir <dir>          Directory to store output results."
        echo ""
        echo "OPtional Parameters: "
        echo "  -c, --chemistry <type>          Chemistry type (optional)."
        echo "  -b, --create_bam <true|false>   Whether to create BAM files (default: true)."
        echo "  --r1-length                     Hard trim the input Read 1 to this length before analysis"
        echo "  --r2-length                     Hard trim the input Read 1 to this length before analysis"
        echo "  --include-introns               posible value true|false, (default: true)"
        echo "  --jobmode <chr>"
        echo "  --localvmem <num>               virtual memory allowed to use (default: 65GB)"
        echo "  --localmem <num>                memory allowed to use (deafult: 65GB)"
        echo "  -n, --localcores <num>          Number of local cores to use (default: 10)."
        echo "  --maxjobs <num>                 number of jobs (default: 1)"
        echo "-h, --help                        Display this help message."
        echo ""
        return 0
    fi

    while [[ $# -gt 0 ]]; do
        case $1 in
            --fastq_dir|-q)
                fastq_dir="$2"
                shift 2
                ;;
            --reference|-g)
                reference="$2"
                shift 2
                ;;
            --output_dir|-o)
                output_dir="$2"
                shift 2
                ;;
            --chemistry|-c)
                chemistry="$2"
                shift 2
                ;;
            --create_bam|-b)
                create_bam="$2"
                shift 2
                ;;
            --localcores|-n)
                localcores="$2"
                shift 2
                ;;
            --localvmem)
                localvmem="$2"
                shift 2
                ;;
            --localmem)
                localmem="$2"
                shift 2
                ;;
            --r1-length)
                r1_length="$2"
                shift 2
                ;;
            --r2-length)
                r2_length="$2"
                shift 2
                ;;
            --include-introns)
                include_introns="$2"
                shift 2
                ;;
            --jobmode)
                jobmode="$2"
                shift 2
                ;;
            --maxjobs)
                maxjobs="$2"
                shift 2
                ;;
            *)
                echo "unknown parameter: $1"
                return 1
                ;;
        esac
    done

    if [[ -z "$fastq_dir" || -z "$reference" || -z "$output_dir" ]]; then
        echo "Error: Missing required parameters."
        return 1
    fi

    # Get unique sample names from FASTQ files
    local sample_names=($(ls "$fastq_dir"/*.fastq.gz | sed 's/.*\/\([^/]*\)_S[0-9]*_L[0-9]*_R[12].*\.fastq\.gz/\1/' | sort -u))

    for sample_name in "${sample_names[@]}"; do
        echo "Processing sample: $sample_name"
        cellranger count --id="$sample_name" \
                         --transcriptome="$reference" \
                         --fastqs="$fastq_dir" \
                         --sample="$sample_name" \
                         --output-dir="$output_dir/$sample_name" \
                         --create-bam="$create_bam" \
                         --localcores="$localcores" \
                         --chemistry="$chemistry" \
                         --localvmem="$localvmem" \
                         --localmem="$localmem" \
                         --jobmode="$jobmode" \
                        #  --r1-length="$r1_length" \
                        #  --r2-length="$r2_length" \
                         --include-introns="$include_introns" 
    done
}

# Example usage
# run_cellranger -q "/path/to/fastq" -g "/path/to/reference" -o "/path/to/output"