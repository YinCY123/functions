run_velocyto() {
    # Usage: run_velocyto -b <BAMFILE1|BAM_DIR> [-b <BAMFILE2|BAM_DIR> ...] -g <GTFFILE> -o <OUTPUT_DIR> [-n <NCORES>] [-s <SAMPLENAME1> ...]
    local BAMFILES=() GTFFILE OUTPUTDIR NCORES SAMPLENAMES=()

    while [[ $# -gt 0 ]]; do
        case $1 in
            -b|--bam)
                if [[ -d "$2" ]]; then
                    # Add all .bam files in the directory
                    for f in "$2"/*.bam; do
                        [[ -e "$f" ]] && BAMFILES+=("$f")
                    done
                else
                    BAMFILES+=("$2")
                fi
                shift 2
                ;;
            -g|--gtf)
                GTFFILE="$2"
                shift 2
                ;;
            -o|--outputDir)
                OUTPUTDIR="$2"
                shift 2
                ;;
            -n|--ncores)
                NCORES="$2"
                shift 2
                ;;
            -s|--sampleName)
                SAMPLENAMES+=("$2")
                shift 2
                ;;
            -h|--help)
                echo "Usage: run_velocyto -b <BAMFILE1|BAM_DIR> [-b <BAMFILE2|BAM_DIR> ...] -g <GTFFILE> -o <OUTPUT_DIR> [-n <NCORES>] [-s <SAMPLENAME1> ...]"
                return 0
                ;;
            *)
                echo "Unknown option: $1"
                return 1
                ;;
        esac
    done

    NCORES=${NCORES:-12}

    if [[ ${#BAMFILES[@]} -eq 0 || -z "$GTFFILE" || -z "$OUTPUTDIR" ]]; then
        echo "Error: At least one BAMFILE, GTFFILE, and OUTPUTDIR are required."
        echo "Use -h for help."
        return 1
    fi

    mkdir -p "$OUTPUTDIR"

    for i in "${!BAMFILES[@]}"; do
        BAM="${BAMFILES[$i]}"
        SAMPLE="${SAMPLENAMES[$i]:-sample_$((i+1))}"
        velocyto run \
            -o "$OUTPUTDIR" \
            -e "$SAMPLE" \
            -@ "$NCORES" \
            "$BAM" \
            "$GTFFILE"
    done
}