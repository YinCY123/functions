run_bwa_bulk() {
    # Default values
    local ref=""
    local in_dir=""
    local out_dir=""
    local threads=8

    if [[ "$#" -eq 1 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
        echo " Function: run_bwa_bulk"
        echo "Usage:"
        echo "run_bwa_bulk --ref=<reference.fa> --in_dir=<input_folder> --out_dir=<output_folder> [--threads=<N>]"
        echo ""
        echo "Notes:"
        echo "  - Detects single-end vs paired-end by filename conventions (_R1/_R2, _1/_2, -1/-2)."
        echo "  - Recognizes .fastq .fq .fastq.gz .fq.gz."
        echo "  - Produces sorted BAMs in <out_dir>."
        echo "  - Requires bwa and samtools installed."
        return 0
    else
        # Parse key-value arguments
        for arg in "$@"; do
            case $arg in
                --ref=*) ref="${arg#*=}" ;;
                --in_dir=*) in_dir="${arg#*=}" ;;
                --out_dir=*) out_dir="${arg#*=}" ;;
                --threads=*) threads="${arg#*=}" ;;
                *)
                    echo "Unknown parameter: $arg"
                    return 1
                    ;;
            esac
        done
    fi

    # Check required parameters
    if [[ -z "$ref" || -z "$in_dir" || -z "$out_dir" ]]; then
        echo "Usage: run_bwa_bulk --ref=<reference.fa> --in_dir=<input_folder> --out_dir=<output_folder> [--threads=<N>]"
        return 1
    fi

    if [[ ! -f "$ref" ]]; then
        echo "Reference file not found: $ref" >&2
        return 1
    fi

    mkdir -p "$out_dir"

    # Index reference if not already done
    if [[ ! -f "${ref}.bwt" ]]; then
        echo "Indexing reference with bwa..."
        bwa index "$ref" || { echo "bwa index failed" >&2; return 1; }
    fi

    # Enable nullglob so globs with no match expand to nothing
    shopt -s nullglob

    # Recognize fastq extensions
    local -a globs=("$in_dir"/*_R1*.{fastq,fq,fastq.gz,fq.gz} "$in_dir"/*_1*.{fastq,fq,fastq.gz,fq.gz} "$in_dir"/*-1*.{fastq,fq,fastq.gz,fq.gz} "$in_dir"/*.{fastq,fq,fastq.gz,fq.gz})

    # Track processed samples to avoid duplicates
    declare -A seen=()

    for r1 in "${globs[@]}"; do
        [[ -f "$r1" ]] || continue
        r1_basename=$(basename "$r1")

        # Skip explicit R2/_2/-2 files
        if [[ "$r1_basename" == *_R2*.* || "$r1_basename" == *"_2."* || "$r1_basename" == *-2.* ]]; then
            continue
        fi

        # Determine sample name and candidate R2 name
        local sample r2 r2_basename tmp
        if [[ "$r1_basename" == *_R1*.* ]]; then
            sample="${r1_basename%%_R1*}"
            r2_basename="${r1_basename/_R1/_R2}"
            r2="$in_dir/$r2_basename"
        elif [[ "$r1_basename" == *_1.* ]]; then
            sample="${r1_basename%%_1.*}"
            r2_basename="${r1_basename/_1./_2.}"
            r2="$in_dir/$r2_basename"
        elif [[ "$r1_basename" == *-1.* ]]; then
            sample="${r1_basename%%-1.*}"
            r2_basename="${r1_basename/-1./-2.}"
            r2="$in_dir/$r2_basename"
        else
            # single-end: strip .gz then extension
            tmp="${r1_basename%.gz}"
            sample="${tmp%.*}"
            r2=""
        fi

        # Skip if already processed
        if [[ -n "${seen[$sample]}" ]]; then
            continue
        fi
        seen[$sample]=1

        # Use pipefail to catch any pipeline errors
        set -o pipefail

        if [[ -n "$r2" && -f "$r2" ]]; then
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing paired-end sample: $sample"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting BWA alignment..."
            # BWA outputs SAM to stdout, progress to stderr (don't redirect stderr to stdout)
            if ! bwa mem -M -t "$threads" "$ref" "$r1" "$r2" | \
                samtools view -b -@ "$threads" - | \
                samtools sort -@ "$threads" -o "$out_dir/${sample}.sorted.bam" -; then
                echo "ERROR: Alignment failed for sample $sample" >&2
                set +o pipefail
                continue
            fi
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Processing single-end sample: $sample"
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting BWA alignment..."
            # BWA outputs SAM to stdout, progress to stderr (don't redirect stderr to stdout)
            if ! bwa mem -M -t "$threads" "$ref" "$r1" | \
                samtools view -b -@ "$threads" - | \
                samtools sort -@ "$threads" -o "$out_dir/${sample}.sorted.bam" -; then
                echo "ERROR: Alignment failed for sample $sample" >&2
                # Debug: check input file
                if [[ "$r1" == *.gz ]]; then
                    lines=$(zcat "$r1" 2>/dev/null | wc -l)
                else
                    lines=$(wc -l < "$r1" 2>/dev/null)
                fi
                reads=$((lines/4))
                echo "Input file: $r1, Reads: $reads" >&2
                set +o pipefail
                continue
            fi
        fi

        set +o pipefail

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Alignment complete. Indexing BAM..."
        if ! samtools index "$out_dir/${sample}.sorted.bam" 2>&1; then
            echo "ERROR: Indexing failed for sample $sample" >&2
            continue
        fi
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed sample: $sample"
    done

    # restore nullglob
    shopt -u nullglob

    echo "Batch alignment complete. Results in $out_dir"
}
