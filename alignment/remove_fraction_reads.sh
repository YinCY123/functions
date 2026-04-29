remove_fraction_reads() {
    # Default values
    local bam=""
    local chr=""
    local start=""
    local end=""
    local fraction=""
    local seed=42
    local out="cleaned.bam"
    local outdir="./"

    # -----------------------------
    # Parse key=value parameters
    # -----------------------------
    for arg in "$@"; do
        case $arg in
            bam=*)      bam="${arg#*=}" ;;
            chr=*)      chr="${arg#*=}" ;;
            start=*)    start="${arg#*=}" ;;
            end=*)      end="${arg#*=}" ;;
            fraction=*) fraction="${arg#*=}" ;;
            seed=*)     seed="${arg#*=}" ;;
            out=*)      out="${arg#*=}" ;;
            outdir=*)   outdir="${arg#*=}" ;;
            *)
                echo "Unknown parameter: $arg"
                return 1
                ;;
        esac
    done

    # -----------------------------
    # Validate required parameters
    # -----------------------------
    if [[ -z "$bam" || -z "$chr" || -z "$start" || -z "$end" || -z "$fraction" ]]; then
        echo "ERROR: Required parameters missing."
        echo "Usage:"
        echo "  remove_fraction_reads bam=FILE chr=CHR start=INT end=INT fraction=FLOAT [seed=INT] [outdir=DIR] [out=FILE]"
        return 1
    fi

    if ! awk -v f="$fraction" 'BEGIN { exit !(f ~ /^([0-9]*[.])?[0-9]+$/ && f >= 0 && f <= 1) }'; then
        echo "ERROR: fraction must be a number between 0 and 1."
        return 1
    fi

    mkdir -p "$outdir" || return 1
    if [[ "$out" != */* ]]; then
        out="${outdir%/}/$out"
    fi

    local tmpd
    tmpd=$(mktemp -d "${outdir%/}/.remove_fraction_reads.XXXXXX") || return 1
    trap 'rm -rf "$tmpd"' RETURN

    # -----------------------------
    # Temporary files
    # -----------------------------
    local region_qnames="$tmpd/region.qnames.tmp.txt"
    local drop_qnames="$tmpd/drop.qnames.tmp.txt"
    local filtered_unsorted="$tmpd/filtered.unsorted.tmp.bam"



    # -----------------------------
    # Step 1: Collect unique read names overlapping region
    # -----------------------------
    echo "[1/4] Collecting read names from region..."
    samtools view "$bam" "${chr}:${start}-${end}" | awk '{print $1}' | sort -u > "$region_qnames"

    # -----------------------------
    # Step 2: Randomly choose read names to drop
    # -----------------------------
    echo "[2/4] Sampling read names to drop (fraction=$fraction, seed=$seed)..."
    awk -v frac="$fraction" -v s="$seed" 'BEGIN { srand(s) } { if (rand() >= frac) print $0 }' "$region_qnames" > "$drop_qnames"

    # -----------------------------
    # Step 3: Remove all alignments from dropped read names
    # -----------------------------
    echo "[3/4] Filtering BAM while preserving mate consistency..."
    samtools view -h "$bam" | awk '
        BEGIN {
            while ((getline < "'"$drop_qnames"'") > 0) {
                drop[$1] = 1
            }
            close("'"$drop_qnames"'")
        }
        /^@/ { print; next }
        !($1 in drop) { print }
    ' | samtools view -b -o "$filtered_unsorted" -

    # -----------------------------
    # Step 4: Sort + index
    # -----------------------------
    echo "[4/4] Sorting and indexing..."
    samtools sort -o "$out" "$filtered_unsorted"
    samtools index "$out"

    echo "Done. Final BAM: $out"
}
