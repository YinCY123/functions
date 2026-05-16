reverse_fraction_region() {
    # --- Usage help ---
    if [[ "$1" == "help" || "$1" == "--help" || "$1" == "-h" ]]; then
        cat <<EOF
reverse_fraction_region — reverse a FRACTION of reads in a genomic region (NO complement)

USAGE:
  reverse_fraction_region \\
      in=input.bam \\
      region=chr:start-end \\
      frac=0.3 \\
      out=final.bam

REQUIRED PARAMETERS:
  in=FILE            Input BAM file (indexed recommended)
  region=STRING      Genomic region (e.g. chr1:100000-200000), 1-based inclusive end
  frac=FLOAT         Fraction of reads to reverse (0–1)
  out=FILE           Output BAM file

DESCRIPTION:
  - Splits reads overlapping vs outside the region (BED + samtools -L/-U)
  - Randomly subsamples a fraction (frac) of templates in the overlap set
  - Reverses SEQ and QUAL for selected reads (bioawk reverse(); no complement)
  - Writes untouched overlap reads and merges with outside reads
  - Final output is coordinate-sorted BAM + index

EXAMPLE:
  reverse_fraction_region in=a.bam region=chr1:1-100000 frac=0.25 out=out.bam

EOF
        return 0
    fi

    # --- Parse key=value arguments ---
    for arg in "$@"; do
        key="${arg%%=*}"
        val="${arg#*=}"
        case "$key" in
            in)      IN="$val" ;;
            region)  REGION="$val" ;;
            frac)    FRAC="$val" ;;
            out)     OUT="$val" ;;
            *) echo "Unknown parameter: $key" >&2; return 1 ;;
        esac
    done

    # --- Check required parameters ---
    if [[ -z "$IN" || -z "$REGION" || -z "$FRAC" || -z "$OUT" ]]; then
        echo "Missing required parameters." >&2
        echo "Run: reverse_fraction_region help" >&2
        return 1
    fi

    # --- region=chr:start-end → BED (0-based start, half-open end) ---
    if [[ ! "$REGION" =~ ^([^:]+):([0-9]+)-([0-9]+)$ ]]; then
        echo "region must look like chr:start-end (1-based inclusive coordinates)." >&2
        return 1
    fi
    R_CHR="${BASH_REMATCH[1]}"
    R_START="${BASH_REMATCH[2]}"
    R_END="${BASH_REMATCH[3]}"

    # --- Temp workspace ---
    TMP=$(mktemp -d)
    trap 'rm -rf "$TMP"' EXIT

    BED_START=$((R_START - 1))
    BED_END="$R_END"
    printf '%s\t%d\t%d\n' "$R_CHR" "$BED_START" "$BED_END" > "$TMP/region.bed"

    echo "[1] Split overlapping vs outside region (BED + -L/-U)..."
    samtools view -b \
        -L "$TMP/region.bed" \
        -U "$TMP/outside_region.bam" \
        -o "$TMP/region.bam" \
        "$IN"

    echo "[2] Subsample fraction $FRAC (by template)..."
    samtools view -b \
        --subsample "$FRAC" \
        --subsample-seed 0 \
        -o "$TMP/region_to_reverse.bam" \
        "$TMP/region.bam"

    echo "[3] Reverse SEQ + QUAL (NO complement); reattach header (bioawk -c sam drops @SQ)..."
    samtools view -H "$TMP/region_to_reverse.bam" > "$TMP/rev_hdr.sam"
    samtools view "$TMP/region_to_reverse.bam" \
        | bioawk -c sam '
        {
            $10 = reverse($10)
            $11 = reverse($11)
            print
        }
    ' > "$TMP/rev_body.sam"
    cat "$TMP/rev_hdr.sam" "$TMP/rev_body.sam" \
        | samtools view -b -o "$TMP/region_reversed.bam" -

    echo "[4] Overlap reads not chosen for reversal (negated qname list)..."
    samtools view "$TMP/region_to_reverse.bam" | cut -f1 | sort -u > "$TMP/rev_qnames.txt"
    samtools view -h "$TMP/region.bam" \
        | samtools view -b -N "^$TMP/rev_qnames.txt" -o "$TMP/region_untouched.bam" -

    echo "[5] Merge all parts..."
    samtools merge -f "$TMP/merged.bam" \
        "$TMP/outside_region.bam" \
        "$TMP/region_untouched.bam" \
        "$TMP/region_reversed.bam"

    echo "[6] Sort + index..."
    samtools sort -N -o "$OUT" "$TMP/merged.bam"
    samtools index "$OUT"

    echo "Done. Output: $OUT"
}
