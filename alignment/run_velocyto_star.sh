# Batch velocyto for BAMs directly under one folder, key=value args
#
# Expected layout:
#   root/
#     sample1.bam
#     sample2.bam
#   barcodes_dir/
#     sample1.tsv
#     sample2.tsv
#
# Example:
#   velocyto_bam_folder_kv \
#     root=/data/bams \
#     barcodes_dir=/data/barcodes \
#     gtf=/ref/gencode.vM31.annotation.gtf \
#     repeat_gtf=/ref/rmsk.gtf \
#     outroot=/data/velocyto_out \
#     threads=8 \
#     bam_glob="*.bam" \
#     barcode_suffix=".tsv"
#
# Required keys: root, barcodes_dir, gtf, repeat_gtf
# Optional keys:
#   outroot=velocyto_out
#   threads=8
#   bam_glob=*.bam
#   barcode_suffix=.tsv

velocyto_bam_folder_kv () {
  set -euo pipefail

  # defaults
  local root=""
  local barcodes_dir=""
  local gtf=""
  local repeat_gtf=""
  local outroot="velocyto_out"
  local threads="8"
  local bam_glob="*.bam"
  local barcode_suffix=".tsv"

  # parse key=value args
  local kv key val
  for kv in "$@"; do
    [[ "$kv" == *=* ]] || { echo "Invalid argument: $kv (expected key=value)" >&2; return 1; }
    key="${kv%%=*}"
    val="${kv#*=}"

    case "$key" in
      root)           root="$val" ;;
      barcodes_dir)   barcodes_dir="$val" ;;
      gtf)            gtf="$val" ;;
      repeat_gtf)     repeat_gtf="$val" ;;
      outroot)        outroot="$val" ;;
      threads)        threads="$val" ;;
      bam_glob)       bam_glob="$val" ;;
      barcode_suffix) barcode_suffix="$val" ;;
      *)
        echo "Unknown key: $key" >&2
        return 1
        ;;
    esac
  done

  # required args
  [[ -n "$root" ]]        || { echo "Missing required key: root" >&2; return 1; }
  [[ -n "$barcodes_dir" ]]|| { echo "Missing required key: barcodes_dir" >&2; return 1; }
  [[ -n "$gtf" ]]         || { echo "Missing required key: gtf" >&2; return 1; }
  [[ -n "$repeat_gtf" ]]  || { echo "Missing required key: repeat_gtf" >&2; return 1; }

  mkdir -p "$outroot"
  shopt -s nullglob

  local bam sample sample_out barcodes
  for bam in "$root"/$bam_glob; do
    [[ -f "$bam" ]] || continue

    sample="$(basename "${bam%.*}")"     # sample1 from sample1.bam
    sample_out="$outroot/$sample"
    barcodes="$barcodes_dir/${sample}${barcode_suffix}"

    if [[ ! -f "$barcodes" ]]; then
      echo "Skipping $sample: missing barcode file $barcodes"
      continue
    fi

    mkdir -p "$sample_out"

    if [[ ! -f "${bam}.bai" ]]; then
      samtools index -@ "$threads" "$bam"
    fi

    echo "Running velocyto for $sample..."
    velocyto run \
      -b "$barcodes" \
      -o "$sample_out" \
      -m "$repeat_gtf" \
      "$bam" \
      "$gtf"
  done
}