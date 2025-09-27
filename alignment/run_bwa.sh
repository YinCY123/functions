#!/usr/bin/env bash

# set -euo pipefail

# run_bwa
#
# Minimal, robust bulk ATAC-seq alignment with BWA (or bwa-mem2 if available).
#
# Pipeline:
#  1) bwa mem (or bwa-mem2 mem) align paired-end FASTQs to reference index
#  2) samtools view: filter unmapped/secondary/supplementary and low MAPQ
#  3) optional mitochondrial removal (default: chrM)
#  4) samtools sort by coordinate
#  5) samtools markdup -r to remove duplicates
#  6) samtools index BAM
#
# Requirements: bwa (or bwa-mem2), samtools >=1.10
#
# Example:
#  source alignment/run_bwa_atac.sh
#  run_bwa \
#    -i /data/idx/hg38/bwa/hg38 \
#    -1 sample_R1.fastq.gz \
#    -2 sample_R2.fastq.gz \
#    -o out/sample \
#    -t 16 \
#    --mapq 30 \
#    --mito chrM
#

run_bwa() {
  local index_prefix=""
  local r1=""
  local r2=""
  local out_prefix=""
  local threads=12
  local mapq=0
  local mito="none"
  local use_mem2=auto
  local require_pp=false
  # Batch-mode options
  local input_dir=""
  local out_dir=""
  local r1_token="_R1"
  local r2_token="_R2"
  local glob_ext=".fastq.gz"

  if [[ $# -eq 0 ]]; then
    _run_bwa_usage
    return 1
  fi

  while [[ $# -gt 0 ]]; do
    case "$1" in
      -i|--index)
        index_prefix="$2"; shift 2 ;;
      -1|--r1)
        r1="$2"; shift 2 ;;
      -2|--r2)
        r2="$2"; shift 2 ;;
      -o|--out-prefix)
        out_prefix="$2"; shift 2 ;;
      -t|--threads)
        threads="$2"; shift 2 ;;
      --mapq)
        mapq="$2"; shift 2 ;;
      --mito)
        mito="$2"; shift 2 ;;
      --require-proper-pair)
        require_pp=true; shift 1 ;;
      --bwa)
        use_mem2="no"; shift 1 ;;
      --bwa-mem2)
        use_mem2="yes"; shift 1 ;;
      -d|--dir)
        input_dir="$2"; shift 2 ;;
      --out-dir)
        out_dir="$2"; shift 2 ;;
      --r1-token)
        r1_token="$2"; shift 2 ;;
      --r2-token)
        r2_token="$2"; shift 2 ;;
      --ext)
        glob_ext="$2"; shift 2 ;;
      -h|--help)
        _run_bwa_usage; return 0 ;;
      *)
        echo "[run_bwa] Unknown argument: $1" >&2
        _run_bwa_usage
        return 1 ;;
    esac
  done

  # Mode detection: batch if input_dir set, else single-sample
  local mode="single"
  if [[ -n "$input_dir" ]]; then
    mode="batch"
  fi

  # Validate inputs
  if [[ -z "$index_prefix" ]]; then
    echo "[run_bwa] Missing index prefix." >&2
    _run_bwa_usage
    return 1
  fi
  if [[ "$mode" == "single" ]]; then
    if [[ -z "$r1" || -z "$r2" || -z "$out_prefix" ]]; then
      echo "[run_bwa] For single-sample: --r1, --r2, and --out-prefix are required." >&2
      _run_bwa_usage
      return 1
    fi
    if [[ ! -f "$r1" ]]; then
      echo "[run_bwa] R1 not found: $r1" >&2; return 1
    fi
    if [[ ! -f "$r2" ]]; then
      echo "[run_bwa] R2 not found: $r2" >&2; return 1
    fi
  else
    if [[ ! -d "$input_dir" ]]; then
      echo "[run_bwa] Input directory not found: $input_dir" >&2; return 1
    fi
    if [[ -z "$out_dir" ]]; then
      echo "[run_bwa] For batch mode: --out-dir is required." >&2; return 1
    fi
    mkdir -p "$out_dir"
  fi

  # Choose aligner
  local aligner="bwa"
  if [[ "$use_mem2" == "yes" ]]; then
    aligner="bwa-mem2"
  elif [[ "$use_mem2" == "auto" ]]; then
    if command -v bwa-mem2 >/dev/null 2>&1; then
      aligner="bwa-mem2"
    elif ! command -v bwa >/dev/null 2>&1; then
      echo "[run_bwa] Neither bwa-mem2 nor bwa found in PATH." >&2; return 1
    fi
  else
    if ! command -v bwa >/dev/null 2>&1; then
      echo "[run_bwa] bwa not found in PATH." >&2; return 1
    fi
  fi
  if ! command -v samtools >/dev/null 2>&1; then
    echo "[run_bwa] samtools not found in PATH." >&2; return 1
  fi

  if [[ "$mode" == "single" ]]; then
    mkdir -p "$(dirname "$out_prefix")"

    local bam_tmp="${out_prefix}.tmp.filtered.sort.bam"
    local bam_dedup="${out_prefix}.bam"
    local bam_index="${bam_dedup}.bai"
    local log_file="${out_prefix}.align.log"

    echo "[run_bwa] Aligner: ${aligner}" | tee "$log_file"
    echo "[run_bwa] Index:   ${index_prefix}" | tee -a "$log_file"
    echo "[run_bwa] R1:      ${r1}" | tee -a "$log_file"
    echo "[run_bwa] R2:      ${r2}" | tee -a "$log_file"
    echo "[run_bwa] Threads: ${threads}" | tee -a "$log_file"
    echo "[run_bwa] MAPQ>=:  ${mapq}" | tee -a "$log_file"
    echo "[run_bwa] Mito:    ${mito}" | tee -a "$log_file"
    echo "[run_bwa] ProperPair required: ${require_pp}" | tee -a "$log_file"

    local sample_id
    sample_id="$(basename "$out_prefix")"
    local rg="@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA"

    # Build filtering arguments dynamically
    local view_args=( -@ "$threads" -h -F 772 )
    if [[ "$require_pp" == true ]]; then
      view_args+=( -f 2 )
    fi
    if [[ "$mapq" -gt 0 ]]; then
      view_args+=( -q "$mapq" )
    fi

    if [[ "$aligner" == "bwa-mem2" ]]; then
      bwa-mem2 mem -t "$threads" -R "$rg" "$index_prefix" "$r1" "$r2" 2>>"$log_file" \
        | samtools view "${view_args[@]}" \
        | awk -v mito="$mito" 'BEGIN{if(mito==""||mito=="none"){skip=0}else{skip=1}} {if($0 ~ /^@/){print; next} if(skip && $3==mito){next} print}' \
        | samtools sort -@ "$threads" -o "$bam_tmp" - 2>>"$log_file"
    else
      bwa mem -t "$threads" -M -R "$rg" "$index_prefix" "$r1" "$r2" 2>>"$log_file" \
        | samtools view "${view_args[@]}" \
        | awk -v mito="$mito" 'BEGIN{if(mito==""||mito=="none"){skip=0}else{skip=1}} {if($0 ~ /^@/){print; next} if(skip && $3==mito){next} print}' \
        | samtools sort -@ "$threads" -o "$bam_tmp" - 2>>"$log_file"
    fi

    samtools markdup -@ "$threads" -r "$bam_tmp" "$bam_dedup" 2>>"$log_file"
    samtools index -@ "$threads" "$bam_dedup" "$bam_index"
    rm -f "$bam_tmp"
    echo "[run_bwa] Done. Output: $bam_dedup" | tee -a "$log_file"
  else
    local r1_list
    shopt -s nullglob
    mapfile -t r1_list < <(find "$input_dir" -maxdepth 1 -type f -name "*${r1_token}*${glob_ext}" | sort)
    shopt -u nullglob

    if [[ ${#r1_list[@]} -eq 0 ]]; then
      echo "[run_bwa] No R1 files found with token '${r1_token}' and ext '${glob_ext}' in: $input_dir" >&2
      return 1
    fi

    echo "[run_bwa] Found ${#r1_list[@]} R1 files. Starting batch..."

    local r1_file r2_file sample base base_noext outp
    for r1_file in "${r1_list[@]}"; do
      r2_file="${r1_file/${r1_token}/${r2_token}}"
      if [[ ! -f "$r2_file" ]]; then
        echo "[run_bwa] Skipping: R2 not found for R1=$r1_file (expected $r2_file)" >&2
        continue
      fi
      base="$(basename "$r1_file")"
      base_noext="$base"
      base_noext="${base_noext%${glob_ext}}"
      sample="${base_noext/${r1_token}/}"
      outp="${out_dir%/}/${sample}"
      run_bwa \
        --index "$index_prefix" \
        --r1 "$r1_file" \
        --r2 "$r2_file" \
        --out-prefix "$outp" \
        --threads "$threads" \
        --mapq "$mapq" \
        --mito "$mito" \
        $( [[ "$require_pp" == true ]] && echo "--require-proper-pair" ) \
        $( [[ "$use_mem2" == "yes" ]] && echo "--bwa-mem2" || ([[ "$use_mem2" == "no" ]] && echo "--bwa") )
    done

    echo "[run_bwa] Batch completed. Outputs in: $out_dir"
  fi
}

_run_bwa_usage() {
  cat >&2 <<EOF
Single sample:
  run_bwa -i INDEX_PREFIX -1 R1.fastq[.gz] -2 R2.fastq[.gz] -o OUT_PREFIX [options]

Batch mode (process all samples in a directory):
  run_bwa -i INDEX_PREFIX -d INPUT_DIR --out-dir OUT_DIR [options]

Required arguments:
  -i, --index         BWA/BWA-MEM2 index prefix (e.g., /ref/hg38/bwa/hg38)
  (single)  -1, --r1            R1 FASTQ(.gz)
  (single)  -2, --r2            R2 FASTQ(.gz)
  (single)  -o, --out-prefix    Output prefix for BAM/logs
  (batch)   -d, --dir           Input directory containing FASTQs
  (batch)       --out-dir       Output directory for per-sample outputs

Options:
  -t, --threads N     Threads (default: 12)
      --mapq N        Minimum MAPQ (default: 0)
      --mito NAME     Mito contig to remove (default: none; e.g., chrM)
      --bwa           Force bwa (disable bwa-mem2)
      --bwa-mem2      Force bwa-mem2
      --require-proper-pair  Require flag 0x2 (properly paired). Off by default
  (batch)   --r1-token STR  Token for R1 filenames (default: _R1)
  (batch)   --r2-token STR  Token for R2 filenames (default: _R2)
  (batch)   --ext STR       File extension to match (default: .fastq.gz)
  -h, --help          Show this help

Outputs:
  OUT_PREFIX.bam      Deduplicated, coordinate-sorted BAM
  OUT_PREFIX.bam.bai  BAM index
  OUT_PREFIX.align.log Log of alignment and filtering steps

Notes:
  - By default, mitochondrial reads are kept (mito=none) and MAPQ is not filtered.
  - Use --mito chrM (or MT) to drop mitochondrial reads.
  - Use --require-proper-pair to enforce proper-pairing; many ATAC reads are not flagged as such.
  - Use --bwa-mem2 if installed for faster alignment.
EOF
}

# If executed directly, run the function with provided args
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  run_bwa "$@"
fi


