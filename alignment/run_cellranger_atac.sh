#!/usr/bin/env bash

# run_cellranger_atac
#
# Align single cell ATAC-seq data with Cell Ranger ATAC.
#
# Pipeline:
#  1) Validate input FASTQs and reference
#  2) Create or use existing Cell Ranger ATAC reference
#  3) Run cellranger-atac count with configurable parameters
#  4) Generate summary reports
#
# Requirements: Cell Ranger ATAC, valid reference genome
#
# Example:
#  source alignment/run_cellranger_atac.sh
#  run_cellranger_atac \
#    --fastqs /path/to/fastqs \
#    --sample Sample1 \
#    --reference /ref/refdata-cellranger-atac-GRCh38-1.2.0 \
#    --output-dir /path/to/output

run_cellranger_atac() {
  local reference=""
  local output_dir=""
  local localcores=16
  local localmem=64
  local expect_cells=""
  local chemistry="auto"
  local nosecondary=false
  local dry_run=false
  # Batch mode options
  local input_dir=""
  local sample_sheet=""
  local use_sample_sheet=false

  if [[ $# -eq 0 ]]; then
    _run_cellranger_atac_usage
    return 1
  fi

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --input-dir)
        input_dir="$2"; shift 2 ;;
      --output-dir)
        output_dir="$2"; shift 2 ;;
      --reference)
        reference="$2"; shift 2 ;;
      --localcores)
        localcores="$2"; shift 2 ;;
      --localmem)
        localmem="$2"; shift 2 ;;
      --expect-cells)
        expect_cells="$2"; shift 2 ;;
      --chemistry)
        chemistry="$2"; shift 2 ;;
      --no-secondary)
        nosecondary=true; shift 1 ;;
      --sample-sheet)
        sample_sheet="$2"; use_sample_sheet=true; shift 2 ;;
      --dry-run)
        dry_run=true; shift 1 ;;
      -h|--help)
        _run_cellranger_atac_usage; return 0 ;;
      *)
        echo "[run_cellranger_atac] Unknown argument: $1" >&2
        _run_cellranger_atac_usage
        return 1 ;;
    esac
  done

  # Validate required arguments
  if [[ -z "$input_dir" || -z "$reference" || -z "$output_dir" ]]; then
    echo "[run_cellranger_atac] Required arguments: --input-dir, --reference, and --output-dir." >&2
    _run_cellranger_atac_usage
    return 1
  fi
  if [[ "$use_sample_sheet" == true && -z "$sample_sheet" ]]; then
    echo "[run_cellranger_atac] Sample sheet specified but not provided." >&2
    _run_cellranger_atac_usage
    return 1
  fi

  # Validate inputs
  if [[ ! -d "$input_dir" ]]; then
    echo "[run_cellranger_atac] Input directory not found: $input_dir" >&2
    return 1
  fi
  if [[ "$use_sample_sheet" == true && ! -f "$sample_sheet" ]]; then
    echo "[run_cellranger_atac] Sample sheet not found: $sample_sheet" >&2
    return 1
  fi

  if [[ ! -d "$reference" ]]; then
    echo "[run_cellranger_atac] Reference directory not found: $reference" >&2
    return 1
  fi

  # Check if cellranger-atac is available
  if ! command -v cellranger-atac >/dev/null 2>&1; then
    echo "[run_cellranger_atac] cellranger-atac not found in PATH." >&2
    echo "[run_cellranger_atac] Please ensure Cell Ranger ATAC is installed and in your PATH." >&2
    return 1
  fi

  # Create output directory
  mkdir -p "$output_dir"
  cd "$output_dir" || return 1

  if [[ "$use_sample_sheet" == true ]]; then
      # Batch mode: process samples using sample sheet
      echo "[run_cellranger_atac] Batch mode: processing samples from $sample_sheet"
      
      # Parse sample sheet (CSV format: sample_id,fastq_path)
      local samples=()
      local fastq_paths=()
      
      while IFS=',' read -r sample_id fastq_path || [[ -n "$sample_id" ]]; do
        # Skip empty lines and header
        if [[ -z "$sample_id" || "$sample_id" == "sample_id" ]]; then
          continue
        fi
        
        # Remove quotes if present
        sample_id=$(echo "$sample_id" | tr -d '"')
        fastq_path=$(echo "$fastq_path" | tr -d '"')
        
        # Use input_dir as base if fastq_path is relative
        if [[ "$fastq_path" != /* ]]; then
          fastq_path="$input_dir/$fastq_path"
        fi
        
        samples+=("$sample_id")
        fastq_paths+=("$fastq_path")
      done < "$sample_sheet"
      
      if [[ ${#samples[@]} -eq 0 ]]; then
        echo "[run_cellranger_atac] No valid samples found in $sample_sheet" >&2
        return 1
      fi
      
      echo "[run_cellranger_atac] Found ${#samples[@]} samples to process:"
      for i in "${!samples[@]}"; do
        echo "  ${samples[$i]} -> ${fastq_paths[$i]}"
      done
      
      # Process each sample
      local processed=0
      local failed=0
      
      for i in "${!samples[@]}"; do
        local sample_name="${samples[$i]}"
        local sample_fastqs="${fastq_paths[$i]}"
        
        echo "[run_cellranger_atac] Processing sample: $sample_name"
        
        # Check if FASTQ files exist
        if [[ ! -d "$sample_fastqs" ]]; then
          echo "[run_cellranger_atac] Warning: FASTQ directory not found: $sample_fastqs, skipping..."
          ((failed++))
          continue
        fi
        
        local fastq_count
        fastq_count=$(find "$sample_fastqs" -name "*.fastq.gz" | wc -l)
        
        if [[ $fastq_count -eq 0 ]]; then
          echo "[run_cellranger_atac] Warning: No FASTQ files found in $sample_fastqs, skipping..."
          ((failed++))
          continue
        fi
        
        echo "[run_cellranger_atac] Found $fastq_count FASTQ files in $sample_fastqs"
        
        # Run cellranger-atac for this sample
        if _run_single_cellranger_atac \
          "$sample_fastqs" \
          "$sample_name" \
          "$reference" \
          "$output_dir" \
          "$localcores" \
          "$localmem" \
          "$expect_cells" \
          "$chemistry" \
          "$nosecondary" \
          "$dry_run"; then
          echo "[run_cellranger_atac] Successfully processed: $sample_name"
          ((processed++))
        else
          echo "[run_cellranger_atac] Failed to process: $sample_name" >&2
          ((failed++))
        fi
      done
      
      echo "[run_cellranger_atac] Batch processing completed:"
      echo "[run_cellranger_atac]   Successfully processed: $processed samples"
      echo "[run_cellranger_atac]   Failed: $failed samples"
      
      if [[ $failed -gt 0 ]]; then
        return 1
      fi
    else
      # Batch mode: auto-discover samples from directory names
      echo "[run_cellranger_atac] Batch mode: auto-discovering samples in $input_dir"
      
      # Find all sample directories
      local sample_dirs
      mapfile -t sample_dirs < <(find "$input_dir" -maxdepth 1 -type d -name "*" | grep -v "^$input_dir$" | sort)
      
      if [[ ${#sample_dirs[@]} -eq 0 ]]; then
        echo "[run_cellranger_atac] No sample directories found in $input_dir" >&2
        return 1
      fi
      
      echo "[run_cellranger_atac] Found ${#sample_dirs[@]} sample directories:"
      for dir in "${sample_dirs[@]}"; do
        echo "  $(basename "$dir")"
      done
      
      # Process each sample
      local processed=0
      local failed=0
      
      for sample_dir in "${sample_dirs[@]}"; do
        local sample_name
        sample_name="$(basename "$sample_dir")"
        
        echo "[run_cellranger_atac] Processing sample: $sample_name"
        
        # Check if sample directory has FASTQ files
        local fastq_count
        fastq_count=$(find "$sample_dir" -name "*.fastq.gz" | wc -l)
        
        if [[ $fastq_count -eq 0 ]]; then
          echo "[run_cellranger_atac] Warning: No FASTQ files found in $sample_dir, skipping..."
          ((failed++))
          continue
        fi
        
        echo "[run_cellranger_atac] Found $fastq_count FASTQ files in $sample_dir"
        
        # Run cellranger-atac for this sample
        if _run_single_cellranger_atac \
          "$sample_dir" \
          "$sample_name" \
          "$reference" \
          "$output_dir" \
          "$localcores" \
          "$localmem" \
          "$expect_cells" \
          "$chemistry" \
          "$nosecondary" \
          "$dry_run"; then
          echo "[run_cellranger_atac] Successfully processed: $sample_name"
          ((processed++))
        else
          echo "[run_cellranger_atac] Failed to process: $sample_name" >&2
          ((failed++))
        fi
      done
      
      echo "[run_cellranger_atac] Batch processing completed:"
      echo "[run_cellranger_atac]   Successfully processed: $processed samples"
      echo "[run_cellranger_atac]   Failed: $failed samples"
      
      if [[ $failed -gt 0 ]]; then
        return 1
      fi
    fi
}

# Helper function to run cellranger-atac for a single sample
_run_single_cellranger_atac() {
  local fastqs="$1"
  local sample="$2"
  local reference="$3"
  local output_dir="$4"
  local localcores="$5"
  local localmem="$6"
  local expect_cells="$7"
  local chemistry="$8"
  local nosecondary="$9"
  local dry_run="${10}"

  # Build cellranger-atac command
  local cmd=(cellranger-atac count)
  cmd+=(--id="$sample")
  cmd+=(--fastqs="$fastqs")
  cmd+=(--reference="$reference")
  cmd+=(--localcores="$localcores")
  cmd+=(--localmem="$localmem")

  if [[ -n "$expect_cells" ]]; then
    cmd+=(--expect-cells="$expect_cells")
  fi

  if [[ "$chemistry" != "auto" ]]; then
    cmd+=(--chemistry="$chemistry")
  fi

  if [[ "$nosecondary" == true ]]; then
    cmd+=(--no-secondary)
  fi

  # Print command for debugging
  echo "[run_cellranger_atac] Command to run:"
  printf '%s ' "${cmd[@]}"
  echo

  if [[ "$dry_run" == true ]]; then
    echo "[run_cellranger_atac] Dry run mode - not executing command."
    return 0
  fi

  # Check FASTQ files
  echo "[run_cellranger_atac] Checking FASTQ files in: $fastqs"
  local fastq_count
  fastq_count=$(find "$fastqs" -name "*.fastq.gz" | wc -l)
  echo "[run_cellranger_atac] Found $fastq_count FASTQ files"

  if [[ $fastq_count -eq 0 ]]; then
    echo "[run_cellranger_atac] No FASTQ files found in $fastqs" >&2
    return 1
  fi

  # Show sample FASTQ files
  echo "[run_cellranger_atac] Sample FASTQ files:"
  find "$fastqs" -name "*.fastq.gz" | head -5

  # Run cellranger-atac
  echo "[run_cellranger_atac] Starting Cell Ranger ATAC analysis..."
  echo "[run_cellranger_atac] This may take several hours depending on data size."
  
  if "${cmd[@]}"; then
    echo "[run_cellranger_atac] Analysis completed successfully!"
    echo "[run_cellranger_atac] Output directory: $output_dir/$sample"
    
    # Check if output files exist
    local output_path="$output_dir/$sample"
    if [[ -d "$output_path" ]]; then
      echo "[run_cellranger_atac] Output files:"
      ls -la "$output_path"
      
      if [[ -f "$output_path/web_summary.html" ]]; then
        echo "[run_cellranger_atac] Web summary: $output_path/web_summary.html"
      fi
      
      if [[ -f "$output_path/summary.csv" ]]; then
        echo "[run_cellranger_atac] Summary CSV: $output_path/summary.csv"
      fi
      
      if [[ -d "$output_path/outs" ]]; then
        echo "[run_cellranger_atac] Main outputs in: $output_path/outs"
        ls -la "$output_path/outs"
      fi
    fi
    return 0
  else
    echo "[run_cellranger_atac] Analysis failed!" >&2
    return 1
  fi
}

_run_cellranger_atac_usage() {
  cat >&2 <<EOF
Usage: 
  Auto-discovery: run_cellranger_atac --input-dir INPUT_DIR --reference REF_DIR --output-dir OUT_DIR [options]
  Sample sheet:   run_cellranger_atac --input-dir INPUT_DIR --sample-sheet CSV_FILE --reference REF_DIR --output-dir OUT_DIR [options]

Required arguments:
  --input-dir INPUT_DIR            Directory containing sample subdirectories or base for FASTQ paths
  --reference REF_DIR              Cell Ranger ATAC reference directory
  --output-dir OUT_DIR             Output directory for results

Options:
  --sample-sheet CSV_FILE          CSV file with sample_id,fastq_path columns (optional)
  --localcores N                   Number of cores to use (default: 16)
  --localmem N                     Memory in GB to use (default: 64)
  --expect-cells N                 Expected number of cells (optional)
  --chemistry NAME                 Chemistry type (auto, ARC-v1, etc.)
  --no-secondary                   Disable secondary analysis
  --dry-run                        Show command without executing
  -h, --help                       Show this help

Batch Mode Options:
  1. Auto-discovery (default): Processes all subdirectories in INPUT_DIR
     - Sample names are inferred from subdirectory names
     - Each subdirectory should contain FASTQ files for one sample
     
  2. Sample sheet: Use --sample-sheet for custom sample/path mapping
     CSV Format:
     sample_id,fastq_path
     Sample1,/path/to/Sample1_fastqs
     Sample2,/path/to/Sample2_fastqs
     Sample3,Sample3_fastqs
     
     Notes:
     - Header row is optional but recommended
     - fastq_path can be absolute or relative to input-dir
     - Each row represents one sample to process

Outputs:
  OUT_DIR/SAMPLE_ID/               Main output directory
    ├── outs/                      Analysis results
    │   ├── filtered_peak_bc_matrix/     Peak matrix
    │   ├── filtered_tf_bc_matrix/       Transcription factor matrix
    │   ├── fragments.tsv.gz             Fragment file
    │   └── peak_annotation.tsv          Peak annotations
    ├── web_summary.html           Web-based summary report
    └── summary.csv                Summary statistics

Notes:
  - FASTQ files should be named according to Cell Ranger conventions
  - Reference must be compatible with Cell Ranger ATAC
  - Use --expect-cells for better cell calling
  - Check web_summary.html for quality metrics
  - Processes samples sequentially
EOF
}

# If executed directly, run the function with provided args
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  run_cellranger_atac "$@"
fi
