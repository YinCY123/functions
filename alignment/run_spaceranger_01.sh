run_spaceranger() {
  # Default values
  local transcriptome=""
  local fastqs_base=""
  local images_dir=""
  local loupe_dir=""
  local slide=""
  local area=""
  local cores=6
  local mem=30

  # Parse options
  while getopts ":t:f:i:l:s:a:c:m:h" opt; do
    case $opt in
      t) transcriptome="$OPTARG" ;;
      f) fastqs_base="$OPTARG" ;;
      i) images_dir="$OPTARG" ;;
      l) loupe_dir="$OPTARG" ;;
      s) slide="$OPTARG" ;;
      a) area="$OPTARG" ;;
      c) cores="$OPTARG" ;;
      m) mem="$OPTARG" ;;
      h)
        echo "Usage: $0 -t transcriptome -f fastqs_base -i images_dir -l loupe_dir -s slide -a area [-c cores] [-m mem]"
        return 0
        ;;
      \?)
        echo "Invalid option: -$OPTARG" >&2
        return 1
        ;;
      :)
        echo "Option -$OPTARG requires an argument." >&2
        return 1
        ;;
    esac
  done
  shift $((OPTIND -1))

  # Validate required parameters
  if [[ -z "$transcriptome" || -z "$fastqs_base" || -z "$images_dir" || -z "$loupe_dir" || -z "$slide" || -z "$area" ]]; then
    echo "Error: Missing required parameters."
    echo "Usage: $0 -t transcriptome -f fastqs_base -i images_dir -l loupe_dir -s slide -a area [-c cores] [-m mem]"
    return 1
  fi

  module load spaceranger

  # Detect samples from fastqs_base subdirectories
  local samples=()
  while IFS= read -r -d $'\0' dir; do
    samples+=("$(basename "$dir")")
  done < <(find "$fastqs_base" -mindepth 1 -maxdepth 1 -type d -print0)

  if [ ${#samples[@]} -eq 0 ]; then
    echo "No sample directories found in $fastqs_base"
    return 1
  fi

  for sample in "${samples[@]}"; do
    echo "Running Space Ranger count for sample: $sample"

    local fastqs="${fastqs_base}/${sample}"
    local image="${images_dir}/${sample}.tif"
    local loupe_alignment="${loupe_dir}/${sample}.json"

    if [[ ! -d "$fastqs" ]]; then
      echo "Warning: FASTQ directory not found for sample $sample: $fastqs"
      continue
    fi
    if [[ ! -f "$image" ]]; then
      echo "Warning: Image file not found for sample $sample: $image"
      continue
    fi
    if [[ ! -f "$loupe_alignment" ]]; then
      echo "Warning: Loupe alignment file not found for sample $sample: $loupe_alignment"
      continue
    fi

    spaceranger count \
      --id="${sample}" \
      --transcriptome="${transcriptome}" \
      --fastqs="${fastqs}" \
      --sample="${sample}" \
      --image="${image}" \
      --slide="${slide}" \
      --area="${area}" \
      --loupe-alignment="${loupe_alignment}" \
      --localcores="${cores}" \
      --localmem="${mem}" \
      --disable-ui

    echo "Finished sample: $sample"
    echo "-------------------------"
  done
}

# example: 
# run_spaceranger \
#   -t /path/to/refdata \
#   -f /path/to/fastqs \
#   -i /path/to/images \
#   -l /path/to/loupe \
#   -s V19J01-123 \
#   -a A1 \
#   -c 8 \
#   -m 32

