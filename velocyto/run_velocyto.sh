#!/bin/bash

# this script is used to process starsolo produced bam files in it's native location

run_velocyto_starsolo(){

    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    cat <<EOF
    --sample_dir          STARsolo prduced sample directory.
    --bcfile              cell barcode file.
    --output_dir          output loom file directory, default to the sample_dir.
    --sample_id            the sample name that will be used to retrieve information form metadatatable.
    --metadatatable       table comtaining metadata of the various samples (csv format, rows are samples and cols are variables).
    --mask                a gtf file containin ggenomic intervals to mask.
    --onefilepercell      if this flag is used every bam file passed is interpreted as an independent cell, otherwise multiple files are interpreted as batch of different cells to analyzed together.
    --logical             The logic to use for the filtering, default: Default
    --without_umi         If this flag is used the data is assumed UMI-less and reads are counted instead of molecules (default: off).
    --umi_extension       In case UMI is too short to gurantee uniqueness (without information from the ampping) set this parameter to chr, Gene or [N]bp.
    --multimap            Consider not unique mappings (not recommended).
    --samtools_threads    The number of threads to use to sort the bam by cellID file using samtools.
    --samtools_memory     The number of MB used for every threads by samtools to sort the bam file.
    --dtype               the dtype of the loom fiel layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation (default: uint32).
EOF
    return 0
    fi

    local sample_dir
    local gtf

    local output_dir
    local mask
    local logical=${logical:-Default}
    local umi_extension=${umi_extension:-no}
    local samtools_threads=${samtools_threads:-10}
    local samtools_memory=${samtools_memory:-1000}
    local dtype=${dtype:-uint32}


    while [[ $# -ge 1 ]]; do
        case $1 in
            --sample_dir)
                sample_dir="$2"
                shift 2
                ;;
            --gtf)
                gtf="$2"
                shift 2
                ;;
            --bcfile)
                bcfile="$2"
                shift 2
                ;;
            --output_dir)
                output_dir="$2"
                shift 2
                ;;
            --sample_id)
                sample_id="$2"
                shift 2
                ;;
            --metadatatable)
                metadatatable="$2"
                shift 2
                ;;
            --mask)
                mask="$2"
                shift 2
                ;;
            --logical)
                logical="$2"
                shift 2
                ;;
            --without_umi)
                without_umi="$2"
                shift 2
                ;;
            --umi_extension)
                umi_extension="$2"
                shift 2
                ;;
            --samtools_threads)
                samtools_threads="$2"
                shift 2
                ;;
            --samtools_memory)
                samtools_memory="$2"
                shift 2
                ;;
            --dtype)
                dtype="$2"
                shift 2
                ;;
            --verbose)
                verbose="$2"
                shift 2
                ;;
            *)
                echo "unknown parameter: $1"
                return 1
                ;;
        esac
    done

    if [ ! -d "$sample_dir" ]; then
        echo "the $sample_dir is not available."
        return 1
    fi

    bam_files=( $( find "$sample_dir" -type f -name "*.bam" ) )

    # processing each sample
    for bam in "${bam_files[@]}"; do
        echo "Processing $bam"
        output_dir=$( dirname "$bam" )
        velocyto run \
            -o $output_dir \
            -m $mask \
            -u $umi_extension \
            -@ $samtools_threads \
            --samtools-memory $samtools_memory \
            -t $dtype \
            "$bam" \
            "$gtf" \
        
        echo "Finished $bam processing..."
    done
}