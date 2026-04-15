#!/bin/bash

run_STAR(){
    # help message
    if [[ "$1" = "--help" || "$1" = "-h" ]]; then
        echo ""
        echo "Usage: run_STAR [options]"
        echo "Mandatory options:"
        echo "    --genomeDir, -g       Path to STAR genome index directory"
        echo "    --fastqDir, -q        Path to directory containing FASTQ files"
        echo "    --outputDir, -o       Path to output directory"
        echo "Optional parameters:"
        echo "    --nthreads, -n        Number of threads (default: 12)"
        echo "    --runMode, -r         STAR run mode (default: alignReads)"
        echo "    --quantMode, -Q       quantification mode"
        echo "    --seq_type, -T        sequence layout, pair or single (default: pair)"
        echo "    --trim3p, -c          trim number of bases from 3' (default: 0)"
        echo "    --trim5p, -C          trim number of bases from 5' (default: 0)"
        echo "    --outSAMtype, -S      output SAM type, None do not produce SAM/BAM file"
        echo "    --outSAMsort, -R      whether sort BAM file (default: SortedByCoordinate)"
        echo " "
        return 0
    fi

    # inititalize variable
    local genomeDir=
    local fastqDir=
    local outputDir=
    local runMode=alignReads
    local quantMode=GeneCounts
    local outSAMtype=BAM
    local outSAMsort=SortedByCoordinate
    local seq_type=pair
    local trim3p=0
    local trim5p=0
    local nthreads=12


    while [[ $# -gt 0 ]]; do
        case $1 in
            --genomeDir|-g)
                genomeDir="$2"
                shift 2
                ;;
            --fastqDir|-q)
                fastqDir="$2"
                shift 2
                ;;
            --outputDir|-o)
                outputDir="$2"
                shift 2
                ;;
            --nthreads|-n)
                nthreads="$2"
                shift 2
                ;;
            --runMode|-r)
                runMode="$2"
                shift 2
                ;;
            --quantMode|-Q)
                quantMode="$2"
                shift 2
                ;;
            --outSAMtype|-S)
                outSAMtype="$2"
                shift 2
                ;;
            --outSAMsort|-R)
                outSAMsort="$2"
                shift 2
                ;;
            --seq_type|-T)
                seq_type="$2"
                shift 2
                ;;
            --trim3p|-c)
                if [[ "$seq_type" == "pair" ]]; then
                    if [[ -n "$3" ]]; then
                        trim3p="$2 $3"
                        shift 3
                    else
                        trim3p="$2 $2"  # Use same value for both reads
                        shift 2
                    fi
                else
                    trim3p="$2"  # Single value for single-end
                    shift 2
                fi
                ;;
            --trim5p|-C)
                if [[ "$seq_type" == "pair" ]]; then
                    if [[ -n "$3" ]]; then
                        trim5p="$2 $3"
                        shift 3
                    else
                        trim5p="$2 $2"  # Use same value for both reads
                        shift 2
                    fi
                else
                    trim5p="$2"  # Single value for single-end
                    shift 2
                fi
                ;;                
        esac
    done

    # process each sample
    for sample in "$fastqDir"/*_R1.fastq.gz; do
        local sample_name=$(basename "$sample" _R1.fastq.gz)
        local read1="$fastqDir/${sample_name}_R1.fastq.gz"
        local read2="$fastqDir/${sample_name}_R2.fastq.gz"
        local sampleDir="$outputDir/$sample_name"
        local readFilesIn=""

        if [[  -d "$sampleDir" ]]; then
            echo "Directory $sampleDir already exists."
        else
            mkdir -p "$sampleDir"
        fi

        # prepare readin files
        if [[ "$seq_type" == "pair" ]]; then
            # if [[ -f "$read2" ]]; then
            readFilesIn="--readFilesIn $read1 $read2"
        elif [[ "$seq_type" == "single" ]]; then
            readFilesIn="--readFilesIn $read1"
        # fi
        # elif [[ "$seq_type" == "single" ]]; then
        #     readFilesIn="--readFilesIn $read1"
        else 
            echo "sequence type should be pair or single."
            exit 1
        fi

        # run STAR
        STAR \
            --runMode "$runMode" \
            --runThreadN "$nthreads" \
            ${readFilesIn} \
            --readFilesCommand zcat \
            --quantMode "$quantMode" \
            --genomeDir "$genomeDir" \
            --outFileNamePrefix "$sampleDir/" \
            --outSAMtype "$outSAMtype" "$outSAMsort" \
            --clip3pNbases "$trim3p" \
            --clip5pNbases "$trim5p"
    done
}