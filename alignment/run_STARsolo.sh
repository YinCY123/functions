#!/bin/bash

run_STARsolo(){
    # help message
    if [[ "$1" == "--help" || "$1" == "-h" ]]; then
        echo ""
        echo "Usage: run_STARsolo [options]"
        echo ""
        echo "Mandatory options:"
        echo "  --genomeDir, -g                 Path to STAR genome index directory"
        echo "  --fastqDir, -q                  Path to directory containing FASTQ files"
        echo "  --outputDir, -o                 Path to output directory"
        echo ""
        echo "Optional parameters:"
        echo "  --seqType                       sequencing type: single or paired (default: paired)"
        echo "  --ncores                        number of threads (default: 12)"
        echo "  --runMode                       STAR run mode (default: alignReads)"
        echo "  --BCread                        which file contains barcode sequence (default: R1)"
        echo "  --clipAdapterType               clip adapter method. (default: Hamming, CellRanger4 for 10x 3' v2)"
        echo "  --clip3pNbases                  number of bases to clip from 3' end of first mate (default: 0)"
        echo "  --clip3pAdapterSeq              sequence to cliped off from 3p"
        echo "  --clip5pNbases                  number of bases to clip from 5' end of first mate (default: 0)"
        echo "  --soloType                      the structure of cell barcode and umi (default: CB_UMI_Simple)"
        echo "  --soloCBtype                    cell barcode types (default: Sequence)"
        echo "  --soloCBwhitelist               cell barcode white list(s), comma-separated or multiple --soloCBwhitelist options (default: None)"
        echo "  --soloCBstart                   the start position of cell barcode (default: 1)"
        echo "  --soloCBlen                     length of cell barcode (deafult: 16)"
        echo "  --soloUMIstart                  start position of UMI (deafult: 17)"
        echo "  --soloUMIlen                    length of UMI (default: 12)"
        echo "  --soloCBposition                position of cell barcode, only used when --soloType=CB_UMI_Complex"
        echo "  --soloUMIposition               position of umi, only used when --soloType=CB_UMI_Complex"
        echo "  --soloAdapterSequence           adapter sequence to anchor barcodes (deafult: -)"
        echo "  --soloAdapterMismatchesNmax     maximum mismatch of adapter seqence (default: 1)"
        echo "  --soloCBmatchWLtype             the maximum mismatch of cell barcode to the white list (default: 1MM_multi)"
        echo "  --soloBarcodeReadLength         wether check read length and barcode length, (default, 1, check), 0 do not check"
        echo "  --soloBarcodeMate               where the cell barcode (default: 0, cell barcode in the second file of --readFilesIn)"
        echo "  --soloStrand                    which strand does the cDNA comes from (deafult: Forward, 3': Forward; 5' Reverse)"
        echo "  --soloFeatures                  features for solo analysis (default: GeneFull)"
        echo "  --soloMultiMappers              counting method for reads mapping to multiple genes (deafult: EM)"
        echo "  --soloUMIdedup                  feature counting for umi (default: 1MM_CR, 1 mismatch distance UMIs are collapsed)"
        echo "  --soloUMIfiltering              UMI filtering (default:  remove UMIs with N and homopolymers and low count UMI that mapped to multiple genes)"
        echo "  --soloCellFilter                cell filtering method (default: EmptyDrops_CR, use algorithms similar to DropemptyDrops)"
        echo "  --soloCellReadStats             output reads stat for each cell barcode"

        echo "  --outSAMunmapped                include unmapped reads in the BAM/SAM file (default: Within)"
        echo "  --limitBAMsortRAM               limit for BAM sort RAM (default: 75G)"
        echo "  --suppressBAM                   suppress SAM/BAM output (default: 0: false; 1: true)"

        echo "  --readFilesCommand              how to read files (default: zcat, for read gzip compressed file) use bzcat for read bz2 files"
        echo ""
        return 0
    fi

    # initialize variable
    local genomeDir=
    local fastqDir=
    local outputDir=

    local seqType=paired
    local ncores=12
    local runMode=alignReads
    local BCread=R1

    local clipAdapterType=CellRanger4
    local outFilterScoreMin=30
    local clip3pNbases=0
    local clip3pAdapterSeq=-
    local clip5pNbases=0

    local soloType=CB_UMI_Simple
    local soloCBtype=Sequence
    local soloCBwhitelist=None
    local soloCBstart=1
    local soloCBlen=16
    local soloUMIstart=17
    local soloUMIlen=12
    local soloCBposition=-
    local soloUMIposition=()
    local soloUMIposition_str="-"
    local soloAdapterSequence=-
    local soloAdapterMismatchesNmax=1
    local soloCBmatchWLtype=1MM
    local soloBarcodeReadLength=0
    local soloBarcodeMate=0
    local soloStrand=Forward
    local soloFeatures=GeneFull
    local soloMultiMappers=EM
    local soloUMIdedup=1MM_CR
    local soloUMIfiltering=MultiGeneUMI_CR
    local soloCellFilter=EmptyDrops_CR
    local soloCellReadStats=None

    local suppressBAM=0
    local outSAMtype=BAM
    local limitBAMsortRAM=75000000000
    local outSAMattributes=(NH HI AS nM NM MD jM jI MC ch XS CR CB CY UR UB UY)

    local readFilesCommand=zcat

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
        --seqType)
            seqType="$2"
            shift 2
            ;;
        --soloCBlen)
            soloCBlen="$2"
            shift 2
            ;;
        --soloUMIstart)
            soloUMIstart="$2"
            shift 2
            ;;
        --soloUMIlen)
            soloUMIlen="$2"
            shift 2
            ;;
        --soloType)
            soloType="$2"
            shift 2
            ;;
        --soloBarcodeReadLength)
            soloBarcodeReadLength="$2"
            shift 2
            ;;
        --soloFeatures)
            soloFeatures="$2"
            shift 2
            ;;
        --soloCellFilter)
            soloCellFilter="$2"
            shift 2
            ;;
        --soloBarcodeMate)
            soloBarcodeMate="$2"
            shift 2
            ;;
        --outSAMattributes)
            outSAMattributes=($2)
            shift 2
            ;;
        --ncores)
            ncores="$2"
            shift 2
            ;;
        --runMode)
            runMode="$2"
            shift 2
            ;;
        --soloCBstart)
            soloCBstart="$2"
            shift 2
            ;;
        --soloCBmatchWLtype)
            soloCBmatchWLtype="$2"
            shift 2
            ;;
        --soloCBposition)
            soloCBposition="$2"
            shift 2
            ;;
        --soloUMIposition)
            # Support multiple occurrences for CB_UMI_Complex mode
            value=$(echo "$2" | xargs)  # trim whitespace
            if [[ -n "$value" ]]; then
                soloUMIposition+=("$value")
            fi
            shift 2
            ;;
        --soloCBwhitelist)
            soloCBwhitelist="$2"
            shift 2
            ;;
        --soloStrand)
            soloStrand="$2"
            shift 2
            ;;
        --outSAMtype)
            outSAMtype="$2"
            shift 2
            ;;
        --outSAMsort)
            outSAMsort="$2"
            shift 2
            ;;
        --outSAMunmapped)
            outSAMunmapped="$2"
            shift 2
            ;;
        --BCread)
            BCread="$2"
            shift 2
            ;;
        --limitBAMsortRAM)
            limitBAMsortRAM="$2"
            shift 2
            ;;
        --clip5pNbases)
            clip5pNbases="$2"
            shift 2
            ;;
        --clip3pNbases)
            clip3pNbases="$2"
            shift 2
            ;;
        --suppressBAM)
            suppressBAM="$2"
            shift 2
            ;;
        --clipAdapterType)
            clipAdapterType="$2"
            shift 2
            ;;
        *)
            echo "Unknown parameter: $1"
            return 1
            ;;
        esac
    done

    # validate mandatory parameters
    if [[ -z "$genomeDir" || -z "$fastqDir" || -z "$outputDir" ]]; then 
        echo "please set the mandatory options."
        return 1
    fi

    # check if output directory exists
    if [[ -d "$outputDir" ]]; then
        echo "$outputDir already exists."
    else
        mkdir -p "$outputDir"
    fi

    # Prepare soloUMIposition string for STAR (space-separated)
    # if [[ ${#soloUMIposition[@]} -gt 0 && "${soloUMIposition[0]}" != "-" ]]; then
    #     soloUMIposition_str=$(IFS=' '; echo "${soloUMIposition[*]}")
    # fi

    # processing each sample
    local fastq_files=( "$fastqDir"/*_R1.fastq.gz )
    if [[ ${#fastq_files[@]} -eq 0 ]]; then
        echo "No FASTQ files matching *_R1.fastq.gz found in $fastqDir"
        return 1
    fi

    for file in "${fastq_files[@]}"; do
        local sampleName
        sampleName=$(basename "$file" _R1.fastq.gz)
        local read1="$fastqDir/${sampleName}_R1.fastq.gz"
        local read2="$fastqDir/${sampleName}_R2.fastq.gz"
        local sampleOutputDir="$outputDir/$sampleName"
        local readFilesInParamArr=()

        # check if input files exist
        if [[ "$seqType" == "paired" ]]; then
            if [[ ! -f "$read1" ]] || [[ ! -f "$read2" ]]; then
                echo "Error: input files not found for sample $sampleName"
                continue
            fi
            if [[ "$BCread" == "R1" ]]; then
                readFilesInParamArr=(--readFilesIn "$read2" "$read1")
            elif [[ "$BCread" == "R2" ]]; then
                readFilesInParamArr=(--readFilesIn "$read1" "$read2")
            fi

        elif [[ "$seqType" == "single" ]]; then
            if [[ ! -f "$read1" ]]; then
                echo "Error: input file not found for sample $sampleName"
                continue
            fi
            readFilesInParamArr=(--readFilesIn "$read1")
        else
            echo "Error: Invalid sequencing type specified: $seqType"
            return 1
        fi

        # check if sample output directory exists
        if [[ -d "$sampleOutputDir" ]]; then
            echo "$sampleOutputDir already exists."
        else
            mkdir -p "$sampleOutputDir"
        fi

        # Set default clipping values based on barcode mate if not specified
        # if [[ "$soloBarcodeMate" == "1" && "$clip5pNbases" == "0" ]]; then
        #     clip5pNbases="$barcodeLength 0"
        # fi

        # prepare outSAMtype parameter (allow suppressing BAM/SAM)
        local outSAMtypeParamArr
        if [[ "$suppressBAM" -eq 1 ]]; then
            outSAMtypeParamArr=(--outSAMtype None)
        else
            # allow user override of outSAMunmapped, default to Within
            outSAMtypeParamArr=(--outSAMtype "$outSAMtype" SortedByCoordinate \
                --outSAMunmapped "${outSAMunmapped:-Within}" \
                --limitBAMsortRAM "$limitBAMsortRAM" \
                --outSAMattributes "${outSAMattributes[@]}")
        fi

        # run STARsolo
        STAR \
            --runThreadN "$ncores" \
            --runMode "$runMode" \
            --genomeDir "$genomeDir" \
            "${readFilesInParamArr[@]}" \
            --readFilesCommand "$readFilesCommand" \
            --soloType "$soloType" \
            --soloCBstart "$soloCBstart" \
            --soloCBlen "$soloCBlen" \
            --soloUMIstart "$soloUMIstart" \
            --soloUMIlen "$soloUMIlen" \
            --soloBarcodeReadLength "$soloBarcodeReadLength" \
            --soloFeatures "$soloFeatures" \
            --soloCBwhitelist "$soloCBwhitelist" \
            --soloCBposition "$soloCBposition" \
            --soloUMIposition "$soloUMIposition_str" \
            --outFileNamePrefix "$sampleOutputDir/" \
            --soloBarcodeMate "$soloBarcodeMate" \
            --soloCBmatchWLtype "$soloCBmatchWLtype" \
            --clip5pNbases "$clip5pNbases" \
            --clip3pNbases "$clip3pNbases" \
            --soloStrand "$soloStrand" \
            --soloCellFilter "$soloCellFilter" \
            "${outSAMtypeParamArr[@]}" \
            --clipAdapterType="$clipAdapterType"

        # check if STAR completed successfully
        if [[ $? -ne 0 ]]; then
            echo "Error: STAR alignment failed for sample $sampleName"
            continue
        fi
    done
}
