#!/bin/bash

run_STARsolo(){
    # help message
    if [[ "$1" == "--help" || "$1" == "-h" ]]; then
        echo ""
        echo "Usage: run_STARsolo [options]"
        echo ""
        echo "Mandatory options:"
        echo "  --genomeDir, -g       Path to STAR genome index directory"
        echo "  --fastqDir, -q        Path to directory containing FASTQ files"
        echo "  --outputDir, -o       Path to output directory"
        echo ""
        echo "Optional parameters:"
        echo "  --seqType, -e         Sequencing type: single or paired (default: paired)"
        echo "  --barcodeLength, -b   Length of the barcode (default: 16)"
        echo "  --umiLength, -u       Length of the UMI (default: 12)"
        echo "  --soloType, -t        Type of solo analysis (default: CB_UMI_Simple)"
        echo "  --soloFeatures, -f    Features for solo analysis (default: GeneFull)"
        echo "  --soloBarcodeMate, -m Barcode mate (default: 0)"
        echo "  --soloBarcodeReadLength, -L, wether check read length and barcode length, (default, 1, check), 0 do not check"
        echo "  --ncores, -n          Number of threads (default: 12)"
        echo "  --runMode, -r         STAR run mode (default: alignReads)"
        echo "  --soloCBstart, -p     Start position of cell barcode (default: 1)"
        echo "  --soloCBwhitelist, -w Whitelist for cell barcodes (default: None)"
        echo "  --outSAMtype, -B      Output SAM type (default: BAM)"
        echo "  --outSAMsort, -S      Output SAM sort order (default: SortedByCoordinate)"
        echo "  --outSAMunmapped -I   Include unmapped reads in the BAM/SAM file (default: Within)"
        echo "  --limitBAMsortRAM, -M Limit for BAM sort RAM (default: 65G)"
        echo "  --BCread, -C          Which file contains barcode sequence (default: R1)"
        echo "  --clip5pNbases, -5    Number of bases to clip from 5' end of first mate (default: 0)"
        echo "  --clip3pNbases, -3    Number of bases to clip from 3' end of first mate (default: 0)"
        echo ""
        return 0
    fi

    # initialize variable
    local genomeDir=
    local fastqDir=
    local outputDir=
    local seqType=paired
    local barcodeLength=16
    local umiLength=12
    local soloType=CB_UMI_Simple
    local soloFeatures=GeneFull
    local soloBarcodeMate=0
    local soloBarcodeReadLength=1
    local ncores=12
    local runMode=alignReads
    local soloCBstart=1
    local soloCBwhitelist=None
    local outSAMtype=BAM
    local outSAMsort=SortedByCoordinate
    local outSAMunmapped=Within
    local limitBAMsortRAM=69793218560
    local BCread=R1
    local clip5pNbases=0
    local clip3pNbases=0

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
        --seqType|-e)
            seqType="$2"
            shift 2
            ;;
        --barcodeLength|-b)
            barcodeLength="$2"
            shift 2
            ;;
        --umiLength|-u)
            umiLength="$2"
            shift 2
            ;;
        --soloType|-t)
            soloType="$2"
            shift 2
            ;;
        --soloBarcodeReadLength|-L)
            soloBarcodeReadLength="$2"
            shift 2
            ;;
        --soloFeatures|-f)
            soloFeatures="$2"
            shift 2
            ;;
        --soloBarcodeMate|-m)
            soloBarcodeMate="$2"
            shift 2
            ;;
        --outSAMattributes|-A)
            outSAMattributes="$2"
            shift 2
            ;;
        --ncores|-n)
            ncores="$2"
            shift 2
            ;;
        --runMode|-r)
            runMode="$2"
            shift 2
            ;;
        --soloCBstart|-p)
            soloCBstart="$2"
            shift 2
            ;;
        --soloCBwhitelist|-w)
            soloCBwhitelist="$2"
            shift 2
            ;;
        --outSAMtype|-B)
            outSAMtype="$2"
            shift 2
            ;;
        --outSAMsort|-S)
            outSAMsort="$2"
            shift 2
            ;;
        --outSAMunmapped|-I)
            outSAMunmapped="$2"
            shift 2
            ;;
        --BCread|-C)
            BCread="$2"
            shift 2
            ;;
        --limitBAMsortRAM|-M)
            limitBAMsortRAM="$2"
            shift 2
            ;;
        --clip5pNbases|-5)
            clip5pNbases="$2"
            shift 2
            ;;
        --clip3pNbases|-3)
            clip3pNbases="$2"
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

    # check soloCBWhitelist exists
    if [[ "$soloCBwhitelist" != "None" ]]; then
        if [[ ! -f "$soloCBwhitelist" ]]; then
            echo "Error: Whitelist file $soloCBwhitelist not found"
            return 1
        fi
    fi

    # processing each sample
    for file in "$fastqDir"/*_1.fastq.gz; do
        local sampleName=$(basename "$file" _1.fastq.gz)
        local read1="$fastqDir/${sampleName}_1.fastq.gz"
        local read2="$fastqDir/${sampleName}_2.fastq.gz"
        local sampleOutputDir="$outputDir/$sampleName"
        local readFilesInParam

        # check if input files exist
        if [[ "$seqType" == "paired" ]]; then
            if [[ ! -f "$read1" ]] || [[ ! -f "$read2" ]]; then
                echo "Error: input files not found for sample $sampleName"
                continue
            fi
            if [[ "$BCread" == "R1" ]]; then
                readFilesInParam="--readFilesIn $read2 $read1"
            elif [[ "$BCread" == "R2" ]]; then
                readFilesInParam="--readFilesIn $read1 $read2"
            fi

        elif [[ "$seqType" == "single" ]]; then
            if [[ ! -f "$read1" ]]; then
                echo "Error: input file not found for sample $sampleName"
                continue
            fi
            readFilesInParam="--readFilesIn $read1"
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
        if [[ "$soloBarcodeMate" == "1" && "$clip5pNbases" == "0" ]]; then
            clip5pNbases="$barcodeLength 0"
        fi

        # run STARsolo
        STAR \
            --runThreadN "$ncores" \
            --runMode "$runMode" \
            --genomeDir "$genomeDir" \
            $readFilesInParam \
            --readFilesCommand zcat \
            --soloType "$soloType" \
            --soloCBstart 1 \
            --soloCBlen "$barcodeLength" \
            --soloUMIstart $((1 + $barcodeLength)) \
            --soloUMIlen "$umiLength" \
            --soloBarcodeReadLength "$soloBarcodeReadLength" \
            --soloFeatures "$soloFeatures" \
            --soloCBwhitelist "$soloCBwhitelist" \
            --outSAMattributes=NH HI AS nM NM MD jM jI MC ch XS CR CB CY UR UB UY \
            --outFileNamePrefix "$sampleOutputDir/" \
            --soloBarcodeMate "$soloBarcodeMate" \
            --clip5pNbases "$clip5pNbases" \
            --clip3pNbases "$clip3pNbases" \
            --outSAMtype $outSAMtype $outSAMsort \
            --limitBAMsortRAM "$limitBAMsortRAM" \
            --outSAMunmapped "$outSAMunmapped"
        
        # check if STAR completed successfully
        if [[ $? -ne 0 ]]; then
            echo "Error: STAR alignment failed for sample $sampleName"
            continue
        fi
    done
}
