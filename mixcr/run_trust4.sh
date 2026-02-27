#!/bin/bash

run_trust4(){
    # help message
    if [[ "$1" == "--help" || "$1" == "-h" ]]; then
        cat << EOF
        Required parameters:
            -b                      path to bam file
            -1                      path to paired-end read files
            -2                      path to paired-end read files
            -u                      path to single-end read file
            -f                      path to the fasta file coordinate and sequence of V/D/J/C genes, default: hg38_bcrtcr.fa

        Optional parameters:
            --ref                   path to detailed V/D/J/C gene reference file from IMGT database, default: human_IMGT+C.fa
            -o                      prefix of output files
            --od                    the directory for output files, default: the current working directory
            -t                      number of threads, default: 10
            -k                      the starting k-mer size for indexing contig, default: 9
            --barcode               if -b, bam field for barcode; if -1, -2/-u, the file containing barcode, default: CR/r2
            --barcodeLevel          barcode is for cell or molecule, default: cell
            --barcodeTranslate      path to the barcode translate file, default: not used
            --UMI                   if -b, BAM field for UMI; if -1, -2/-, the file containing UMI, default: UR/r2
            --readFormat            format for read, barcode and UMI files (example: r1:0:-1,r2:0:-1,bc:0:15,um:16:-1 for paired-end files with barcode and UMI)
            --repseq                the data is from bulk, non-UMI-based TCR-seq or BCR-seq, default: not set
            --contigMinCov          ignore contigs that have bases covered by fewer that UMI INT reads, default: 0
            --minHitLen             the minimal hit length for a valid overlap, default: auto
            --mateIdSuffixLen       the suffix length in read id for mate, default: not used
            --skipMateExtension     do not extend assemblies with mate information, useful for SMART-seq, default: not used
            --skipReadRealign       do not realign reads in annotator, useful for reducing computation cost of barcode/UMI-based repseq, default: not used
            --cgenEnd               skip reads aligned after the first 200bp of the C gene, default: 200
            --abnormalUnmapFlag     the flag in BAM for the unmapped read-pair is nonconcordant, default: not set
            --noExtraction          directly use the files from provided -1 -2/-u to assemble, default: not set
            --imgtAdditionalGap     description for additional gap codon position in IMGT (0-based), e.g. "TRAV:7,83"
            --assembleWithRef       conduct the assemble with --ref file, default: use -f file
            --outputReadAssignment  output read assignment results to the prefix_assign.out file, default: no output
            --stage                 start TRUST4 on specified stage, default: 0
                                        0: start from begining, candidate read extraction
                                        1: start from assembly
                                        2: start from annotation
                                        3: start from generating the report table
            --clean                 clean up files. 0: no clean, 1: clean intermediate files, 2: only keep AIRR files, default: 0
        EOF

        return 0
    fi

    # initialize variable
    local bam
    local fq
    local bam_pattern="SRR[0-9]+"
    local fq_pattern="_1.fastq.gz"
    local f=/home/yincy/BioHome/datasets/trust4/human/hg38_bcrtcr.fa
    local ref=/home/yincy/BioHome/datasets/trust4/human/human_IMGT+C.fa
    local o=
    local od="./"
    local t=10
    local k=9
    local barcode=CB
    local barcodeLevel
    local barcodeWhitelist
    local barcodeTranslate
    local UMI=UB
    local readFormat
    local repseq
    local contigMinCov
    local minHitLen
    local mateIdSuffixLen
    local skipMateExtension
    local skipReadRealign
    local cgeneEnd
    local abnormalUnmapFlag
    local noExtraction
    local imgtAdditionalGap
    local assembleWithRef
    local outputReadAssignment
    local stage=0
    local clean=0

    if [[ -n "$bam" ]]; then
        echo "process from bam files..."

        # processing each bam files
        local bam_files=($( find "$bam" -type f -iname "*.bam" ))
        
        for sample in "${bam_files[@]}"; do
            sample_name=$( echo $sample | grep -oP "$bam_pattern" )

            run-trust4 \
                -b "$sample" \
                -f "$f" \
                --ref "$ref" \
                -o "$o"/"$sample_name"_ \
                --barcode "$barcode" \
                --UMI "$UMI"
        done
    fi

    if [[ -n "$fq" && "$paired" == "True" ]]; then
        local files=$( find "$fq" -type f \( -name "*_1.fastq.gz" -o -name "*_R1.fastq.gz" \) )

        for sample in "${files[@]}"; do
            sample_name=$(basename "$sample" "$fq_pattern")
            
            r1="$fq"/"$sample_name"/"$fq_pattern"
            r2="$fq"/"$sample_name"/"$( sed -i 's/1/2' $fq_pattern)"

            run-trust4 \
                -1 r1 \
                -2 r2 \
                -f "$f" \
                --ref "$ref" \
                --barcode "$barcode" \
                --UMI "$UMI" \
                --readFormat "$readFormat" \
                -o "$o"/"$sample_name"_
        done
    fi
}