#!/bin/bash
# Define a function to process scATAC FASTQs into a BAM
# Usage:
#   run_bwa SAMPLE=sample1 R1=sample1_R1.fastq.gz R2=sample1_R2.fastq.gz REF=hg38.fa THREADS=16

run_bwa() {

  # -----------------------------
  # Parse key=value parameters
  # -----------------------------
  for arg in "$@"; do
    # Check for help flags
    if [[ "$arg" == "-h" || "$arg" == "--help" ]]; then
      cat << 'EOF'
run_bwa - Process scATAC-seq FASTQ files into a filtered BAM file

Usage:
  run_bwa SAMPLE=<sample_name> R1=<R1.fastq.gz> R2=<R2.fastq.gz> REF=<reference.fa> THREADS=<num_threads> OUTDIR=<output_directory>

Description:
  This function performs a typical scATAC-seq alignment pipeline:
    1. Adapter trimming using cutadapt.
    2. Alignment to the reference genome using BWA-MEM2.
    3. BAM sorting and indexing with samtools.
    4. Duplicate marking with Picard.
    5. Filtering BAM for mapping quality >= 30.

Parameters:
  SAMPLE    - Sample name (used for output file naming and read group).
  R1        - Path to read 1 FASTQ file (gzipped).
  R2        - Path to read 2 FASTQ file (gzipped).
  REF       - Path to reference genome FASTA file.
  THREADS   - Number of threads to use for parallel processing.
  OUTDIR    - Output directory for all result files (default: current directory).

Outputs:
  <OUTDIR>/<SAMPLE>_R1.trim.fastq.gz, <OUTDIR>/<SAMPLE>_R2.trim.fastq.gz   - Adapter-trimmed FASTQ files.
  <OUTDIR>/<SAMPLE>.sorted.bam                                            - Sorted BAM file.
  <OUTDIR>/<SAMPLE>.dedup.bam                                             - Deduplicated BAM file.
  <OUTDIR>/<SAMPLE>.q10.bam                                               - Final filtered BAM file (MAPQ >= 10).
  <OUTDIR>/<SAMPLE>.q10.bam.bai                                           - BAM index file.
  <OUTDIR>/<SAMPLE>.dup_metrics.txt                                       - Duplicate marking metrics.

Example:
  run_bwa SAMPLE=sample1 R1=sample1_R1.fastq.gz R2=sample1_R2.fastq.gz REF=hg38.fa THREADS=16 OUTDIR=/path/to/output

Notes:
  - Requires cutadapt, bwa-mem2, samtools, and Picard.
  - All output files are saved to the specified OUTDIR directory.
  - If OUTDIR is not specified, files are saved to the current directory.
EOF
      return 0
    fi
    
    key=$(echo $arg | cut -f1 -d=)
    val=$(echo $arg | cut -f2- -d=)
    case "$key" in
      SAMPLE)   SAMPLE=$val ;;
      R1)       R1=$val ;;
      R2)       R2=$val ;;
      REF)      REF=$val ;;
      THREADS)  THREADS=$val ;;
      OUTDIR)   OUTDIR=$val ;;
      *) echo "Unknown parameter: $key"; return 1 ;;
    esac
  done

  # Set default OUTDIR if not provided
  : "${OUTDIR:=.}"

  # Check required params
  : "${SAMPLE:?Missing SAMPLE}"
  : "${R1:?Missing R1}"
  : "${R2:?Missing R2}"
  : "${REF:?Missing REF}"
  : "${THREADS:?Missing THREADS}"

  # Create output directory if it doesn't exist
  mkdir -p "${OUTDIR}"

  # -----------------------------
  # Pipeline
  # -----------------------------
  echo "[`date`] Trimming adapters..."
  cutadapt -j ${THREADS} \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -o ${OUTDIR}/${SAMPLE}_R1.trim.fastq.gz \
    -p ${OUTDIR}/${SAMPLE}_R2.trim.fastq.gz \
    ${R1} ${R2}

  echo "[`date`] Aligning with BWA-MEM2..."
  bwa-mem2 mem -t ${THREADS} -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}" \
    ${REF} \
    ${OUTDIR}/${SAMPLE}_R1.trim.fastq.gz ${OUTDIR}/${SAMPLE}_R2.trim.fastq.gz \
    | samtools view -b -@ ${THREADS} -o ${OUTDIR}/${SAMPLE}.unsorted.bam

  echo "[`date`] Sorting BAM..."
  samtools sort -@ ${THREADS} -o ${OUTDIR}/${SAMPLE}.sorted.bam ${OUTDIR}/${SAMPLE}.unsorted.bam
  samtools index ${OUTDIR}/${SAMPLE}.sorted.bam
  rm ${OUTDIR}/${SAMPLE}.unsorted.bam

  echo "[`date`] Marking duplicates..."
  /home/yincy/miniconda3/bin/picard MarkDuplicates \
    I=${OUTDIR}/${SAMPLE}.sorted.bam \
    O=${OUTDIR}/${SAMPLE}.dedup.bam \
    M=${OUTDIR}/${SAMPLE}.dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    VALIDATION_STRINGENCY=SILENT
  samtools index ${OUTDIR}/${SAMPLE}.dedup.bam

  echo "[`date`] Filtering BAM..."
  samtools view -b -q 10 ${OUTDIR}/${SAMPLE}.dedup.bam > ${OUTDIR}/${SAMPLE}.q10.bam
  samtools index ${OUTDIR}/${SAMPLE}.q10.bam

  echo "[`date`] Done. Final BAM: ${OUTDIR}/${SAMPLE}.q10.bam"
}
