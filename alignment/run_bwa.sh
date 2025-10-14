#!/bin/bash
run_bwa() {
  # Parse key=value parameters
  for arg in "$@"; do
    if [[ "$arg" == "-h" || "$arg" == "--help" ]]; then
      cat <<EOF
Usage: run_bwa SAMPLE=sample_name R1=read1.fastq.gz R2=read2.fastq.gz REF=reference.fa THREADS=num_threads [OUTDIR=output_dir] [CB_LEN=cell_barcode_length]

Required parameters:
  SAMPLE    Sample name (used for output files)
  R1        Read 1 FASTQ file (gzipped)
  R2        Read 2 FASTQ file (gzipped)
  REF       Reference genome FASTA file
  THREADS   Number of threads for alignment (default: 8)

Optional parameters:
  OUTDIR    Output directory (default: .)
  CB_LEN    Cell barcode length (default: 16)

Example:
  run_bwa SAMPLE=test R1=test_R1.fastq.gz R2=test_R2.fastq.gz REF=ref.fa THREADS=8 OUTDIR=results CB_LEN=16

This function:
  1. Extracts cell barcode from R1, trims it, and saves barcode info.
  2. Copies R2 unchanged.
  3. Aligns reads using BWA-MEM2.
  4. Sorts, indexes, marks duplicates, filters BAM.
  5. Adds CB tags to BAM based on the extracted barcodes.

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
      CB_LEN)   CB_LEN=$val ;;
      *) echo "Unknown parameter: $key"; return 1 ;;
    esac
  done

  : "${OUTDIR:=.}"
  : "${SAMPLE:?Missing SAMPLE}"
  : "${R1:?Missing R1}"
  : "${R2:?Missing R2}"
  : "${REF:?Missing REF}"
  : "${THREADS:=8}"
  : "${CB_LEN:=16}"

  mkdir -p "${OUTDIR}"

  # 1. Extract CB from R1, trim first CB_LEN bases, and save barcode info
  echo "[`date`] Extracting cell barcode..."
  BARCODE_MAP=${OUTDIR}/${SAMPLE}.r1_barcodes.tsv
  : > "${BARCODE_MAP}"

  # Trim R1 and write mapping (readname<TAB>barcode) to BARCODE_MAP
  zcat "${R1}" | awk -v cb_len=${CB_LEN} -v map="${BARCODE_MAP}" '
    NR%4==1{
      header=$0
      getline; seq=$0
      getline; plus=$0
      getline; qual=$0
      # extract read name (first token after @)
      h_no = header
      sub(/^@/,"",h_no)
      split(h_no, a, " ")
      readname = a[1]
      cb = substr(seq,1,cb_len)
      trimmed_seq = substr(seq, cb_len+1)
      trimmed_qual = substr(qual, cb_len+1)
      print header
      print trimmed_seq
      print plus
      print trimmed_qual
      print readname "\t" cb >> map
    }' | gzip > "${OUTDIR}/${SAMPLE}_R1.trim.fastq.gz"

  # R2 unchanged
  cp "${R2}" "${OUTDIR}/${SAMPLE}_R2.trim.fastq.gz"

  # 2. Alignment
  echo "[`date`] Aligning with BWA-MEM2..."
  bwa-mem2 mem -t ${THREADS} -M \
    -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}" \
    ${REF} \
    ${OUTDIR}/${SAMPLE}_R1.trim.fastq.gz ${OUTDIR}/${SAMPLE}_R2.trim.fastq.gz \
    | samtools view -b -@ ${THREADS} -o ${OUTDIR}/${SAMPLE}.unsorted.bam

  # 3. Sorting, indexing, deduplication, filtering (unchanged)
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

  # 4. Add CB tags to BAM using the barcode mapping file
  # The mapping file contains: <readname>\t<barcode>
  echo "[`date`] Adding CB tags to BAM..."
  samtools view -h "${OUTDIR}/${SAMPLE}.dedup.bam" | \
    awk -v map="${BARCODE_MAP}" '
      BEGIN {
        while ((getline < map) > 0) {
          split($0, f, "\t")
          bar[f[1]] = f[2]
        }
      }
      /^@/ { print; next }
      {
        qname = $1
        # remove /1 or /2 suffix if present
        sub(/\/[12]$/, "", qname)
        if (qname in bar) {
          print $0 "\tCB:Z:" bar[qname]
        } else {
          print
        }
      }' | samtools view -b -@ ${THREADS} -o "${OUTDIR}/${SAMPLE}.dedup.cb.bam" -

  samtools index "${OUTDIR}/${SAMPLE}.dedup.cb.bam"

  echo "[`date`] Filtering BAM..."
  samtools view -b -q 10 "${OUTDIR}/${SAMPLE}.dedup.cb.bam" > ${OUTDIR}/${SAMPLE}.q10.bam
  samtools index ${OUTDIR}/${SAMPLE}.q10.bam

  echo "[`date`] Done. Final BAM with CB tags: ${OUTDIR}/${SAMPLE}.q10.bam"
}