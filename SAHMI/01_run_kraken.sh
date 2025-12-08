#!/bin/bash

# output_dir="/home/yincy/disk14/qs/20250428/GSE264408/SAHMI/" # do not forget the last back slash
# fq_dir="/home/yincy/disk14/qs/20250428/GSE264408/fastq"

# samples=("SRR28746891" "SRR28746892" "SRR28746893" \
# 		"SRR28746894" "SRR28746895" "SRR28746896" \
# 		"SRR28746897" "SRR28746898" "SRR28746899" \
# 		"SRR28746900")

fq_dir="/home/yincy/disk14/partime/2025/20250418/GSE131882/fastq"
output_dir="/home/yincy/disk14/partime/2025/20250418/GSE131882/SAHMI/"
samples=("SRR9141212" "SRR9141213" "SRR9141214" "SRR9141215" \
"SRR9141216" "SRR9141217" "SRR9141218" "SRR9141219" \
"SRR9141220" "SRR9141221")


for sample in "${samples[@]}"; do
	echo "processing $sample..."
	Rscript /home/yincy/git/github/SAHMI/functions/run_kraken.r \
		--sample "$sample" \
		--fq1 $fq_dir/${sample}_1.fastq.gz \
		--fq2 $fq_dir/${sample}_2.fastq.gz \
		--out_path $output_dir \
		--ncbi_blast_path /usr/lib/ncbi-blast+ \
		--Kraken2Uniq_path /home/yincy/bin/kraken2 \
		--kraken_database_path /home/yincy/ssd/kraken2/database/k2_minusb_20250402 \
		--paired T \
		--kreport2mpa_path /home/yincy/git/github/SAHMI/functions/kreport2mpa.py
done


