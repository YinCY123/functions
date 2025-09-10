#!/bin/bash

# samples=("SRR28746891" "SRR28746892" "SRR28746893" \
# 		"SRR28746894" "SRR28746895" "SRR28746896" \
# 		"SRR28746897" "SRR28746898" "SRR28746899" \
# 		"SRR28746900")
# base_dir="/home/yincy/disk14/qs/20250428/GSE264408/SAHMI"

samples=("SRR9141212" "SRR9141213" "SRR9141214" "SRR9141215" \
	"SRR9141216" "SRR9141217" "SRR9141218" "SRR9141219" \
	"SRR9141220" "SRR9141221")
base_dir="/home/yincy/disk14/partime/2025/20250418/GSE131882/SAHMI"


for sample in "${samples[@]}"; do
	echo "Processing $sample..."
	Rscript /home/yincy/git/github/SAHMI/functions/taxa_counts.r \
		--sample_name $sample \
		--fa1 $base_dir/${sample}_1.fa \
		--fa2 $base_dir/${sample}_2.fa \
		--taxa "taxa.tsv" \
		--kraken_report "$base_dir"/${sample}.kraken.report.txt \
		--mpa_report "$base_dir"/${sample}.kraken.report.mpa.txt \
		--out_path "$base_dir/" \
		--cb_len 16 \
		--umi_len 10
	done
