#!/bin/bash

# samples=("SRR9141212" "SRR9141213" "SRR9141214" "SRR9141215" \
# 	"SRR9141216" "SRR9141217" "SRR9141218" "SRR9141219" \
# 	"SRR9141220" "SRR9141221")
# samples="SRR9141219"
# paired=TRUE
# base_dir="/home/yincy/disk14/partime/2025/20250418/GSE131882/SAHMI"

# paired=TRUE
# base_dir="/home/yincy/disk14/qs/20250428/GSE264408/SAHMI"
# samples=("SRR28746891" "SRR28746892" "SRR28746893" \
# 		"SRR28746894" "SRR28746895" "SRR28746896" \
# 		"SRR28746897" "SRR28746898" "SRR28746899" \
# 		"SRR28746900")

fq_dir="/home/yincy/disk14/yr/20250223-WX-耀联/SRA/GSE246779/fastq"
output_dir="/home/yincy/disk14/yr/20250223-WX-耀联/SRA/GSE246779/SAHMI"

# samples=("SRR9141212" "SRR9141213" "SRR9141214" "SRR9141215" "SRR9141216" "SRR9141217" \
#     "SRR9141218" "SRR9141219" "SRR9141220" "SRR9141221")

# remove SRR26630464, SRR26630485, SRR26630482
samples=("SRR26630462" "SRR26630463" "SRR26630465" "SRR26630466" \
	"SRR26630467" "SRR26630480" "SRR26630481" "SRR26630483" \
	"SRR26630484")

paired="TRUE"


for sample in "${samples[@]}"; do
	if [[ "$paired" == TRUE ]]; then
		echo ""
		echo "processing paired sample: $sample"
		# echo ""
		Rscript /home/yincy/git/github/SAHMI/functions/sckmer.r \
			--sample_name $sample \
			--fa1 ${base_dir}/${sample}_1.fa \
			--fa2 ${base_dir}/${sample}_2.fa \
			--microbiome_output_file ${base_dir}/${sample}.microbiome.output.txt \
			--kraken_report ${base_dir}/${sample}.kraken.report.txt \
			--mpa_report ${base_dir}/${sample}.kraken.report.mpa.txt \
			--out_path "$base_dir/" \
			--cb_len 16 \
			--umi_len 10 \
			--ranks 'G'
	else
		echo "processing single end sample: $sample"
		Rscript /home/yincy/git/github/SAHMI/functions/sckmer_unpaired.r \
			--sample_name $sample \
			--fa1 ${base_dir}/${sample}_1.fa \
			--microbiome_output_file ${base_dir}/${sample}.microbiome.output.txt \
			--kraken_report ${base_dir}/${sample}.kraken.report.txt \
			--mpa_report ${base_dir}/${sample}.kraken.report.mpa.txt \
			--out_path $base_dir \
			--cb_len 16 \
			--umi_len 10 \
			--ranks 'G'
	fi
done