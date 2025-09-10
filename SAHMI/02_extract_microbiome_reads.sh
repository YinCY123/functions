#!/bin/bash

base_dir="/home/yincy/disk14/partime/2025/20250418/GSE131882/SAHMI"
samples=("SRR9141212" "SRR9141213" "SRR9141214" "SRR9141215" "SRR9141216" \
	"SRR9141217" "SRR9141218" "SRR9141219" "SRR9141220" "SRR9141221")
# samples="SRR9141219"
# base_dir="/home/yincy/disk14/qs/20250428/GSE264408/SAHMI"
# samples=("SRR28746891" "SRR28746892" "SRR28746893" "SRR28746894" "SRR28746895" \
# 	"SRR28746896" "SRR28746897" "SRR28746898" "SRR28746899" "SRR28746900")

# extract microbiome reads
echo "extract microbiome reads..."
for sample in ${samples[@]}; do
	echo "Processing the first file of $sample..."
	Rscript /home/yincy/git/github/SAHMI/functions/extract_microbiome_reads.r \
		--sample_name "$sample" \
		--fq $base_dir/${sample}_1.fq \
		--kraken_report $base_dir/${sample}.kraken.report.txt \
		--mpa_report $base_dir/${sample}.kraken.report.mpa.txt \
		--out_path "$base_dir/"
	mv $base_dir/${sample}.fa  $base_dir/${sample}_1.fa

	echo "processing the second file of $sample..."

	Rscript /home/yincy/git/github/SAHMI/functions/extract_microbiome_reads.r \
		--sample_name "$sample" \
		--fq $base_dir/${sample}_2.fq \
		--kraken_report $base_dir/${sample}.kraken.report.txt \
		--mpa_report $base_dir/${sample}.kraken.report.mpa.txt \
		--out_path "$base_dir/"
	mv $base_dir/${sample}.fa  $base_dir/${sample}_2.fa 

done

# extract microbiome output

echo "extract microbiome output..."
for sample in "${samples[@]}"; do
	echo ""
	echo "Processing $sample..."
	echo ""
	kraken_report="$base_dir/${sample}.kraken.report.txt"
	mpa_report="$base_dir/${sample}.kraken.report.mpa.txt"
	if [[ -f "$kraken_report" && -f "$mpa_report" ]]; then
		Rscript /home/yincy/git/github/SAHMI/functions/extract_microbiome_output.r \
			--sample_name "$sample" \
			--output_file $base_dir/${sample}.kraken.output.txt \
			--kraken_report "$kraken_report" \
			--mpa_report "$mpa_report" \
			--out_path "$base_dir/"
	else
		echo "Missing input files for $sample. Skipping."
	fi
	echo ""
	echo "Sample $sample has been processed..."
	echo ""
done
