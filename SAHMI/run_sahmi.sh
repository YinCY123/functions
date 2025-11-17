# fq_dir="/home/yincy/disk14/works/partime/2025/250418/GSE131882/unmapped"
# output_dir="/home/yincy/disk14/works/partime/2025/250418/GSE131882/SAHMI"
# samples=("SRR9141212_sorted" "SRR9141213_sorted" \
#   "SRR9141214_sorted" "SRR9141215_sorted"\
#   "SRR9141216_sorted" "SRR9141217_sorted" "SRR9141218_sorted" \
#   "SRR9141219_sorted" "SRR9141220_sorted" "SRR9141221_sorted")
# paired="TRUE"

run_shami(){

fq_dir="/home/yincy/disk14/works/partime/2025/250418/GSE131882/unmapped"
output_dir="/home/yincy/disk14/works/partime/2025/250418/GSE131882/SAHMI"
# samples=("SRR9141212" "SRR9141213" \
#   "SRR9141214" "SRR9141215"\
#   "SRR9141216" "SRR9141217" "SRR9141218" \
#   "SRR9141219" "SRR9141220" "SRR9141221")

samples=("SRR9141218" \
  "SRR9141219" "SRR9141220" "SRR9141221")


paired="TRUE"


# 01 run kraken
echo "run kraken..."
for sample in "${samples[@]}"; do
	echo "processing $sample..."
	Rscript functions/run_kraken.r \
		--sample "$sample" \
		--fq1 $fq_dir/${sample}_R1.fastq.gz \
		--fq2 $fq_dir/${sample}_R2.fastq.gz \
		--out_path $output_dir/ \
		--ncbi_blast_path /usr/lib/ncbi-blast+ \
		--Kraken2Uniq_path /home/yincy/bin/kraken2 \
		--kraken_database_path /home/yincy/ssd/kraken2/database/k2_minusb_20250402 \
		--paired T \
		--kreport2mpa_path functions/kreport2mpa.py
done

# 02 extrat microbiome reads
echo "extract microbiome reads..."
for sample in ${samples[@]}; do
	echo "Processing the first file of $sample..."
	Rscript functions/extract_microbiome_reads.r \
		--sample_name "$sample" \
		--fq $output_dir/${sample}_1.fq \
		--kraken_report $output_dir/${sample}.kraken.report.txt \
		--mpa_report $output_dir/${sample}.kraken.report.mpa.txt \
		--out_path $output_dir/
	mv $output_dir/${sample}.fa  $output_dir/${sample}_1.fa

    if [[ "$paired" == "TRUE" ]]; then
        echo "processing the second file of $sample..."

        Rscript functions/extract_microbiome_reads.r \
            --sample_name "$sample" \
            --fq $output_dir/${sample}_2.fq \
            --kraken_report $output_dir/${sample}.kraken.report.txt \
            --mpa_report $output_dir/${sample}.kraken.report.mpa.txt \
            --out_path $output_dir/
        mv $output_dir/${sample}.fa  $output_dir/${sample}_2.fa 
    fi
done

# extract microbiome output

echo "extract microbiome output..."
for sample in "${samples[@]}"; do
	echo ""
	echo "Processing $sample..."
	echo ""
	kraken_report="$output_dir/${sample}.kraken.report.txt"
	mpa_report="$output_dir/${sample}.kraken.report.mpa.txt"
	if [[ -f "$kraken_report" && -f "$mpa_report" ]]; then
		Rscript functions/extract_microbiome_output.r \
			--sample_name "$sample" \
			--output_file $output_dir/${sample}.kraken.output.txt \
			--kraken_report "$kraken_report" \
			--mpa_report "$mpa_report" \
			--out_path $output_dir/
	else
		echo "Missing input files for $sample. Skipping."
	fi
	echo ""
	echo "Sample $sample has been processed..."
	echo ""
done

# 03 single cell kmer analysis
for sample in "${samples[@]}"; do
	if [[ "$paired" == TRUE ]]; then
	if [[ "$1" == "--help" || "$1" == "-h" ]]; then
		echo "Usage: run_shami [options]"
		echo ""
		echo "Required parameters:"
		echo "  --fq_dir, -q          Path to input fastq directory"
		echo "	--output_dir, -o      Path to output directory"
		echo "	--samples, -s         Path to file containing sample names (one per line)"
		echo "	--paired, -p          Whether the data is paired-end (TRUE/FALSE)"
		echo ""
		return 0
	fi

	local fq_dir
	local output_dir
	local paired
	local sample_file
	declare -a samples=()


	while [[ #$ -gt 1 ]]; do
		case $1 in
			--fq_dir|-q)
				fq_dir=$2
				shift 2
				;;
			--output_dir|-o)
				output_dir=$2
				shift 2
				;;
			--samples|-s)
				sample_file=$2 # append the sample to the array
				shift 2
				;;
			--paired|-p)
				paired=$2
				shift 2
				;;
			*)
				echo "Unknown option: $1"
				return 1
				;;
		esac
	done

	# load sample from file
	if [[ -n $sample_file ]]; then
		if [[ -f $sample_file ]]; then
			mapfile -t samples < $sample_file
		else 
			echo "Error: sample file $sample_file not found."
			return 1
		fi
	else 
		echo "Error: --samples option is required."
		return 1
	fi

	# 01 run kraken
	echo "run kraken..."
	for sample in "${samples[@]}"; do
		echo "processing $sample..."
		Rscript functions/run_kraken.r \
			--sample "$sample" \
			--fq1 $fq_dir/${sample}_R1.fastq.gz \
			--fq2 $fq_dir/${sample}_R2.fastq.gz \
			--out_path $output_dir/ \
			--ncbi_blast_path /usr/lib/ncbi-blast+ \
			--Kraken2Uniq_path /home/yincy/bin/kraken2 \
			--kraken_database_path /home/yincy/ssd/kraken2/database/k2_minusb_20250402 \
			--paired T \
			--kreport2mpa_path functions/kreport2mpa.py
	done

	# 02 extrat microbiome reads
	echo "extract microbiome reads..."
	for sample in ${samples[@]}; do
		echo "Processing the first file of $sample..."
		Rscript functions/extract_microbiome_reads.r \
			--sample_name "$sample" \
			--fq $output_dir/${sample}_1.fq \
			--kraken_report $output_dir/${sample}.kraken.report.txt \
			--mpa_report $output_dir/${sample}.kraken.report.mpa.txt \
			--out_path $output_dir/
		mv $output_dir/${sample}.fa  $output_dir/${sample}_1.fa

		if [[ "$paired" == "TRUE" ]]; then
			echo "processing the second file of $sample..."

			Rscript functions/extract_microbiome_reads.r \
				--sample_name "$sample" \
				--fq $output_dir/${sample}_2.fq \
				--kraken_report $output_dir/${sample}.kraken.report.txt \
				--mpa_report $output_dir/${sample}.kraken.report.mpa.txt \
				--out_path $output_dir/
			mv $output_dir/${sample}.fa  $output_dir/${sample}_2.fa 
		fi
	done

	# extract microbiome output

	echo "extract microbiome output..."
	for sample in "${samples[@]}"; do
		echo ""
		echo "Processing $sample..."
		echo ""
		kraken_report="$output_dir/${sample}.kraken.report.txt"
		mpa_report="$output_dir/${sample}.kraken.report.mpa.txt"
		if [[ -f "$kraken_report" && -f "$mpa_report" ]]; then
			Rscript functions/extract_microbiome_output.r \
				--sample_name "$sample" \
				--output_file $output_dir/${sample}.kraken.output.txt \
				--kraken_report "$kraken_report" \
				--mpa_report "$mpa_report" \
				--out_path $output_dir/
		else
			echo "Missing input files for $sample. Skipping."
		fi
		echo ""
		echo "Sample $sample has been processed..."
		echo ""
	done

	# 03 single cell kmer analysis
	for sample in "${samples[@]}"; do
		if [[ "$paired" == TRUE ]]; then
			echo ""
			echo "processing paired sample: $sample"
			# echo ""
			Rscript functions/sckmer.r \
				--sample_name $sample \
				--fa1 ${output_dir}/${sample}_1.fa \
				--fa2 ${output_dir}/${sample}_2.fa \
				--microbiome_output_file ${output_dir}/${sample}.microbiome.output.txt \
				--kraken_report ${output_dir}/${sample}.kraken.report.txt \
				--mpa_report ${output_dir}/${sample}.kraken.report.mpa.txt \
				--out_path $output_dir/ \
				--cb_len 16 \
				--umi_len 12 \
				--ranks 'G'
		else
			echo "processing single end sample: $sample"
			Rscript functions/sckmer_unpaired.r \
				--sample_name $sample \
				--fa1 ${output_dir}/${sample}_1.fa \
				--microbiome_output_file ${output_dir}/${sample}.microbiome.output.txt \
				--kraken_report ${output_dir}/${sample}.kraken.report.txt \
				--mpa_report ${output_dir}/${sample}.kraken.report.mpa.txt \
				--out_path $output_dir/ \
				--cb_len 16 \
				--umi_len 12 \
				--ranks 'G'
		fi
	done

	# 04 taxa counts
	for sample in "${samples[@]}"; do
		echo "Processing $sample..."
		Rscript functions/taxa_counts.r \
			--sample_name $sample \
			--fa1 $output_dir/${sample}_1.fa \
			--fa2 $output_dir/${sample}_2.fa \
			--taxa "taxa.tsv" \
			--kraken_report "$output_dir"/${sample}.kraken.report.txt \
			--mpa_report "$output_dir"/${sample}.kraken.report.mpa.txt \
			--out_path $output_dir/ \
			--cb_len 16 \
			--umi_len 12
		done
}
