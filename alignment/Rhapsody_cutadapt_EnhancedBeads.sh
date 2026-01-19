#!/bin/bash -i

CMDNAME=`basename $0`
if [ $# -ne 3 ]; then
 echo "Usage: sh $CMDNAME Read1.fastq.gz Read2.fastq.gz t<num_threads>"
 echo "read1 is cell barcode reads and read2 is mRNA reads" 
 exit 1
fi

file1=`basename $1 .fastq.gz`
file2=`basename $2 .fastq.gz`
threads=`echo $3`
samplename=`basename $2 | sed -e 's/_R2_001.fastq.gz//g'`


#trimming adapters by cutadapt v4, deal with phase-shift design (remove phase-shift bases) of BD Rhapsody Enhanced beads
#also check polyT existence

cutadapt --cores=$threads -Z --no-trim --report=minimal \
-g zero="XNNNNNNNNNGTG;min_overlap=12;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g one="XNNNNNNNNNNGTG;min_overlap=13;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g two="XNNNNNNNNNNNGTG;min_overlap=13;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g three="XNNNNNNNNNNNNGTG;min_overlap=14;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
--untrimmed-output untrim_${file1}.fastq.gz --untrimmed-paired-output untrim_${file2}.fastq.gz \
-o {name}_${file1}.fastq.gz -p {name}_${file2}.fastq.gz \
${file1}.fastq.gz ${file2}.fastq.gz 1> cutadapt_summary.log

#second round extraction. fuzzy match 1st spacers, and anchored by 2nd spacers
cutadapt --cores=$threads -Z --no-trim --report=minimal \
-g one="XNNNNNNNNNNGVGANNNNNNNNNGA;min_overlap=25;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g two="XNNNNNNNNNNNGVGANNNNNNNNNGA;min_overlap=26;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g three="XNNNNNNNNNNNNGVGANNNNNNNNNGA;min_overlap=27;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g one1="XNNNNNNNNNNHTGANNNNNNNNNGA;min_overlap=25;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g two1="XNNNNNNNNNNNHTGANNNNNNNNNGA;min_overlap=26;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g three1="XNNNNNNNNNNNNHTGANNNNNNNNNGA;min_overlap=27;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g one2="XNNNNNNNNNNGTHANNNNNNNNNGA;min_overlap=25;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g two2="XNNNNNNNNNNNGTHANNNNNNNNNGA;min_overlap=26;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
-g three2="XNNNNNNNNNNNNGTHANNNNNNNNNGA;min_overlap=27;max_error_rate=0.1...TTTTTTTTX;min_overlap=6;max_error_rate=0.2" \
--untrimmed-output untrim2_${file1}.fastq.gz --untrimmed-paired-output untrim2_${file2}.fastq.gz \
-o {name}_untrim_${file1}.fastq.gz -p {name}_untrim_${file2}.fastq.gz \
untrim_${file1}.fastq.gz untrim_${file2}.fastq.gz 1> /dev/null

cat one_untrim_${file1}.fastq.gz one1_untrim_${file1}.fastq.gz one2_untrim_${file1}.fastq.gz one_${file1}.fastq.gz > one3_${file1}.fastq.gz
cat one_untrim_${file2}.fastq.gz one1_untrim_${file2}.fastq.gz one2_untrim_${file2}.fastq.gz one_${file2}.fastq.gz > one3_${file2}.fastq.gz
cat two_untrim_${file1}.fastq.gz two1_untrim_${file1}.fastq.gz two2_untrim_${file1}.fastq.gz two_${file1}.fastq.gz > two3_${file1}.fastq.gz
cat two_untrim_${file2}.fastq.gz two1_untrim_${file2}.fastq.gz two2_untrim_${file2}.fastq.gz two_${file2}.fastq.gz > two3_${file2}.fastq.gz
cat three_untrim_${file1}.fastq.gz three1_untrim_${file1}.fastq.gz three2_untrim_${file1}.fastq.gz three_${file1}.fastq.gz > three3_${file1}.fastq.gz
cat three_untrim_${file2}.fastq.gz three1_untrim_${file2}.fastq.gz three2_untrim_${file2}.fastq.gz three_${file2}.fastq.gz > three3_${file2}.fastq.gz


cat untrim2_${file1}.fastq.gz zero_${file1}.fastq.gz > zero3_${file1}.fastq.gz
cat untrim2_${file2}.fastq.gz zero_${file2}.fastq.gz > zero3_${file2}.fastq.gz

rm *untrim*${file1}.fastq.gz
rm *untrim*${file2}.fastq.gz
rm one_${file1}.fastq.gz
rm one_${file2}.fastq.gz
rm two_${file1}.fastq.gz
rm two_${file2}.fastq.gz
rm three_${file1}.fastq.gz
rm three_${file2}.fastq.gz
rm zero_${file1}.fastq.gz
rm zero_${file2}.fastq.gz


cat one3_${file1}_trim.fastq.gz two3_${file1}_trim.fastq.gz three3_${file1}_trim.fastq.gz zero3_${file1}_trim.fastq.gz > ${file1}_trim.fastq.gz
cat one3_${file2}_trim.fastq.gz two3_${file2}_trim.fastq.gz three3_${file2}_trim.fastq.gz zero3_${file2}_trim.fastq.gz > ${file2}_trim.fastq.gz

rm one3_${file1}_trim.fastq.gz
rm one3_${file2}_trim.fastq.gz
rm two3_${file1}_trim.fastq.gz
rm two3_${file2}_trim.fastq.gz
rm three3_${file1}_trim.fastq.gz
rm three3_${file2}_trim.fastq.gz
rm zero3_${file1}_trim.fastq.gz
rm zero3_${file2}_trim.fastq.gz

exit 0
