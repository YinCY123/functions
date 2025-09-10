library(optparse)

option_list = list(
  make_option(c("--sample"), action="store", help = "sample name"),
  make_option(c("--fq1"), action="store", help = "path to fastq 1 file"),
  make_option(c("--fq2"), action="store", help = "path to fastq 2 file"),
  make_option(c("--out_path"), action="store", help = "output directory path"),
  make_option(c("--ncbi_blast_path"), action="store", help = "path to ncbi-blast"),
  make_option(c("--Kraken2Uniq_path"), action="store", help = "path to Kraken2 main 'kraken2' function"),
  make_option(c("--kraken_database_path"), action="store", help = "path to kraken database"),
  make_option(c("--kreport2mpa_path"), action="store", help = "path to kreport2mpa.py' function"),
  make_option(c("--paired"), action="store", default=T, help = "paired-end fastq files (T) or single-end (F)")
)
opt = parse_args(OptionParser(option_list = option_list))
opt$paired = as.logical(opt$paired)

library(stringr)

# Create output directory if it doesn't exist
dir.create(opt$out_path, showWarnings = FALSE, recursive = TRUE)

if(opt$paired == T){
  # run Kraken paired end
  kraken_cmd = paste0(
    'export PATH=$PATH:', opt$ncbi_blast_path, ' && ',
    opt$Kraken2Uniq_path,
    ' --db ', opt$kraken_database_path,
    ' --threads 24',
    ' --paired',
    ' --use-names',
    ' --report-minimizer-data',
    ' --classified-out ', file.path(opt$out_path, paste0(opt$sample, '#.fq')),
    ' --unclassified-out ', file.path(opt$out_path, paste0(opt$sample, '.unclassified#.fq')),
    ' --output ', file.path(opt$out_path, paste0(opt$sample, '.kraken.output.txt')),
    ' --report ', file.path(opt$out_path, paste0(opt$sample, '.kraken.report.txt')),
    ' ', opt$fq1,
    ' ', opt$fq2
  )
} else {
  # run Kraken unpaired
  kraken_cmd = paste0(
    'export PATH=$PATH:', opt$ncbi_blast_path, ' && ',
    opt$Kraken2Uniq_path,
    ' --db ', opt$kraken_database_path,
    ' --threads 24',
    ' --use-names',
    ' --report-minimizer-data',
    ' --classified-out ', file.path(opt$out_path, paste0(opt$sample, '#.fq')),
    ' --unclassified-out ', file.path(opt$out_path, paste0(opt$sample, '.unclassified#.fq')),
    ' --output ', file.path(opt$out_path, paste0(opt$sample, '.kraken.output.txt')),
    ' --report ', file.path(opt$out_path, paste0(opt$sample, '.kraken.report.txt')),
    ' ', opt$fq1
  )
}

# Execute Kraken command
system(kraken_cmd)

# Create MPA style report and standard Kraken report
report_cmd = paste0(
  'cut -f1-3,6-8 ',
  file.path(opt$out_path, paste0(opt$sample, '.kraken.report.txt')),
  ' > ',
  file.path(opt$out_path, paste0(opt$sample, '.kraken.report.std.txt')),
  ' && ',
  opt$kreport2mpa_path,
  ' -r ', file.path(opt$out_path, paste0(opt$sample, '.kraken.report.std.txt')),
  ' -o ', file.path(opt$out_path, paste0(opt$sample, '.kraken.report.mpa.txt')),
  ' --intermediate-ranks'
)

# Execute report generation command
system(report_cmd)
