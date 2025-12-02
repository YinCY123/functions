ld_clump_local <- function(dat, clump_kb, clump_r2, clump_p, bfile, plink_bin){
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn <- tempfile()
  write.table(data.frame(SNP=dat[["rsid"]], P=dat[["pval"]]), file=fn, row.names=F, col.names=T, quote=F)
  
  fun2 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --clump ", shQuote(fn, type=shell), 
    " --clump-p1 ", clump_p, 
    " --clump-r2 ", clump_r2, 
    " --clump-kb ", clump_kb, 
    " --threads 1 ",
    " --out ", shQuote(fn, type=shell)
  )
  system(fun2)
  
  # --- 关键修改从这里开始 ---
  clumped_file_path <- paste(fn, ".clumped", sep="")
  
  if (file.exists(clumped_file_path)) {
    # 如果文件存在，说明PLINK成功执行并产生了结果
    res <- read.table(clumped_file_path, header=T)
    unlink(paste(fn, "*", sep="")) # 清理临时文件
    
    # 检查读取的结果是否为空
    if (nrow(res) == 0) {
      message("PLINK clumping resulted in an empty .clumped file. Returning empty data frame.")
      return(data.frame())
    }
    
    # 原有的后续逻辑
    y <- subset(dat, !dat[["rsid"]] %in% res[["SNP"]])
    if(nrow(y) > 0)
    {
      message("Removing ", length(y[["rsid"]]), " of ", nrow(dat), " variants due to LD with other variants or absence from LD reference panel")
    }
    return(subset(dat, dat[["rsid"]] %in% res[["SNP"]]))
    
  } else {
    # 如果文件不存在，说明PLINK没有生成输出文件（通常是因为没有显著结果）
    message(paste("PLINK did not generate a .clumped file. This is likely because no SNPs passed the clumping criteria. Returning empty data frame."))
    unlink(paste(fn, "*", sep="")) # 同样要清理临时文件
    return(data.frame()) # 返回一个空数据框，而不是报错
  }
  # --- 修改结束 ---
}
