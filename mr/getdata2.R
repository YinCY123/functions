#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿

####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com

setwd("G:\\eqtl_pqtl\\1")

library(data.table)
library(foreach)
source("ld_clump.R")
source("ld_matrix.R")
source("afl2.r")
source("api.R")
source("backwards.R")
source("query.R")
source("utils-pipe.R")
source("variants.R")
source("zzz.R")

genelist=read.table("drug_eQTLGen_id.txt",header = T,sep = "\t")

genelist1=as.vector(genelist$id)

foreach(i=genelist1, .errorhandling = "pass") %do%{

expo_rt=fread(file=paste0("eQTLGen/",i,".txt.gz"), header = T)

#head(expo_rt)

expo_rt=expo_rt[expo_rt$p<5e-8,]#5e-8

expo_rt2=expo_rt[,c("SNP","p")]
colnames(expo_rt2)=c("rsid", "pval")
#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿

####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com

clumdf <- ld_clump_local(dat = expo_rt2, clump_kb = 10000, clump_r2 = 0.1,clump_p=1,
                         bfile ="G:/eqtl_pqtl/1/data_maf0/data_maf0.01_rs_ref", 
                         plink_bin = "G:/eqtl_pqtl/1/plink_win64_20231018/plink.exe")

expo_rt3=expo_rt[which(expo_rt$SNP%in%clumdf$rsid),]

write.table(expo_rt3,file = paste0("eqtlclump/",i,"_expo_rt.txt"),row.names = F,sep = "\t",quote = F)
}

#数据与代码声明
#如果没有购买SCI狂人团队或者生信狂人团队的正版会员
#没有经过我们的同意，擅自使用我们整理好的数据与代码发文章
#如果被我们发现你的文章用了我们的数据与代码，我们将使用一切手段让你的文章撤稿

####关注微信公众号生信狂人团队
###遇到代码报错等不懂的问题可以添加微信scikuangren进行答疑
###作者邮箱：sxkrteam@shengxinkuangren.com


