#' @title 提供本地聚类分析的R包，但是需要联系顶刊研习社获取参考面板数据或者微信SFM19950928
#' @param keyssh 输入密钥
#' @param temp_dat 需要聚类分析的变量名
#' @param filepos 参考面板数据（.bed/.fam/.bim）的文件位置
#' @param pop 人种选择，仅支持EUR/EAS
#' @param clump_kb 聚类的kb
#' @param clump_r2 聚类的r2
#' @export

local_clump_data<-function(keyssh,temp_dat,filepos,pop,clump_kb,clump_r2){
  if(!require("tidyr",quietly=T))
      install.packages("tidyr")
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  if(!require("TwoSampleMR",quietly=T))
    devtools::install_github("MRCIEU/TwoSampleMR")
  if(!require("ieugwasr",quietly=T))
  devtools::install_github("mrcieu/ieugwasr")
  if(!require("sofmkey",quietly=T))
    devtools::install_github("shaoming1995/sofmkey")
  library(plinkbinr)
  library(TwoSampleMR)
  library(ieugwasr)
  library(tidyr)
  library(key)
  K<-key(keyssh)
  if(K==9544){
 remove.packages("key")
  # 补充ld_clump需要的三列
  temp_dat$id <- temp_dat$id.exposure # id列用于区分是单个或多个来源的gwas数据，clump按照单个gwas数据进行
  temp_dat$rsid <- temp_dat$SNP
  temp_dat$pval <- temp_dat$pval.exposure
  # 执行ld_clump
  # bfile指定参考文件的路径

  temp_dat <- ld_clump(temp_dat,
                       # get_plink_exe()
                       plink_bin = get_plink_exe(),
                       bfile = filepos,
                       clump_kb = 500, clump_r2 = 0.01,pop=pop)
  # 运行结果：
  # Removing 69 of 92 variants due to LD with other variants or absence from LD reference panel
  # 移除重复的列：rsid，pval
  temp_dat$rsid <- NULL
  temp_dat$pval <- NULL
  print("已经完成聚类分析")}else{cat("keyssh不正确,请联系管理员微信SFM19950928获取密钥")
    }
}
