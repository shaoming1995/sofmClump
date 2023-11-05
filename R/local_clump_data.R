#' @title 提供本地聚类分析的R包，但是需要联系顶刊研习社获取参考面板数据或者微信SFM19950928
#' @param keyssh 输入密钥
#' @param temp_dat 需要聚类分析的变量名
#' @param filepath 参考面板数据（.bed/.fam/.bim）的文件位置
#' @param pop 人种选择，仅支持EUR/EAS
#' @param clump_kb 聚类的kb
#' @param clump_r2 聚类的r2
#' @export

local_clump_data<-function(keyssh,temp_dat,filepath,pop,clump_kb,clump_r2){
  if (!is.null(keyssh)) {
    warning("该密钥可以联系顶刊研习社获取或者抖音ID医小研~")
  }
  if(!require("tidyr",quietly=T))
      install.packages("tidyr")
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  library(plinkbinr)
  library(tidyr)
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  library(key)
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  if(Sys.info()["nodename"]==keyssh){
  temp_dat$id <- temp_dat$id.exposure
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  temp_dat$rsid <- temp_dat$SNP
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  temp_dat$pval <- temp_dat$pval.exposure
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  filepath1<-paste0(filepath,"/",pop)
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  ld_sofm<-function (dat = NULL, clump_kb = 10000, clump_r2 = 0.001, clump_p = 0.99,
            pop = "EUR", access_token = NULL, bfile = NULL, plink_bin = NULL)
  {
    stopifnot("rsid" %in% names(dat))
    stopifnot(is.data.frame(dat))
    if (is.null(bfile)) {
      message("Please look at vignettes for options on running this locally if you need to run many instances of this command.")
    }
    if (!"pval" %in% names(dat)) {
      if ("p" %in% names(dat)) {
        warning("No 'pval' column found in dat object. Using 'p' column.")
        dat[["pval"]] <- dat[["p"]]
      }
      else {
        warning("No 'pval' column found in dat object. Setting p-values for all SNPs to clump_p parameter.")
        dat[["pval"]] <- clump_p
      }
    }
    if (!"id" %in% names(dat)) {
      dat$id <- random_string(1)
    }
    if (is.null(bfile)) {
      access_token = check_access_token()
    }
    ids <- unique(dat[["id"]])
    res <- list()
    for (i in 1:length(ids)) {
      x <- subset(dat, dat[["id"]] == ids[i])
      if (nrow(x) == 1) {
        message("Only one SNP for ", ids[i])
        res[[i]] <- x
      }
      else {
        message("Clumping ", ids[i], ", ", nrow(x), " variants, using ",
                pop, " population reference")
        if (is.null(bfile)) {
          res[[i]] <- ld_clump_api(x, clump_kb = clump_kb,
                                   clump_r2 = clump_r2, clump_p = clump_p, pop = pop,
                                   access_token = access_token)
        }
        else {
          res[[i]] <- ld_clump_local(x, clump_kb = clump_kb,
                                     clump_r2 = clump_r2, clump_p = clump_p, bfile = bfile,
                                     plink_bin = plink_bin)
        }
      }
    }
    res <- dplyr::bind_rows(res)
    return(res)
  }
  if(!require("TwoSampleMR",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  temp_dat$rsid <- temp_dat$SNP
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  temp_dat$pval <- temp_dat$pval.exposure
  if(!require("plinkbinr",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  filepath1<-paste0(filepath,"/",pop)

  if(!require("TwoSampleMR",quietly=T))
    devtools::install_github("explodecomputer/plinkbinr")
  temp_dat <- ld_sofm(temp_dat,
                       plink_bin = get_plink_exe(),
                       bfile = filepath1,
                       clump_kb = clump_kb, clump_r2 = clump_r2,pop=pop)
  temp_dat$rsid <- NULL
  temp_dat$pval <- NULL
  return(temp_dat)
  cat("已经完成聚类分析")}
  else{warning("keyssh不正确,请联系管理员微信SFM19950928获取密钥")}
}























