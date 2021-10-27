#!/usr/local/bin/Rscript

library(argparser)
library(tidyverse)
library(stringr)
library(ggsci)
library(data.table)
library(httr)
library(sqldf)
library(parallel)



'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

cosmic.run <- function(df, sample.column) {
  
  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  #X.1 <- parallel::mclapply(ls(X1), function(x) {
  X.1 <- lapply(ls(X1), function(x) {
    df.chr <- X1[[x]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP) ## need HGVSp_VEP 63:HGVSp_VEP, 55:SYMBOL_VEP 247:key
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    cosmic <- read.table(paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.final.more.minimal.test.",x,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
    colnames(cosmic) <- c("ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")
    
    cosmic.test <- sqldf("SELECT l.*, r.ID, r.var_key as var_key_cosmic, r.CosmicCount, r.ID, r.heme_cosmic_count, r.myeloid_cosmic_count
              FROM `df.chr` as l
              LEFT JOIN `cosmic` as r
              on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
    cosmic.test <- cosmic.test[,c(colnames(df.chr),"CosmicCount","heme_cosmic_count","myeloid_cosmic_count","ID")]
    cosmic.test$CosmicCount <- fillna(cosmic.test$CosmicCount,0)
    cosmic.test$heme_cosmic_count <- fillna(cosmic.test$heme_cosmic_count,0)
    cosmic.test$myeloid_cosmic_count <- fillna(cosmic.test$myeloid_cosmic_count,0)
    cosmic.test <- cosmic.test %>% 
      group_by(var_key,cosmic.test[[sample.column]]) %>%
      slice_max(order_by = heme_cosmic_count, n = 1, with_ties = F)
    cosmic.test <- data.frame(cosmic.test)
    
    print("done")
    df.chr$var_key_sample <- with(df.chr, paste(var_key,SAMPLE,sep=":"))
    cosmic.test$var_key_sample <- with(cosmic.test, paste(var_key,SAMPLE,sep=":"))
    row.names(cosmic.test) <- cosmic.test$var_key_sample
    cosmic.test <- cosmic.test[df.chr$var_key_sample,]
    
    if (all(cosmic.test[,c("var_key",sample.column)]==df.chr[,c("var_key",sample.column)])) {
      print(paste("good inner", x))
      return(cosmic.test)
    } else {
      print(paste("bad inner", x))
    }
  })
  proc.time() - ptm
  
  names(X.1) <- ls(X1)
  all.match <- all(sapply(ls(X1), function(x) {
    dim(X1[[x]])[1] == dim(X.1[[x]])[1]
  }))
  if (all.match) {
    cosmic.final2 <- bind_rows(X.1)
  }
  if (all(cosmic.final2[,c("var_key",sample.column)]==df[,c("var_key",sample.column)])) {
    print("good")
    MUTS <- cbind(df, cosmic.final2 %>% dplyr::select(myeloid_cosmic_count, ID))
    return(MUTS)
  } else {
    print("bad")
  }
}

parser <- arg_parser("Process NGS variant calls")
parser <- add_argument(parser, "--input", type="character", help="in file")
parser <- add_argument(parser, "--cols", type="character", help="in file")
parser <- add_argument(parser, "--output", type="character", help="out file")


args <- parse_args(parser)
# args$input="myeloid.input.txt"

final.test <- read.table(args$input, sep = "\t", header = T, quote = "", comment.char = "")
col_names <- read.table(args$cols, sep = "\t", comment.char = "")$V1
colnames(final.test) <- col_names
colnames(final.test)[247] <- "var_key"
###########annotate with COSMIC########
#source("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/cosmic.parallel.compute1.R")
## need to add a parallel
#MUTS <- cosmic.run(final.test[grep("^chr22:",final.test$var_key),], "SAMPLE")
MUTS <- cosmic.run(final.test, "SAMPLE")
#MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
#MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  
  
 
write.table(MUTS, args$output, row.names = F, sep = "\t", quote = F)
# ./RE_ANNOTATE_compute1.myeloid_count.R -i 6022994_23153_0_0.pon.pass.tsv -o 6022994_23153_0_0.pon.pass.myeloid.tsv -c 1218265_23153_0_0.columns.txt

