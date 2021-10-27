#!/usr/local/bin/Rscript
library(parallel)
library(dplyr)
library(sqldf)


fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

cosmic.run <- function(df, sample.column) {

  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  X.1 <- parallel::mclapply(ls(X1), function(x) {
    df.chr <- X1[[x]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP)
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    cosmic <- read.table(paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.final.more.minimal.test.",x,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
    colnames(cosmic) <- c("ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")
    
    cosmic.test <- sqldf("SELECT l.*, r.ID, r.var_key as var_key_cosmic, HGVSP, r.CosmicCount, 
              r.heme_cosmic_count, r.myeloid_cosmic_count
              FROM `df.chr` as l
              LEFT JOIN `cosmic` as r
              on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
    cosmic.test <- cosmic.test[,c(colnames(df.chr),"CosmicCount","heme_cosmic_count","myeloid_cosmic_count")]
    cosmic.test$CosmicCount <- fillna(cosmic.test$CosmicCount,0)
    cosmic.test$heme_cosmic_count <- fillna(cosmic.test$heme_cosmic_count,0)
    cosmic.test$myeloid_cosmic_count <- fillna(cosmic.test$myeloid_cosmic_count,0)
    cosmic.test <- cosmic.test %>% 
      group_by(var_key,cosmic.test[[sample.column]]) %>%
      slice_max(order_by = heme_cosmic_count, n = 1, with_ties = F)
    cosmic.test <- data.frame(cosmic.test)
  
    row.names(cosmic.test) <- cosmic.test$var_key
    cosmic.test <- cosmic.test[df.chr$var_key,]
    
    if (all(cosmic.test[,c("CHROM","POS","REF","ALT",sample.column)]==df.chr[,c("CHROM","POS","REF","ALT",sample.column)])) {
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
  if (all(cosmic.final2[,c("CHROM","POS","REF","ALT",sample.column)]==MUTS[,c("CHROM","POS","REF","ALT",sample.column)])) {
    print("good")
    MUTS <- cbind(MUTS, cosmic.final %>% dplyr::select(CosmicCount,heme_cosmic_count,myeloid_cosmic_count))
    return(MUTS)
  } else {
    print("bad")
  }
}




