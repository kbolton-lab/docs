library(parallel)
library(dplyr)
library(sqldf)


'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

#'@returns dataframe with equal number of rows as input with binded counts columns
#'@param df dataframe: this is the MUTS dataframe; needs to have `HGVSp_VEP` column
#'@param sample.column string: column name for sample name/id
cosmic.run <- function(df, sample.column) {
  
  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  X.1 <- parallel::mclapply(ls(X1), function(x) {
  # X.1 <- lapply(ls(X1)[7:8], function(x) {
    df.chr <- X1[[x]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP)
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    cosmic <- read.table(paste0("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.final.more.minimal.test.",x,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
    colnames(cosmic) <- c("COSMIC_ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")
    
    cosmic.test <- sqldf("SELECT l.*, r.COSMIC_ID, r.var_key as var_key_cosmic, HGVSP, r.CosmicCount, 
              r.heme_cosmic_count, r.myeloid_cosmic_count
              FROM `df.chr` as l
              LEFT JOIN `cosmic` as r
              on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
    cosmic.test <- cosmic.test[,c(colnames(df.chr),"COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count")]
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
  
  if (length(X1) == length(X.1)) { 
    names(X.1) <- ls(X1) 
  } else {
      stop("WRONG")
    }
  
  all.match <- all(sapply(ls(X1), function(x) {
    dim(X1[[x]])[1] == dim(X.1[[x]])[1]
  }))
  if (all.match) {
    cosmic.final2 <- bind_rows(X.1)
  }
  if (all(cosmic.final2[,c("CHROM","POS","REF","ALT",sample.column)]==df[,c("CHROM","POS","REF","ALT",sample.column)])) {
    print("good")
    MUTS <- cbind(df, cosmic.final2 %>% dplyr::select(COSMIC_ID,CosmicCount,heme_cosmic_count,myeloid_cosmic_count))
    return(MUTS)
  } else {
    print("bad")
  }
}



MUTS <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/12/update/transplant_match/match/1218352_23153_0_0.final.tsv", 
                   sep = "\t", header = T, comment.char = "")
## remove original cosmic columns
dim(MUTS)
MUTS <- MUTS[,-c(grep("CosmicCount",colnames(MUTS)),
                 grep("heme_cosmic",colnames(MUTS)))
             ]
dim(MUTS)
MUTS$var_key <- with(MUTS, paste(CHROM,POS,REF,ALT,sep = ":"))
MUTS$HGVSp_VEP <- gsub("%3D","=",MUTS$HGVSp_VEP)
# df <- MUTS
MUTS <- cosmic.run(MUTS, "SAMPLE")



