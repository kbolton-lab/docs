#!/usr/bin/Rscript
library(parallel)
# library(dplyr)
library(sqldf)


fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}
i <- commandArgs(trailingOnly = TRUE)
x <- paste0("chr", i)
data <- read.table("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/transplant.ch_pd2.for.cosmic.for.join.tsv",
                   header = T, sep = "\t", comment.char = "", quote="")

X1 <- split(data, data$CHROM)
ptm2 <- proc.time()
df <- X1[[x]]
df$HGVSp_VEP <- gsub("\\.[0-9]{1,2}","",df$HGVSp_VEP)
df$HGVSp_VEP <- gsub("%3D","=",df$HGVSp_VEP)
cosmic <- read.table(paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.final.minimal.",x,".tsv"),
                     header = T, sep = "\t", comment.char = "", quote="")
cosmic$HGVSP <- gsub("\\.[0-9]{1,2}","",cosmic$HGVSP)
colnames(cosmic) <- c("ID","tumor_id","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count")
cosmic.test <- sqldf("SELECT l.*, r.ID, r.var_key as var_key_cosmic, HGVSP, r.CosmicCount as CosmicCount2, 
          r.heme_cosmic_count as heme_cosmic_count2, r.myeloid_cosmic_count
       FROM `df` as l
       LEFT JOIN `cosmic` as r
       on l.var_key = r.var_key OR l.HGVSp_VEP = r.HGVSP")
cosmic.test <- cosmic.test[,c("eid","var_key","HGVSp_VEP","var_key_cosmic","ID","HGVSP","CosmicCount2","heme_cosmic_count2","myeloid_cosmic_count")]
cosmic.test$CosmicCount2 <- fillna(cosmic.test$CosmicCount2,0)
cosmic.test$heme_cosmic_count2 <- fillna(cosmic.test$heme_cosmic_count2,0)
cosmic.test$myeloid_cosmic_count <- fillna(cosmic.test$myeloid_cosmic_count,0)
cosmic.test <- cosmic.test %>% 
  group_by(var_key,eid) %>%
  slice_max(order_by = heme_cosmic_count2, n = 1, with_ties = F)
print(paste(x, proc.time() - ptm2, sep = ": "))
final.cosmic.annotated <- data.frame(cosmic.test)

row.names(final.cosmic.annotated) <- with(final.cosmic.annotated, paste(var_key,eid, sep = ":"))
final <- final.cosmic.annotated[with(df, paste(var_key,eid, sep=":")),]

write.table(final,paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/cosmic_count.",x,".tsv"),
            col.names = T, row.names = F, sep = "\t", quote = F)

