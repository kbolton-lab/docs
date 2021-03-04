# Mutational Signatures
getwd()
source('../mutational-descriptive-characteristics/plotSignature.R')  
library(BSgenome.Hsapiens.UCSC.hg19)

## Define Mutational Signature Patterns
subs.df <- filter(MSK.IMPACT_VCF, Variant_Type=='SNP' & nchar(MSK.IMPACT_VCF$Alt)==1  &  nchar(MSK.IMPACT_VCF$Ref)==1 )
subs.df$type <- paste0(subs.df$Ref, '>', subs.df$Alt)
nrow(subs.df)
table(subs.df$type)

sl <- seqlengths(Hsapiens)
subs.df$hasMargin <- subs.df$Start <  sl[paste0('chr',subs.df$Chrom)]
table( subs.df$hasMargin)
subs.df <- subset(subs.df,hasMargin==TRUE )

subs.df$triplets <- as.character(getSeq(Hsapiens, paste0('chr',subs.df$Chrom), start=subs.df$Start-1, end=subs.df$Start+1))

subs.df_clinical <- merge(MSK.IMPACT_FULL, subs.df, by="MRN", all.y=TRUE)

nrow(subs.df)

preBase <- rep(c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)),6)
refBase <- c(rep('C', 48 ), rep('T', 48))
altBase <- c(rep('A', 16 ), rep('G', 16), rep('T', 16),rep('A', 16 ), rep('C', 16), rep('G', 16))
postBase <- rep(c('A', 'C', 'G', 'T'), 96/4)
mut.order <- paste0(preBase, '[',refBase, '>', altBase, ']', postBase )

#Plot
subs.ann <- plotSignature(subs.df, mut.order, 'All Mutations', plot = TRUE)
