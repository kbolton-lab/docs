#!/usr/bin/env Rscript 

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  stop("Must supply (output file).n", call.=FALSE)
}

df = read.table(args[1], header=F)
#https://statisticsglobe.com/split-data-frame-variable-into-multiple-columns-in-r
df <- cbind(df, data.frame(do.call("rbind", strsplit(as.character(df$V7), ",", fixed = TRUE))))[,-7]
if (length(colnames(df)) != 8) {
  stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALTt%INFO/PON_RefDepth\t%INFO/PON_AltDepth\\t[%AD]", call.=FALSE)
}

df$fisher.exact.pval <- apply(df, 1, function(x) {
  x <- as.numeric(x[-c(1,2,3,4)])
  if (x[2]==0 & x[1]!=0) {
    return(0)
  } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
<<<<<<< HEAD
    return(1) 
  } else if (x[2]/(x[1]+x[2]) >= x[4]/(x[3]+x[4])) {
=======
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
    return(1)
  } else {
    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
  }
})


write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")