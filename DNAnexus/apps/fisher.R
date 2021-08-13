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
stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_AltCounts\t%INFO/PON_RefDepth\t[%AD]", call.=FALSE)
}

df$fisher.exact.pval <- apply(df, 1, function(x) {
x <- as.numeric(x[-c(1,2,3,4)])
if ((x[1]+x[2])==0) {
    return(0)
} else {
    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
}
})
write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")