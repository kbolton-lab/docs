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
  stop("Must supply file with 7 columns to split to 8: %CHROM\t%POS\t%REF\t%ALTt%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]", call.=FALSE)
}

df$fisher.exact.pval <- apply(df, 1, function(x) {
  x <- as.numeric(x[-c(1,2,3,4)])
  if (x[2]==0 & x[1]!=0) {
    return(0)
  } else if ((x[1]==0 & x[2]!=0) | (x[3]==0 & x[4]!=0)) {
    return(1)
  } else {
    return(fisher.test(matrix(c(x[1], x[2], x[3], x[4]), ncol=2))$p.value)
  }
})


# df$X1 <- as.integer(df$X1)
# df$X2 <- as.integer(df$X2)
# Sys.setenv(JULIA_NUM_THREADS = "6")
# library(JuliaConnectoR)
# juliaEval("Threads.nthreads()")
# juliaEval("using HypothesisTests")

# fisherstestjulia <- juliaEval('function fisherstestjulia(a,b,c,d)
#   if typeof(d[1]) == String
#     d=parse.(Int,d)
#   end
#   ls = zeros(length(a))
#   for i in 1:length(ls)
#     # if (a[i]+b[i]==0)
#     if (a[i]==0 && b[i]==0)
#       continue # this is most likely because the reference allele is actually the minor allele
#     else if (a[i]==0 && c[i]==0) # all normals and sample are homozygous for alternate allele, i.e. the Reference is the minor allele
#       ls[i] = 1
#     else
#       try
#         ls[i] = pvalue( FisherExactTest(a[i], b[i], c[i], d[i] ), method = :minlike )
#       catch
#         ls[i] = pvalue( FisherExactTest(a[i], b[i], c[i], d[i]), method = :minlike )
#       end
#     end
#   end
#   return ls
# end')

# df$fisher.exact.pval=fisherstestjulia(df[,5],df[,6],df[,7],df[,8])
write.table(df, file=args[2], row.names = F, quote = F, col.names = F, sep = "\t")