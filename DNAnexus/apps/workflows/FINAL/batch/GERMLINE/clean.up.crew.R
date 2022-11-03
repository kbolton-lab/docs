library(chron)
library(dplyr)
###################################################
# alignment.qc.txt
# dx ls project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" > test.txt
# cat test.txt | cut -d_ -f1 | sort | uniq -c | awk '$3=($1/4)-1' | awk '$3>=1' | awk '$4=$3*4' | sort -k1,1rn | awk '{print $4}' | SUM
###################################################
germline.alignment.qc <- read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/alignment.qc.txt", header = T,
                                    sep = ",")
germline.alignment.qc <- germline.alignment.qc[,1:5]
germline.alignment.qc$eid <- gsub("\\..*","",germline.alignment.qc$Name)
germline.alignment.qc$file.name <- gsub(".*bqsr.{1}","",germline.alignment.qc$Name)
# gsub("\\..*","",germline.alignment.qc$Name[1])
# gsub(".*bqsr.{1}","",germline.alignment.qc$Name)

dtparts = as.data.frame(t(strsplit(germline.alignment.qc$Last.modified,' ')))
dtparts <- t(as.data.frame(unname(strsplit(germline.alignment.qc$Last.modified,' '))))
row.names(dtparts) <- NULL
thetimes = chron(dates=dtparts[,1], times=dtparts[,2], format=c('y-m-d','h:m:s'))
germline.alignment.qc$date.time <- thetimes
chron(germline.alignment.qc$Last.modified[1:5], format = c(dates = "y-m-d", times = "h:m:s"))

test <- germline.alignment.qc %>%
  group_by(eid, file.name) %>%
  arrange(desc(date.time)) %>%
  top_n(1, date.time)

file.ids.to.delete <- germline.alignment.qc$ID[!(germline.alignment.qc$ID %in% test$ID)]
length(file.ids.to.delete) + nrow(test) == dim(germline.alignment.qc)[1]

write.table(file.ids.to.delete, "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/alignment.qc.delete.txt",
            col.names = F, row.names = F, quote = F)


test <- read.table("/Users/brian/test/4390335_23153_0_0.final.tsv", header = T,
                                            sep = "\t", comment.char = "")



###################################################
# FilterMutectCalls.txt
###################################################

germline.FilterMutectCalls <- read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/FilterMutectCalls.txt", header = F,
                                    sep = ",")
germline.FilterMutectCalls <- germline.FilterMutectCalls[,1:5]
colnames(germline.FilterMutectCalls) <- colnames(germline.alignment.qc)[1:5]
germline.FilterMutectCalls$eid <- gsub("\\..*","",germline.FilterMutectCalls$Name)
germline.FilterMutectCalls$file.name <- gsub(".*bqsr.{1}","",germline.FilterMutectCalls$Name)

dtparts = as.data.frame(t(strsplit(germline.FilterMutectCalls$Last.modified,' ')))
dtparts <- t(as.data.frame(unname(strsplit(germline.FilterMutectCalls$Last.modified,' '))))
row.names(dtparts) <- NULL
thetimes = chron(dates=dtparts[,1], times=dtparts[,2], format=c('y-m-d','h:m:s'))
germline.FilterMutectCalls$date.time <- thetimes


test2 <- germline.FilterMutectCalls %>%
  group_by(eid, Name) %>%
  arrange(desc(date.time)) %>%
  top_n(1, date.time)

file.ids.to.delete <- germline.FilterMutectCalls$ID[!(germline.FilterMutectCalls$ID %in% test2$ID)]
length(file.ids.to.delete) + nrow(test2) == dim(germline.FilterMutectCalls)[1]

write.table(file.ids.to.delete, "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/FilterMutectCalls.delete.txt",
            col.names = F, row.names = F, quote = F)


# test <- read.table("/Users/brian/test/4390335_23153_0_0.final.tsv", header = T,
#                    sep = "\t", comment.char = "")


  