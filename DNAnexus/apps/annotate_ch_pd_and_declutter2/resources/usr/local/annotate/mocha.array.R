chrOrder <-c((1:22),"X","Y","M")
mocha.ukbb.array.map <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/array.sample.eid.map",
                                   header = F, sep = "\t", comment.char = "", quote="")

colnames(mocha.ukbb.array.map) <- c("sample_id", "EID")
results <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/mocha_results/ukbb.calls.tsv",
                      header = T, sep = "\t", comment.char = "", quote="")
filtered_results <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/mocha_results/ukbb.calls.filtered.tsv",
                             header = T, sep = "\t", comment.char = "", quote="")
germline.eids <- read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/cases_1111.uniq.eids",
                            header = F, sep = "\t", comment.char = "", quote="")

write.table(filtered_results, "/Volumes/bolton/Active/projects/mocha/UKBB/arrays/mocha_results/bolton.ukbb.calls.filtered.with.eid.tsv",
            row.names = F, quote = F, sep = "\t")
bolton.filtered_results <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/mocha_results/bolton.ukbb.calls.filtered.with.eid.tsv",
                                      header = F, sep = "\t", comment.char = "", quote="")
# results <- results %>% left_join(mocha.ukbb.array.map, by="sample_id")

filtered_results <- filtered_results %>% left_join(mocha.ukbb.array.map, by="sample_id")
filtered_results <- filtered_results[,c(1,23,2:22)]
colnames(filtered_results)[5:7] <- c("START_MB", "END_MB", "SIZE_MB")
filtered_results$START_MB <- round(filtered_results$START_MB/1000000, 2)
filtered_results$END_MB <- round(filtered_results$END_MB/1000000, 2)
filtered_results$SIZE_MB <- round(filtered_results$SIZE_MB/1000000, 2)
colnames(filtered_results)[4] <- "CHR"
# colnames(filtered_results)[5] <- "START_GRCh38"
# colnames(filtered_results)[6] <- "END_GRCh38"
filtered_results <- filtered_results[order(filtered_results$EID, filtered_results$CHR, 
                                           filtered_results$START_MB, filtered_results$type),]
# filtered_results$CHR <- as.numeric(stringr::str_replace(filtered_results$CHR, "chr", ""))
filtered_results$CHR <- stringr::str_replace(filtered_results$CHR, "chr", "")
filtered_results$CHR <- factor(filtered_results$CHR,levels=chrOrder, ordered=TRUE)
filtered_results <- filtered_results[do.call(order, filtered_results[,c("EID", "CHR", "START_MB", "type")]), ]
filtered_results <- filtered_results[filtered_results$EID > 0,]
filtered_results.g <- filtered_results %>%
  group_by(EID) %>%
  # mutate(num=n()) %>%
  summarise(num=n()) %>%
  slice_max(num, n=3, with_ties = T)
bolton.higest.cnv.count <- filtered_results[filtered_results$EID=="1025117",]
# bolton.higest.cnv.count$CHR <- as.numeric(stringr::str_replace(bolton.higest.cnv.count$CHR, "chr", ""))
# bolton.higest.cnv.count <- bolton.higest.cnv.count[order(bolton.higest.cnv.count$CHR, bolton.higest.cnv.count$START_GRCh38, bolton.higest.cnv.count$type),]
rownames(bolton.higest.cnv.count) <- 1:nrow(bolton.higest.cnv.count)
# bolton.higest.cnv.count[,c("EID", "CHR", "START_GRCh38", "END_GRCh38", "type")]
bolton.higest.cnv.count[,c("EID", "CHR", "START_MB", "END_MB", "SIZE_MB", "type")]
# bolton.higest.cnv.count2 <- bolton.higest.cnv.count
# colnames(bolton.higest.cnv.count2)[5:6] <- colnames(Poh.return.3094)[3:4]
# colnames(bolton.higest.cnv.count2)[7] <- "SIZE_MB"
# bolton.higest.cnv.count2$START_MB <- round(bolton.higest.cnv.count2$START_MB/1000000, 2)
# bolton.higest.cnv.count2$END_MB <- round(bolton.higest.cnv.count2$END_MB/1000000, 2)
# bolton.higest.cnv.count2$SIZE_MB <- round(bolton.higest.cnv.count2$SIZE_MB/1000000, 2)
# 5500172 = 19 
bolton.higest.cnv.count2 <- filtered_results[filtered_results$EID=="5500172",]
rownames(bolton.higest.cnv.count2) <- 1:nrow(bolton.higest.cnv.count2)
bolton.higest.cnv.count2[,c("EID", "CHR", "START_MB", "END_MB", "SIZE_MB", "type")]

Poh.return.3094 <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/mocha_results/File_for_Retman/CNV_data_19808.2.txt",
                          header = T, sep = "\t", comment.char = "", quote="")
bridge19808 <- fread("~/Downloads/ukb55288bridge19808.txt")
colnames(bridge19808) <- c("eid", "ID")
ourData19808 <-
  inner_join(Poh.return.3094, bridge19808, by = "ID") %>% dplyr::select(-ID)
## from Application 55288
mocha.ukbb.array.EIDs <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/ukb22418_c1_b0_v2.fam",
                                    header = F, comment.char = "", quote="")
Poh.sample.EIDs <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/arrays/mocha_results/File_for_Retman/ukb4777.samples_analyzed.txt",
                              header = F, sep = "\t", comment.char = "", quote="")
sum(mocha.ukbb.array.EIDs$V1 %in% Poh.sample.EIDs$V1)
# Poh.return.3094 <- Poh.return.3094[order(Poh.return.3094$ID, Poh.return.3094$CHR, 
#                                          Poh.return.3094$START_MB, Poh.return.3094$COPY_CHANGE),]
Poh.return.3094$CHR <- factor(Poh.return.3094$CHR,levels=chrOrder, ordered=TRUE)
Poh.return.3094 <- Poh.return.3094[do.call(order, Poh.return.3094[,c("ID", "CHR", "START_MB", "COPY_CHANGE")]), ]
colnames(Poh.return.3094)[1] <- "EID"

Poh.return.3094.g <- Poh.return.3094 %>%
  group_by(EID) %>%
  # mutate(num=n()) %>%
  summarise(num=n()) %>%
  slice_max(num, n=5, with_ties = T)


poh.higest.cnv.count <- Poh.return.3094[Poh.return.3094$EID=="4851293",]
# poh.higest.cnv.count[do.call(order, poh.higest.cnv.count[,c("CHR", "START_MB")]), c("EID", "CHR", "START_MB", "END_MB", "SIZE_MB","COPY_CHANGE")]
rownames(poh.higest.cnv.count) <- 1:nrow(poh.higest.cnv.count)
poh.higest.cnv.count[,c("EID", "CHR", "START_MB", "END_MB", "SIZE_MB","COPY_CHANGE")]
# poh.higest.cnv.count$CHR <- as.numeric(poh.higest.cnv.count$CHR)


filtered_results.g
bolton.higest.cnv.count2 <- filtered_results[filtered_results$EID=="1929950",]
rownames(bolton.higest.cnv.count2) <- 1:nrow(bolton.higest.cnv.count2)
bolton.higest.cnv.count2[,c("EID", "CHR", "START_MB", "END_MB", "SIZE_MB", "type")]

Poh.return.3094.g
poh.higest.cnv.count2 <- Poh.return.3094[Poh.return.3094$EID=="4913544",]
rownames(poh.higest.cnv.count2) <- 1:nrow(poh.higest.cnv.count2)
poh.higest.cnv.count2[,c("EID", "CHR", "START_MB", "END_MB", "SIZE_MB","COPY_CHANGE")]

bolton.males <- filtered_results[filtered_results$computed_gender=="M",]
bolton.females <- filtered_results[filtered_results$computed_gender=="F",]
bolton.males$EID

nrow(poh.Y.males)
length(unique(poh.Y.males$EID))

ukbb <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbbExome_09232021.rds")
ukbb.males <- ukbb[ukbb$sex == "Male",]
ukbb.females <- ukbb[ukbb$sex == "Female",]

poh.Y.males <- Poh.return.3094[Poh.return.3094$CHR=="Y",]

sum(unique(ukbb.males$eid) %in% unique(poh.Y.males$EID))
sum(unique(ukbb.females$eid) %in% unique(poh.Y.males$EID))

intersect(bolton.males$EID, poh.Y.males$EID)
"1000047" %in% filtered_results$EID
intersect(filtered_results$EID, Poh.return.3094$EID)

sum(Poh.return.3094$EID %in% filtered_results$EID)
sort(filtered_results$EID)[1:20]
par(mfrow=c(2,1))


########################################################################################################################
########################################################################################################################


germline <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/exome/germline/HaplotypeCaller2/stats.txt",
                       header = F, sep = "\t", comment.char = "", quote="")
colnames(germline) <-c("QD","FS","SOR","MQ","MQRankSum","ReadPosRankSum")
germline$MQRankSum <- as.numeric(germline$MQRankSum)
germline$ReadPosRankSum
for (colname in colnames(germline)) {
  print(colname)
  plot(density(germline[[colname]]), main=colname)
}



CCOC <- readRDS("/Users/brian/Bolton/CCOC/CCOC.RNAseq.normalized.rds")
length(rownames(CCOC))

head(rownames(CCOC))

test2 <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/lung/update/5104995_23153_0_0.final.tsv",
                    header = T, sep = "\t", comment.char = "", quote="")

test3 <- read.table("/Users/brian/test/5104995_23153_0_0.final.tsv",
                    header = T, sep = "\t", comment.char = "", quote="")





library(chron)
library(tidyverse)
germline <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/brian.germline.txt",
                       header = F, sep = "\t", comment.char = "", quote="") 
germline <- read.table("/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/kelly.germline.txt",
                       header = F, sep = "\t", comment.char = "", quote="") 
germline <- germline %>% separate(V1, c("status", "nothing", "date", "time", "size", "size_unit", "blank", "id"), 
                       sep=" ",extra = "merge", fill = "right") %>%
 select(-c("nothing", "blank"))


germline$id <- trimws(germline$id , which = "left")
germline <- germline %>% separate(id, c("file_name", "id"), 
                                  sep=" ", extra = "merge", fill = "right")
germline$id <- str_replace(germline$id, "\\(", "")
germline$id <- str_replace(germline$id, "\\)", "")
germline <- germline[order(germline$file_name), ]
germline$date_time <- chron(dates=germline$date, times=germline$time, format=c('y-m-d','h:m:s'))
germline.t <- germline %>%
  group_by(file_name) %>%
  filter(date_time == max(date_time))

write.table(germline.t$file_name, "/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/brian.eids.txt",
            col.names = F, row.names = F, quote = F)

write.table(germline.t$file_name, "/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/kelly.eids.txt",
            col.names = F, row.names = F, quote = F)
  