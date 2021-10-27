library(dplyr)
set.seed(seed = 14412)

cancer <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbCancerRegistry_09232021.rds")
lung <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbLungCancerCaCo_10132021.rds")
ukbb <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbbExome_09232021.rds")


transplant <- read.csv("/Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv")
sort(table(transplant$eid), decreasing = T)[sort(table(transplant$eid), decreasing = T) >1]
# transplant <- transplant %>% dplyr::select(eid, ageBaseline)
## DNAnexus only uses floor age
#transplant$ageBaseline <- round(transplant$ageBaseline) 
transplant$ageBaseline <- floor(transplant$ageBaseline)
transplant <- transplant[order(transplant$ageBaseline),]
transplant <- distinct(transplant, eid, .keep_all = TRUE)
transplant$match_id <- 1:332
# tt <- read.table("/Users/brian/Bolton/UKBB/tranplant.test.txt")
transplant.match <- read.table("~/Bolton/UKBB/results/12/update/transplant_match/match.txt", sep = "\t")
transplant.case.control.eids <- c(unique(transplant$eid),unique(transplant.match$V1))
sum(transplant.case.control.eids %in% lung$eid)
1239160 %in% transplant.case.control.eids

all.samples <- read.csv("~/Bolton/UKBB/data.csv")
colnames(all.samples)[2] <- "ageBaseline"
all.samples <- all.samples %>% dplyr::select(-p23153)


samples.ran.including.transplant <- read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt")
## also ran 479 samples from batch folder 13
batch.13 <- read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/13/samples.completed.txt")
colnames(samples.ran.including.transplant) <- "eid"
colnames(batch.13) <- "eid"
all.completed <- c(samples.ran.including.transplant$eid, batch.13$eid)

samples.ran.including.transplant <- samples.ran.including.transplant %>% left_join(all.samples, by = "eid") 
length(setdiff(samples.ran.including.transplant$eid, transplant.case.control.eids))+length(transplant.case.control.eids) # 2153
length(unique(samples.ran.including.transplant$eid)) # 2153
sum(unique(all.completed) %in% lung$eid) # 65

sum(all.completed %in% lung$eid)
## difference is [1] 1228453 1234566 1236681
# ls /Volumes/bolton/Active/projects/mocha/UKBB/exome/results/12/update/transplant_match/1228453*.tsv
# setdiff(unique(samples.ran.including.transplant$eid)[unique(samples.ran.including.transplant$eid) %in% lung$eid],
#         unique(transplant.case.control.eids)[transplant.case.control.eids %in% lung$eid])
lung <- lung %>% left_join(samples.ran.including.transplant %>% dplyr::select(-ageBaseline))
lung$ran <- ifelse(lung$eid %in% all.completed, 1, 0)
lung$ran[is.na(lung$ran)] <- 0
write.table(lung, "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/lung/all_lung.txt", sep = '\t', row.names = F, quote = F)
write.table(lung[!lung$ran,"eid"], "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/FINAL_COMBINE_JOBS_LUNG/lung_eids_remaining.txt", sep = '\t', 
            row.names = F, quote = F, col.names = F)
length(unique(stringr::str_extract(lung$eid[!lung$ran],"[0-9]{2}")))

# total.12 <- rbind(read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/12/brian_project.no_uniq"),
#                   read.table("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/12/kelly_project.no_uniq"))
# intersect(total.12$V1, transplant$eid)

samples.ran.excluding.transplant <- samples.ran.including.transplant[!(samples.ran.including.transplant$eid %in% transplant$eid),]
samples.ran.excluding.transplant <- samples.ran.excluding.transplant[order(samples.ran.excluding.transplant$ageBaseline),]
## excludes transplant
samples.ran.excluding.transplant <- samples.ran.excluding.transplant %>% dplyr::select(eid, ageBaseline)


buckets.by.age <- table(transplant$ageBaseline) * 1
final <- data.frame(buckets.by.age)
final$total.match.needed <- final$Freq * 5
final$total.match.have <- 0
final$total.match.have.eid.list <- 0


for (age in unique(transplant$ageBaseline)) {
  age.samples.completed <- samples.ran.excluding.transplant %>% 
    filter(ageBaseline==age)
  final[final$Var1==age,]$total.match.have <- nrow(age.samples.completed)
  age.samples.completed <- samples.ran.excluding.transplant %>% 
    filter(ageBaseline==age) %>%
    top_n(-final[final$Var1==age,]$total.match.needed, eid)
  final[final$Var1==age,]$total.match.have.eid.list <- list(age.samples.completed$eid)
}
sum(apply(final,1,function(x) {
  print(length(x[["total.match.have.eid.list"]]))
}))
matched_eids <- unlist(apply(final,1,function(x) {
  x[["total.match.have.eid.list"]]
}))
# sum(rep(final$Var1[21], final$total.match.needed)==60)
matched_eids <- data.frame(eid=matched_eids)
matched_eids <- matched_eids %>% dplyr::left_join(all.samples, by = c("eid"="eid")) 
matched_eids$match_id <- rep(1:332,each=5)
write.table(matched_eids, "~/Bolton/UKBB/results/12/update/transplant_match/match.txt", sep = "\t", col.names = F, row.names = F)

sum(transplant$eid %in% cancer$eid)/length(transplant$eid)
sum(matched_eids$eid %in% cancer$eid)/length(matched_eids$eid)


eids <- sort(unique(c(matched_eids$eid, transplant$eid)))
eids2 <- sort(unique(gsub("_.*","",total.passed.pon.nsamp$eid)))
all(eids == eids2)
transplant$transplant <- 1
matched_eids$transplant <- 0
write.table(rbind(transplant, matched_eids), "~/Bolton/UKBB/transplant_match.txt", sep = "\t")
# P <- read.table("~/Bolton/UKBB/results/12/update/transplant_match/match.txt")
# 
# final$have.enough <- final$total.match.have >= final$total.match.needed
# final$number.of.age.needed <- ifelse(final$total.match.have >= final$total.match.needed, 0, final$total.match.needed-final$total.match.have)
# colnames(final)[1] <- "transplant_age_bin"
# 
# sum(all.samples$eid %in% samples.ran.including.transplant$eid)
# all.samples.leftToRun <- all.samples %>% filter(!(eid %in% samples.ran.including.transplant$eid))
# 
# samples.remain.12 <- all.samples.leftToRun %>% filter(grepl("^12",eid)) %>% select(-p23153)
# final.needed <- final %>% filter(!have.enough)
# 
# run.eids <- data.frame()
# for (i in 1:nrow(final.needed)) {
#   #age.filter <- samples.remain.12 %>% filter(ageBaseline==age)
#   fn <- final.needed[i,]
#   age.filter.remain <- samples.remain.12 %>% filter(ageBaseline==fn$transplant_age_bin)
#   ## default replace = FALSE
#   ## already scrambled identifiers so no need for random sample
#   # run.eids <- rbind(run.eids, dplyr::sample_n(age.filter.remain, fn$number.of.age.needed))
#   # top_n(-2) # lowest values would be increasing so this is the lowest batchest of group 12
#   run.eids <- rbind(run.eids, dplyr::top_n(age.filter.remain, -fn$number.of.age.needed, eid))
# }
# 
# final.needed2 <- final.needed %>% select(transplant_age_bin, number.of.age.needed)
# r2 <- run.eids %>% group_by(ageBaseline) %>% summarise(n=n())
# r2==final.needed2


##############################################################################
#############  TSV Outputs
##############################################################################
library(tibble)

missing.cols <- c("Vardict_calpos",
                  "Vardict_calref",
                  "Vardict_calalt",
                  "Vardict_calpos_best",
                  "Vardict_calref_best",
                  "Vardict_calalt_best")
transplant.257 <- read.table("~/Bolton/UKBB/results/12/update/transplant/257.txt")
columns.263 <- colnames(read.table("~/Bolton/UKBB/results/12/update/transplant/1006140_23153_0_0.final.tsv",
                                   header = T, sep = "\t", comment.char = "", quote=""))
columns.2632 <- colnames(read.table("~/Bolton/UKBB/results/12/update/transplant_match/match/1200019_23153_0_0.final.tsv",
                                   header = T, sep = "\t", comment.char = "", quote=""))

# missing Alex's columns
# for (file in transplant.257$V1) {
#   sample.257 <-  read.table(paste0("~/Bolton/UKBB/results/12/update/transplant/", file), 
#                             header = T, sep = "\t", comment.char = "", quote="")
#   for (column in missing.cols) {
#     sample.257[[column]] <- NA
#   }
#   print(all(colnames(sample.257[,columns.263]) == columns.263))
#   write.table(sample.257[,columns.263], paste0("~/Bolton/UKBB/results/12/update/transplant/", file), 
#               col.names = T, quote = F, row.names = F, sep = "\t")
# }
transplant.rearrange <- read.table("~/Bolton/UKBB/results/12/update/transplant_match/match/rearrange.txt")
i=0
for (file in transplant.rearrange$V1) {
  sample.rearrange <-  read.table(paste0("~/Bolton/UKBB/results/12/update/transplant_match/match/", file),
                            header = T, sep = "\t", comment.char = "", quote="")
  print(all(colnames(sample.rearrange[,columns.263]) == columns.263))
  write.table(sample.rearrange[,columns.263], paste0("~/Bolton/UKBB/results/12/update/transplant_match/match/", file),
              col.names = T, quote = F, row.names = F, sep = "\t")
  i <- i+1
  print(i)
}


transplant.output.passed.pon <- read.table("~/Bolton/UKBB/results/12/update/transplant/combined.UKBB.12.pon.tsv", 
                                      header = T, sep = "\t", comment.char = "", quote="")
transplant.output <- read.table("~/Bolton/UKBB/results/12/update/transplant/combined.UKBB.12.tsv", 
                                           header = T, sep = "\t", comment.char = "", quote="")
transplant.output.passed.pon <- transplant.output %>% filter(PON_FISHER < 0.05/38793002)
transplant.output$received.transplant <- 1
bf.correction <- 0.05/38793002
transplant.output$eid2
left_join(all.samples, by = "eid") 
transplant.output %>% left_join(all.samples, by = c("eid2"=="eid"))


# for file in $(ls *23153_0_0.final.tsv); do NF=$(head -n1 $file | awk -F'\t' '{print NF}'); printf "$file\t"; echo $NF; done > colnums.txt
# awk '{print $2}' colnums.txt | sort | uniq
# ggrep -P "\t257" colnums.txt > 257.txt
# for file in $(ls ~/Bolton/UKBB/results/12/update/transplant_match/); do NF=$(head -n1 $file | awk -F'\t' '{print NF}'); printf "$file\t"; echo $NF; done > colnums.txt, file
transplant_match.257 <- read.table("~/Bolton/UKBB/results/12/update/transplant_match/257.txt")
for (file in transplant_match.257$V1) {
  print(paste0("~/Bolton/UKBB/results/12/update/transplant_match/", file))
  sample.257 <-  read.table(paste0("~/Bolton/UKBB/results/12/update/transplant_match/", file),
                            header = T, sep = "\t", comment.char = "", quote="")
  for (column in missing.cols) {
    sample.257[[column]] <- NA
  }
  print(all(colnames(sample.257[,columns.263]) == columns.263))
  write.table(sample.257[,columns.263], paste0("~/Bolton/UKBB/results/12/update/transplant_match/", file),
              col.names = T, quote = F, row.names = F, sep = "\t")
}
transplant.match.output.passed.pon <- read.table("~/Bolton/UKBB/results/12/update/transplant_match/match/combined.UKBB.12.pon.tsv", 
                                      header = T, sep = "\t", comment.char = "", quote="")


transplant.failed.pon <- transplant.output %>% filter(PON_FISHER >= bf.correction)
transplant.output.passed.pon <- transplant.output %>% filter(PON_FISHER < bf.correction)
sum(is.na(transplant.output$PON_FISHER))
fisher.na <- transplant.output[is.na(transplant.output$PON_FISHER) & (transplant.output$Mutect2_CALLER),]

#!/usr/local/bin/Rscript

library(argparser)
library(dplyr)
library(data.table)

parser2 <- arg_parser("Filter NSamples")
parser2 <- add_argument(parser2, "--number", type="integer", help="filter out > number")
parser2 <- add_argument(parser2, "--column", type="character", help="which column to group by, usually column with <chrom:pos:ref:alt>")
parser2 <- add_argument(parser2, "--input", type="character", help="in")
parser2 <- add_argument(parser2, "--output", type="character", help="out")
args2 <- parse_args(parser2)
args2 <- parse_args(parser2, list(c("-n","27"),
                                c("-c","key")))


# args2$input <- "/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.tsv"
args2$input <- "/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.tsv"
# args2$input <- "~/Bolton/HIV/final_all_rows_all_filters.tsv"
args2$output <- "/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.nsamp.tsv"

transplant.match.output.passed.pon <- fread(file=args2$input, header = T, sep = "\t", quote="")
transplant.match.output.passed.pon$key <- with(transplant.match.output.passed.pon, paste(CHROM,POS,REF,ALT,sep=":"))

total.passed.pon.nsamp <- transplant.match.output.passed.pon %>% 
  dplyr::group_by(key) %>% 
  dplyr::mutate(num=dplyr::n()) %>%
  ungroup() %>%
  filter(num<=args2$number)
  

write.table(total.passed.pon.nsamp, args2$output, sep = "\t", row.names = F, quote = F)

total.passed.pon.nsamp1 <- fread(file="/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.nsamp.tsv",
                                header = T, sep = "\t", quote="")
total.passed.pon.nsamp2 <- fread(file="/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.nsamp.updated.sorted.tsv",
                                header = T, sep = "\t", quote="")
total.passed.pon.nsamp2$pan_my <- apply(total.passed.pon.nsamp2[,c("MDS","AML","MPN")], 1, sum)
hotspots <- total.passed.pon.nsamp2[total.passed.pon.nsamp2$called & # called by both callers
                                      !total.passed.pon.nsamp2$Mutect2_PASS & # didn't pass mutect
                                      (total.passed.pon.nsamp2$heme_cosmic_count >= 10 | # either heme count >= 10
                                         total.passed.pon.nsamp2$myeloid_cosmic_count >= 5 | # or myeloid count >=5
                                         total.passed.pon.nsamp2$n.loci.vep >= 5 | # or bick/bolton loci >= 5
                                         total.passed.pon.nsamp2$pan_my >= 5),] # or pan myeloid >=5

transplant.ch_pd2 <- total.passed.pon.nsamp2[(total.passed.pon.nsamp2$ch_pd2 & # annotated ch_pd2 from the annotate_PD script
                           total.passed.pon.nsamp2$Mutect2_PASS & # Passed by Mutect
                           total.passed.pon.nsamp2$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
                           total.passed.pon.nsamp2$max.over.0.02 & # At least 2% VAF
                           (total.passed.pon.nsamp2$max.under.0.35 | # Under 35% VAF 
                              (!total.passed.pon.nsamp2$max.under.0.35 & (!is.na(total.passed.pon.nsamp2$n.loci.vep)))) & # or if over 35% is hotspot
                           !total.passed.pon.nsamp2$Vardict_PON_2AT2_percent & # Passed PON2
                           !total.passed.pon.nsamp2$Mutect2_PON_2AT2_percent),] # Passed PON2
transplant.ch_pd2 <- transplant.ch_pd2[, c(1:252,272,270,271,257,256,259:269,253:255,258)]
hotspots <- hotspots[, c(1:252,272,270,271,257,256,259:269,253:255,258)]

# for (i in 1:length(colnames(transplant.ch_pd2))) { print(paste0(i,": ",colnames(transplant.ch_pd2)[i]))}

transplant.ch_pd2 <- sqldf("SELECT l.*, r.`source.totals`
            FROM `transplant.ch_pd2` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
transplant.ch_pd2 <- transplant.ch_pd2[!duplicated(transplant.ch_pd2),]
transplant.ch_pd2$n.loci.vep<- fillna(transplant.ch_pd2$n.loci.vep, 0)
transplant.ch_pd2 <- sqldf("SELECT l.*, r.`n.HGVSp`
            FROM `transplant.ch_pd2` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSp_VEP = r.HGVSp_VEP)")
transplant.ch_pd2 <- transplant.ch_pd2[!duplicated(transplant.ch_pd2),]
transplant.ch_pd2 <- sqldf("SELECT l.*, r.`n.HGVSc`
            FROM `transplant.ch_pd2` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSc_VEP = r.HGVSc_VEP)")
transplant.ch_pd2 <- transplant.ch_pd2[!duplicated(transplant.ch_pd2),]

hotspots <- sqldf("SELECT l.*, r.`source.totals`
            FROM `hotspots` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
hotspots <- hotspots[!duplicated(hotspots),]
hotspots$n.loci.vep<- fillna(hotspots$n.loci.vep, 0)
hotspots <- sqldf("SELECT l.*, r.`n.HGVSp`
            FROM `hotspots` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSp_VEP = r.HGVSp_VEP)")
hotspots <- hotspots[!duplicated(hotspots),]
hotspots <- sqldf("SELECT l.*, r.`n.HGVSc`
            FROM `hotspots` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSc_VEP = r.HGVSc_VEP)")
hotspots <- hotspots[!duplicated(hotspots),]


transplant.ch_pd2$gnomAD_MAX.Stringent.005 <- (transplant.ch_pd2$MAX_gnomAD_AF_VEP < 0.005 &
                                                 transplant.ch_pd2$MAX_gnomADe_AF_VEP < 0.005 &
                                                 transplant.ch_pd2$MAX_gnomADg_AF_VEP < 0.005)
hotspots$gnomAD_MAX.Stringent.005 <- (hotspots$MAX_gnomAD_AF_VEP < 0.005 &
                                        hotspots$MAX_gnomADe_AF_VEP < 0.005 &
                                        hotspots$MAX_gnomADg_AF_VEP < 0.005)

## if False means all 3 groups has a Max VAF > .0007 so this is less stringent for filtering out because
## all 3 needs to have a Max above .005
transplant.ch_pd2$gnomAD_MAX.lessStringent.005 <- (transplant.ch_pd2$MAX_gnomAD_AF_VEP < 0.005 |
                                                     transplant.ch_pd2$MAX_gnomADe_AF_VEP < 0.005 |
                                                     transplant.ch_pd2$MAX_gnomADg_AF_VEP < 0.005)
hotspots$gnomAD_MAX.lessStringent.005 <- (hotspots$MAX_gnomAD_AF_VEP < 0.005 |
                                            hotspots$MAX_gnomADe_AF_VEP < 0.005 |
                                            hotspots$MAX_gnomADg_AF_VEP < 0.005)

# test.before <- fread(file="/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB_BW_FINAL2.txt", 
#                      header = T, sep = "\t", quote="")
# setdiff(test.before$CHROM.POS.REF.ALT, transplant.ch_pd2$key)
# total.passed.pon.nsamp2[1453254,]
# total.passed.pon.nsamp2$ch_pd2 & # annotated ch_pd2 from the annotate_PD script
#  total.passed.pon.nsamp2$Mutect2_PASS & # Passed by Mutect
#  total.passed.pon.nsamp2$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#  total.passed.pon.nsamp2$max.over.0.02 & # At least 2% VAF
#  (total.passed.pon.nsamp2$max.under.0.35 | # Under 35% VAF 
#     (!total.passed.pon.nsamp2$max.under.0.35 & (!is.na(total.passed.pon.nsamp2$n.loci.vep)))) & # or if over 35% is hotspot
#  !total.passed.pon.nsamp2$Vardict_PON_2AT2_percent & # Passed PON2
#  !total.passed.pon.nsamp2$Mutect2_PON_2AT2_percent # Passed PON2



# transplant.ch_pd2 <- total.passed.pon.nsamp[(total.passed.pon.nsamp$ch_pd2 & # annotated ch_pd2 from the annotate_PD script
#                                                total.passed.pon.nsamp$Mutect2_PASS & # Passed by Mutect
#                                                total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#                                                total.passed.pon.nsamp$max.over.0.02 & # At least 2% VAF
#                                                (total.passed.pon.nsamp$max.under.0.35 | # Under 35% VAF 
#                                                   (!total.passed.pon.nsamp$max.under.0.35 & (!is.na(total.passed.pon.nsamp$n.loci.vep))))  # or if over 35% is hotspot
#                                                ),]
# 
# transplant.output.passed.pon.keys <- with(transplant.output.passed.pon,
#                                           paste(CHROM, POS, REF, ALT, sep = ":"))
# transplant.output.passed.pon.CPRA.keys <- with(transplant.output.passed.pon.CPRA,
#                                                paste(CHROM, POS, REF, ALT, sep = ":"))
# 
# transplant.output.passed.pon.nsamples <- transplant.output.passed.pon[transplant.output.passed.pon.keys %in% transplant.output.passed.pon.CPRA.keys,]
# # transplant.output.passed.pon.nsamples[(!transplant.output.passed.pon.nsamples$ch_pd2.y &
# #                                          (transplant.output.passed.pon.nsamples$called & transplant.output.passed.pon.nsamples$alt_strand_counts_min_2_callers)
# #                                        ),]
# transplant.output.passed.pon.nsamples$Vardict_CALLER <- as.numeric(transplant.output.passed.pon.nsamples$Vardict_CALLER)
# transplant.output.passed.pon.nsamples$Mutect2_PON_2AT2_percent <- as.numeric(transplant.output.passed.pon.nsamples$Mutect2_PON_2AT2_percent)
# test <- transplant.output.passed.pon.nsamples[(!transplant.output.passed.pon.nsamples$ch_pd2.y &
#                                          ((transplant.output.passed.pon.nsamples$passed & transplant.output.passed.pon.nsamples$alt_strand_counts_min_2_callers) |
#                                          (transplant.output.passed.pon.nsamples$Mutect2_PASS & !transplant.output.passed.pon.nsamples$Vardict_CALLER & transplant.output.passed.pon.nsamples$Mutect2_SB)) &
#                                          (transplant.output.passed.pon.nsamples$complexity_filters & transplant.output.passed.pon.nsamples$min.under.0.35 & transplant.output.passed.pon.nsamples$max.over.0.02 &
#                                             !transplant.output.passed.pon.nsamples$Vardict_PON_2AT2_percent & !transplant.output.passed.pon.nsamples$Mutect2_PON_2AT2_percent)
#                                           
#                                            
# ),]
# 
# transplant.no_ch_pd2 <- transplant.output.passed.pon.nsamples[(!transplant.output.passed.pon.nsamples$ch_pd2.y & 
#                                                  transplant.output.passed.pon.nsamples$passed_everything_mutect & 
#                                                  !(transplant.output.passed.pon.nsamples$blacklist) &
#                                                  !(transplant.output.passed.pon.nsamples$dups) &
#                                                  !(transplant.output.passed.pon.nsamples$simpleTandemRepeat) &
#                                                  !(transplant.output.passed.pon.nsamples$gnomAD_MAX.Stringent.0007)),]
# 
# transplant.ch_pd2 <- transplant.output.passed.pon.nsamples[(transplant.output.passed.pon.nsamples$ch_pd2.y &
#                                                  ((transplant.output.passed.pon.nsamples & transplant.output.passed.pon.nsamples$alt_strand_counts_min_2_callers) |
#                                                     (transplant.output.passed.pon.nsamples$Mutect2_CALLER & !transplant.output.passed.pon.nsamples$Vardict_CALLER & transplant.output.passed.pon.nsamples$Mutect2_SB)) &
#                                                  transplant.output.passed.pon.nsamples$max.over.0.02 &
#                                                     (!transplant.output.passed.pon.nsamples$Vardict_PON_2AT2_percent & !transplant.output.passed.pon.nsamples$Mutect2_PON_2AT2_percent)
# ),]

## NSAMPLES
# total.passed.pon <- read.table("/Users/brian/Bolton/UKBB/results/12/update/combined/combined.UKBB.12.tsv.key", 
#                     header = T, sep = "\t", comment.char = "", quote="")
# test3 <- data.frame(total.passed.pon$CHROM.POS.REF.ALT)
# test4 <- test3 %>% group_by(total.passed.pon.CHROM.POS.REF.ALT) %>% 
#   summarise(n=n()) %>%
#   filter(n<=27)
# total.passed.pon.nsamp <- total.passed.pon[test3$total.passed.pon.CHROM.POS.REF.ALT %in% test4$total.passed.pon.CHROM.POS.REF.ALT,]



















##########################################
# total.passed.pon.nsamp2 <- read.table("/Users/brian/Bolton/UKBB/results/12/update/combined/combined.UKBB.pon.nsamp.tsv",
#                                                         header = T, sep = "\t", comment.char = "", quote="")
# total.passed.pon.nsamp2$VariantClass <- ifelse(total.passed.pon.nsamp2$VariantClass=="Missense","Missense_Mutation",total.passed.pon.nsamp2$VariantClass)
# total.passed.pon.nsamp2 <- total.passed.pon.nsamp2 %>%
#   dplyr::mutate(loci = str_extract(AAchange.y, "[A-Z]\\d+"),
#                 AAchange = str_remove(AAchange.y, "p."),
#                 gene_aachange = paste(Gene.y,AAchange,sep = "_"),
#                 gene_loci = paste(Gene.y,loci,sep = "_"),
#                 ch_pd2 = case_when(ch_pd==1 ~ 1,
#                                    (gene_loci %in% topmed.loci.n$gene_loci | gene_aachange %in% topmed.mutation.1$gene_aachange) & VariantClass=="Missense_Mutation" ~ 1,
#                                    gene_loci %in% topmed.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
#                                    (gene_loci %in% kelly.loci.n$gene_loci | gene_aachange %in% kelly.mutation.1$gene_aachange) & VariantClass=="Missense_Mutation" ~ 1,
#                                    gene_loci %in% kelly.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
#                                    TRUE ~ 0))
# total.passed.pon.nsamp2$Vardict_CALLER <- as.numeric(total.passed.pon.nsamp2$Vardict_CALLER)
# total.passed.pon.nsamp2$Mutect2_PON_2AT2_percent <- as.numeric(total.passed.pon.nsamp2$Mutect2_PON_2AT2_percent)
# total.passed.pon.nsamp2$alt_strand_counts_min_1_caller_only <- as.logical(total.passed.pon.nsamp2$alt_strand_counts_min_1_caller_only)
# total.passed.pon.nsamp2$max.over.0.02 <- as.logical(total.passed.pon.nsamp2$max.over.0.02)
# total.passed.pon.nsamp2$min.under.0.35 <- as.logical(total.passed.pon.nsamp2$min.under.0.35)
# total.passed.pon.nsamp2$max.under.0.35 <- apply(total.passed.pon.nsamp2[,getVAFs(total.passed.pon.nsamp2)],1,max) < 0.35
# total.passed.pon.nsamp2$complexity_filters <- as.logical(total.passed.pon.nsamp2$complexity_filters)
# 
# transplant.ch_pd2 <- total.passed.pon.nsamp2[(total.passed.pon.nsamp2$ch_pd2.y & # annotated ch_pd2 from the annotate_PD script
#                                                 total.passed.pon.nsamp2$Mutect2_PASS & # Passed by Mutect
#                                                 total.passed.pon.nsamp2$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#                                                 total.passed.pon.nsamp2$max.over.0.02 & # At least 2% VAF
#                                                (total.passed.pon.nsamp2$min.under.0.35 | # Under 35% VAF 
#                                                   (!total.passed.pon.nsamp2$min.under.0.35 & total.passed.pon.nsamp2$CHROM.POS.REF.ALT %in% vars$key)) & # or if over 35% is hotspot
#                                                !total.passed.pon.nsamp2$Vardict_PON_2AT2_percent & # Passed PON2
#                                                !total.passed.pon.nsamp2$Mutect2_PON_2AT2_percent),] # Passed PON2
# test <- read.table("/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB3.tsv",
#                    header = T, sep = "\t", comment.char = "", quote="")
# paste(transplant.ch_pd2$CHROM.POS.REF.ALT, transplant.ch_pd2$eid,sep = ":")
# paste(test$CHROM.POS.REF.ALT, test$eid,sep = ":")[!paste(test$CHROM.POS.REF.ALT, test$eid,sep = ":") %in% paste(transplant.ch_pd2$CHROM.POS.REF.ALT, transplant.ch_pd2$eid,sep = ":")]
##########################################

# total.passed.pon.nsamp2 <- read.table("/Users/brian/Bolton/UKBB/results/12/update/combined/combined.UKBB.pon.nsamp.tsv",
#                                       header = T, sep = "\t", comment.char = "", quote="")
# total.passed.pon.nsamp <- total.passed.pon.nsamp2[!(total.passed.pon.nsamp2$PON_AltDepth/(total.passed.pon.nsamp2$PON_AltDepth+total.passed.pon.nsamp2$PON_RefDepth) >= 
#                                                       total.passed.pon.nsamp2$Mutect2_gt_AF),]
# total.passed.pon.nsamp$VariantClass <- ifelse(total.passed.pon.nsamp$VariantClass=="Missense","Missense_Mutation",total.passed.pon.nsamp$VariantClass)
# source("/Users/brian/Bolton/R_stuff/hotspots.R")
# 
# total.passed.pon.nsamp <- total.passed.pon.nsamp %>%
#   dplyr::mutate(loci = str_extract(AAchange.y, "[A-Z]\\d+"),
#                 AAchange = str_remove(AAchange.y, "p."),
#                 gene_aachange = paste(Gene.y,AAchange,sep = "_"),
#                 gene_loci = paste(Gene.y,loci,sep = "_"),
#                 ch_pd2 = case_when(ch_pd==1 ~ 1,
#                                    (gene_loci %in% topmed.loci.n$gene_loci | gene_aachange %in% topmed.mutation.1$gene_aachange) & VariantClass=="Missense_Mutation" ~ 1,
#                                    gene_loci %in% topmed.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
#                                    (gene_loci %in% kelly.loci.n$gene_loci | gene_aachange %in% kelly.mutation.1$gene_aachange) & VariantClass=="Missense_Mutation" ~ 1,
#                                    gene_loci %in% kelly.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
#                                    TRUE ~ 0))
# 
# total.passed.pon.nsamp$Vardict_CALLER <- as.numeric(total.passed.pon.nsamp$Vardict_CALLER)
# total.passed.pon.nsamp$Mutect2_PON_2AT2_percent <- as.numeric(total.passed.pon.nsamp$Mutect2_PON_2AT2_percent)
# total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only <- as.logical(total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only)
# total.passed.pon.nsamp$max.over.0.02 <- as.logical(total.passed.pon.nsamp$max.over.0.02)
# total.passed.pon.nsamp$min.under.0.35 <- as.logical(total.passed.pon.nsamp$min.under.0.35)
# total.passed.pon.nsamp$max.under.0.35 <- apply(total.passed.pon.nsamp[,getVAFs(total.passed.pon.nsamp)],1,max) < 0.35
# total.passed.pon.nsamp$complexity_filters <- as.logical(total.passed.pon.nsamp$complexity_filters)
# total.passed.pon.nsamp <- total.passed.pon.nsamp %>% dplyr::left_join(vars %>% dplyr::select(n,gene_loci, source),
#                                                                                                     by = "gene_loci")
## 1 passes SB and 0 fails SB
# Notes:
#MUTS$ch_pd = ifelse(MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic", 1,0)
#MUTS$ch_pd = ifelse(MUTS$Gene %in% TSG & MUTS$VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"), 1, MUTS$ch_pd)
# transplant.ch_pd2 <- total.passed.pon.nsamp[(total.passed.pon.nsamp$Mutect2_PASS & # Passed by Mutect
#                                                total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#                                                total.passed.pon.nsamp$max.over.0.02 & # At least 2% VAF
#                                                (total.passed.pon.nsamp$max.under.0.35 | # Under 35% VAF 
#                                                   (!total.passed.pon.nsamp$max.under.0.35 & total.passed.pon.nsamp$CHROM.POS.REF.ALT %in% vars$key)) & # or if over 35% is hotspot
#                                                !total.passed.pon.nsamp$Vardict_PON_2AT2_percent & # Passed PON2
#                                                !total.passed.pon.nsamp$Mutect2_PON_2AT2_percent & # Passed PON2
#                                                total.passed.pon.nsamp$complexity_filters),] 
# 
# 
# transplant.ch_pd2.called <- total.passed.pon.nsamp[(total.passed.pon.nsamp$Mutect2_CALLER & # called by Mutect
#                                                       !total.passed.pon.nsamp$Mutect2_PASS & # NOT passed by Mutect
#                                                       total.passed.pon.nsamp$Vardict_CALLER & # called by Vardict
#                                                       total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#                                                       total.passed.pon.nsamp$max.over.0.02 & # At least 2% VAF
#                                                       total.passed.pon.nsamp$gene_loci %in% vars$gene_loci & # or if over 35% is hotspot
#                                                       !total.passed.pon.nsamp$Vardict_PON_2AT2_percent & # Passed PON2
#                                                       !total.passed.pon.nsamp$Mutect2_PON_2AT2_percent),] # Passed PON2
# transplant.ch_pd2.called <- transplant.ch_pd2.called %>% dplyr::left_join(vars %>% dplyr::select(n,gene_loci, source),
#                                          by = "gene_loci")
# transplant.ch_pd2.called$key2 <- with(transplant.ch_pd2.called,paste(CHROM.POS.REF.ALT,eid,sep = ":"))
# transplant.ch_pd2.called <- transplant.ch_pd2.called[!duplicated(transplant.ch_pd2.called$key2),]
# 
# transplant.ch_pd2.for.cosmic <- total.passed.pon.nsamp[(total.passed.pon.nsamp$Mutect2_CALLER & # Passed by Mutect
#                                                           !total.passed.pon.nsamp$Mutect2_PASS &
#                                                         total.passed.pon.nsamp$Vardict_CALLER & # Passed by Mutect
#                                                           total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#                                                           total.passed.pon.nsamp$max.over.0.02 & # At least 2% VAF
#                                                           (total.passed.pon.nsamp$max.under.0.35 | # Under 35% VAF 
#                                                              (!total.passed.pon.nsamp$max.under.0.35 & total.passed.pon.nsamp$CHROM.POS.REF.ALT %in% vars$key)) & # or if over 35% is hotspot
#                                                           !total.passed.pon.nsamp$Vardict_PON_2AT2_percent & # Passed PON2
#                                                           !total.passed.pon.nsamp$Mutect2_PON_2AT2_percent),] # Passed PON2
# sum(transplant.ch_pd2.for.cosmic$PON_AltDepth/(transplant.ch_pd2.for.cosmic$PON_AltDepth+transplant.ch_pd2.for.cosmic$PON_RefDepth) >= 
#       transplant.ch_pd2.for.cosmic$Mutect2_gt_AF)
# transplant.ch_pd2.for.cosmic <- transplant.ch_pd2.for.cosmic[!((grepl("synonymous_variant", transplant.ch_pd2.for.cosmic$Consequence_VEP) & 
#                                                         !(grepl("splice", transplant.ch_pd2.for.cosmic$Consequence_VEP)) & 
#                                                         !(grepl("NMD_transcript", transplant.ch_pd2.for.cosmic$Consequence_VEP)))),]
# transplant.ch_pd2.for.cosmic$key2 <- with(transplant.ch_pd2.for.cosmic,paste(CHROM.POS.REF.ALT,eid,sep = ":"))
# transplant.ch_pd2.called <- transplant.ch_pd2.called %>% dplyr::left_join(vars %>% dplyr::select(n,gene_loci, source),
#                                                                           by = "gene_loci")
# write.table(transplant.ch_pd2.for.cosmic,"/Volumes/bolton/Active/projects/annotation_files/cosmic/transplant.ch_pd2.for.cosmic.big.tsv",
#             col.names = T, row.names = F, sep = "\t", quote = F)
# 
# 
# cosmic2.annotations <- read.table("/Users/brian/Bolton/UKBB/cosmic/cosmic_count.called.total.tsv",
#                                   header=T,sep = "\t", quote = "")   
# cosmic2.annotations$key2 <-  with(cosmic2.annotations,paste(var_key,eid,sep = ":"))
# row.names(cosmic2.annotations) <- cosmic2.annotations$key2
# transplant.ch_pd2.for.cosmic.final <- cbind(transplant.ch_pd2.for.cosmic,cosmic2.annotations[transplant.ch_pd2.for.cosmic$key2,c(4,5,6,7,8,9)])
# test2 <- transplant.ch_pd2.for.cosmic.final[transplant.ch_pd2.for.cosmic.final$heme_cosmic_count2>=10,] %>%
#   dplyr::select(CHROM,POS,REF,ALT,eid,HGVSp_VEP,Existing_variation_VEP,CHROM.POS.REF.ALT,ID,var_key_cosmic,HGVSP,CosmicCount2,heme_cosmic_count2,myeloid_cosmic_count)
# 
# test <- transplant.ch_pd2.for.cosmic.final %>%
#   filter(grepl("COSV", Existing_variation_VEP),!grepl("COSV", ID)) %>%
#   dplyr::select(CHROM,POS,REF,ALT,eid,HGVSp_VEP,Existing_variation_VEP,CHROM.POS.REF.ALT,ID,var_key_cosmic,HGVSP)





# sum(transplant.ch_pd2.for.cosmic$CHROM.POS.REF.ALT %in% transplant.ch_pd2.called$CHROM.POS.REF.ALT)
# test <- total.passed.pon.nsamp[(!(grepl("synonymous_variant", total.passed.pon.nsamp$Consequence_VEP) & 
#                                     !(grepl("splice", total.passed.pon.nsamp$Consequence_VEP)))),]

# transplant.ch_pd2.for.cosmic.for.join <- transplant.ch_pd2.for.cosmic %>%
#   dplyr::select(CHROM,eid,CHROM.POS.REF.ALT,HGVSp_VEP) %>%
#   rename_at(vars(CHROM.POS.REF.ALT), function(x) x="var_key")
# write.table(transplant.ch_pd2.for.cosmic.for.join, "/Volumes/bolton/Active/projects/annotation_files/cosmic/transplant.ch_pd2.for.cosmic.for.join.tsv",
#             col.names = T, row.names = F, sep = "\t", quote = F)
# transplant.ch_pd2.for.cosmic.for.join <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/transplant.ch_pd2.for.cosmic.for.join.tsv",
#                                                     header=T,sep = "\t", quote = "")   
# 
# 
# transplant.ch_pd2.2 <- total.passed.pon.nsamp[(total.passed.pon.nsamp$ch_pd2.y & # annotated ch_pd2 from the annotate_PD script
#                                              total.passed.pon.nsamp$Mutect2_PASS & # Passed by Mutect
#                                              total.passed.pon.nsamp$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
#                                              total.passed.pon.nsamp$max.over.0.02 & # At least 2% VAF
#                                                (total.passed.pon.nsamp$min.under.0.35 | # Under 35% VAF 
#                                                   (!total.passed.pon.nsamp$min.under.0.35 & total.passed.pon.nsamp$CHROM.POS.REF.ALT %in% vars$key)) & # or if over 35% is hotspot
#                                              !total.passed.pon.nsamp$Vardict_PON_2AT2_percent & # Passed PON2
#                                              !total.passed.pon.nsamp$Mutect2_PON_2AT2_percent),] # Passed PON2
# paste(transplant.ch_pd2$CHROM.POS.REF.ALT, transplant.ch_pd2$eid,sep = ":")[!paste(transplant.ch_pd2$CHROM.POS.REF.ALT, transplant.ch_pd2$eid,sep = ":") %in% paste(transplant.ch_pd2.2$CHROM.POS.REF.ALT, transplant.ch_pd2.2$eid,sep = ":")]
# sum(paste(transplant.ch_pd2.2$CHROM.POS.REF.ALT, transplant.ch_pd2.2$eid,sep = ":") %in% paste(transplant.ch_pd2.for.cosmic$CHROM.POS.REF.ALT, transplant.ch_pd2.for.cosmic$eid,sep = ":"))

# transplant.complex <- total.passed.pon.nsamp[total.passed.pon.nsamp$Vardict_calpos_best != "",]
# sum(total.passed.pon.nsamp$Mutect2_CALLER)
get_oncokb = function(mut) {
  
  apiKey = "a83627cd-47e4-4be0-82dd-8f4cc9e4d6d0"
  request_url = paste0("https://oncokb.org//api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
                       mut['SYMBOL_VEP'], "&alteration=", mut['AAchange'], "&consequence=", mut['VariantClass'])
  
  res=httr::content(httr::GET(request_url, httr::add_headers(Authorization = paste("Bearer", apiKey))))
  return(res$oncogenic)
}
h = curl::new_handle()
curl::handle_setopt(h, http_version = 2)
httr::set_config(httr::config(http_version = 0))
cl = parallel::makeCluster(8)
transplant.ch_pd2$oncoKB <- parallel::parApply(cl, transplant.ch_pd2, 1, get_oncokb)
hotspots$oncoKB <- parallel::parApply(cl, hotspots, 1, get_oncokb)
parallel::stopCluster(cl)
transplant.ch_pd2$isOncogenic <- transplant.ch_pd2$oncoKB=="Oncogenic" | transplant.ch_pd2$oncoKB=="Likely Oncogenic"
hotspots$isOncogenic <- hotspots$oncoKB=="Oncogenic" | hotspots$oncoKB=="Likely Oncogenic"

###create list of tumor suppressor genes####
## file 1
gene_census = read.table(args$TSG_file, sep="\t", header=T)
## file 2
oncoKB_curated = read.table(args$oncoKB_curated, sep="\t", header=T)
oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
transplant.ch_pd2$isTSG <- transplant.ch_pd2$SYMBOL_VEP %in% TSG
hotspots$isTSG <- hotspots$SYMBOL_VEP %in% TSG

ct.heme <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.high.heme.tsv", 
                      sep = "\t", quote = "", header = T)
ct.myeloid <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.high.myeloid.tsv", 
                         sep = "\t", quote = "", header = T)
ct <- rbind(ct.heme,ct.myeloid) %>% group_by(GENOMIC_MUTATION_ID) %>% mutate(n=n()) %>% distinct()
ct$gene <- gsub("_.*", "", ct$Gene_HGVSp_VEP)
ct$aa.pos <- gsub(".*_", "", ct$Gene_HGVSp_VEP)
ct$aa.pos <- as.numeric(str_extract(ct$aa.pos, "\\d+"))
transplant.ch_pd2$aa.pos <- as.numeric(str_extract(gsub(".*_", "",transplant.ch_pd2$gene_loci_vep), "\\d+"))
hotspots$aa.pos <- as.numeric(str_extract(gsub(".*_", "",hotspots$gene_loci_vep), "\\d+"))
colnames(ct)[2] <- "key"
ct$CHROM.POS <- unlist(lapply(ct$key, function(x) paste(str_split(x,":")[[1]][1],str_split(x,":")[[1]][2],sep = ":")))
ct$GENE.AA.POS <- with(ct, paste(gene, aa.pos, sep=":"))

library(jsonlite)
transplant.ch_pd2$near.cosmic.HS <- apply(transplant.ch_pd2[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange")], 1, function(x) {
  p = c(-3:-1,1:9)
  n = c(-9:-1,1:9)
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% ct$GENE.AA.POS
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% ct$CHROM.POS
  if (any(any.in.p)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p], sep="|"))
    )
  } else if (any(any.in.n)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$CHROM.POS[ct$CHROM.POS %in% vector.n], ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n], sep="|"))
    )
  } else {
    return()
  }
})
transplant.ch_pd2$near.cosmic.HS <- sapply(transplant.ch_pd2$near.cosmic.HS, toJSON)
transplant.ch_pd2$near.cosmic.HS <- ifelse(transplant.ch_pd2$near.cosmic.HS == "{}","",transplant.ch_pd2$near.cosmic.HS)

hotspots$near.cosmic.HS <- apply(hotspots[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange")], 1, function(x) {
  p = c(-3:-1,1:9)
  n = c(-9:-1,1:9)
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% ct$GENE.AA.POS
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% ct$CHROM.POS
  if (any(any.in.p)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p], sep="|"))
    )
  } else if (any(any.in.n)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$CHROM.POS[ct$CHROM.POS %in% vector.n], ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n], sep="|"))
    )
  } else {
    return()
  }
})
hotspots$near.cosmic.HS <- sapply(hotspots$near.cosmic.HS, toJSON)
hotspots$near.cosmic.HS <- ifelse(hotspots$near.cosmic.HS == "{}","",hotspots$near.cosmic.HS)





# transplant.ch_pd2$M2_SB_Percent <- round(transplant.ch_pd2$Mutect2_gt_alt_fwd/(
#   transplant.ch_pd2$Mutect2_gt_alt_fwd+transplant.ch_pd2$Mutect2_gt_alt_rev
#   ),2)
# transplant.ch_pd2$`M2_SB_Alt_Forw,Rev` <- paste(transplant.ch_pd2$Mutect2_gt_alt_fwd,transplant.ch_pd2$Mutect2_gt_alt_rev, sep = ",")
# transplant.ch_pd2$CHROM.POS.REF.ALT <- with(transplant.ch_pd2, paste(CHROM,POS,REF,ALT, sep = ":"))
# transplant.ch_pd2$isHotSpot <- transplant.ch_pd2$CHROM.POS.REF.ALT %in% vars$key
# transplant.ch_pd2 <- transplant.ch_pd2 %>% left_join(vars %>% dplyr::select(key,source), by = c("CHROM.POS.REF.ALT"="key"))
# 
# alex <- rbind(read.table("/Users/brian/Bolton/UKBB/results/12/update/transplant_match/match/combined.UKBB.12.alex.complex.tsv",
#                    header = T, sep = "\t", comment.char = "", quote=""),
#               read.table("/Users/brian/Bolton/UKBB/results/12/update/transplant/combined.UKBB.12.alex.complex.tsv",
#                          header = T, sep = "\t", comment.char = "", quote="")
# )
# test.alex <- alex %>% filter(Vardict_gt_AF>=.02 & Vardict_gt_AF<=.35 & Vardict_PASS)
# test.alex <- test.alex[,c(1:4,46:51,6:45)]
# apply(transplant.ch_pd2[,c("Vardict_gt_AF","Mutect2_gt_AF","Mutect2_gt_AD_alt","Mutect2_gt_AD_ref")], 1, function(x) {
#   c(x[["Vardict_gt_AF"]],x[["Mutect2_gt_AF"]],x[["Mutect2_gt_AD_alt"]]/(x[["Mutect2_gt_AD_ref"]]+x[["Mutect2_gt_AD_alt"]]))
# })
# transplant.ch_pd2$VAF_list <- apply(transplant.ch_pd2[,c("Vardict_gt_AF","Mutect2_gt_AF","Mutect2_gt_AD_alt","Mutect2_gt_AD_ref")], 1, function(x) {
# list(x[["Vardict_gt_AF"]],x[["Mutect2_gt_AF"]],x[["Mutect2_gt_AD_alt"]]/(x[["Mutect2_gt_AD_ref"]]+x[["Mutect2_gt_AD_alt"]]))
# })
# transplant.ch_pd2$Vardict_VAF <- transplant.ch_pd2$Vardict_gt_AF
# transplant.ch_pd2$M1_VAF <- transplant.ch_pd2$Mutect2_gt_AF
# transplant.ch_pd2$M2_VAF <- transplant.ch_pd2$Mutect2_gt_AD_alt/(transplant.ch_pd2$Mutect2_gt_AD_alt+transplant.ch_pd2$Mutect2_gt_AD_ref)
#   
# test10 <- transplant.ch_pd2 %>% filter(passed_everything_mutect=="True" & ch_my_pd.y > 0)
# test10 <- transplant.ch_pd2 %>% filter(passed_everything_mutect=="True" & ch_my_pd.y > 0 & gnomAD_MAX.lessStringent.0007=="True")
# test10 <- transplant.ch_pd2 %>% filter(passed_everything_mutect=="True" & ch_my_pd.y > 0 & gnomAD_MAX.Stringent.0007=="True")
# 
# 
# 
# transplant.ch_pd2$Vardict_CALLER <- as.numeric(transplant.ch_pd2$Vardict_CALLER)
# transplant.ch_pd2$Mutect2_PON_2AT2_percent <- as.numeric(transplant.ch_pd2$Mutect2_PON_2AT2_percent)
# transplant.ch_pd2$alt_strand_counts_min_1_caller_only <- as.logical(transplant.ch_pd2$alt_strand_counts_min_1_caller_only)
# transplant.ch_pd2$max.over.0.02 <- as.logical(transplant.ch_pd2$max.over.0.02)
# transplant.ch_pd2$min.under.0.35 <- as.logical(transplant.ch_pd2$min.under.0.35)
# transplant.ch_pd2$max.over.0.02 <- as.logical(transplant.ch_pd2$max.over.0.02)
# transplant.ch_pd2$passed_everything <- as.logical(transplant.ch_pd2$passed_everything)
# transplant.ch_pd2$passed_everything_mutect <- as.logical(transplant.ch_pd2$passed_everything_mutect)
# transplant.ch_pd2$alt_strand_counts_min_2_callers <- as.logical(transplant.ch_pd2$alt_strand_counts_min_2_callers)
# transplant.ch_pd2$passed <- as.logical(transplant.ch_pd2$passed)
# transplant.ch_pd2$complexity_filters <- as.logical(transplant.ch_pd2$complexity_filters)
# transplant.ch_pd2$gnomAD_MAX.Stringent.0007 <- as.logical(transplant.ch_pd2$gnomAD_MAX.Stringent.0007)
# transplant.ch_pd2$gnomAD_MAX.lessStringent.0007 <- as.logical(transplant.ch_pd2$gnomAD_MAX.lessStringent.0007)
# 
# final.passed$passed_everything <- (!final.passed$Vardict_PON_2AT2_percent &
#                                      !final.passed$Mutect2_PON_2AT2_percent &
#                                      final.passed$alt_strand_counts_min_2_callers &
#                                      final.passed$min.under.0.35 &
#                                      final.passed$max.over.0.02 &
#                                      final.passed$passed &
#                                      final.passed$complexity_filters &
#                                      final.passed$PON_FISHER <= bf.correction)
# 
# transplant.ch_pd2_single <- transplant.ch_pd2[1,]
# 
# 
# 
# with(transplant.ch_pd2_single,
#      list(
#        paste("Vardict_PON2:", Vardict_PON_2AT2_percent),
#        paste("Mutect2_PON2:", Mutect2_PON_2AT2_percent),
#        paste("SB both Callers:", as.numeric(!alt_strand_counts_min_2_callers)),
#        paste("min.under.0.35:",as.numeric(!min.under.0.35)),
#        paste("max.over.0.02:",as.numeric(!max.over.0.02)),
#        paste("Complexity filters:",as.numeric(!complexity_filters)),
#        paste("Passed 2 Callers:",as.numeric(!passed)),
#        paste("PoN:",as.numeric(!(PON_FISHER <= bf.correction)))
#      ))
# 
# transplant.ch_pd2 %>% 
#   rowwise() %>%
#   dplyr::mutate(failure_reasons = list(
#     paste("Vardict_PON2:", Vardict_PON_2AT2_percent),
#     paste("Mutect2_PON2:", Mutect2_PON_2AT2_percent),
#     paste("SB both Callers:", as.numeric(!alt_strand_counts_min_2_callers)),
#     paste("min.under.0.35:",as.numeric(!min.under.0.35)),
#     paste("max.over.0.02:",as.numeric(!max.over.0.02)),
#     paste("Complexity filters:",as.numeric(!complexity_filters)),
#     paste("Passed 2 Callers:",as.numeric(!passed)),
#     paste("PoN:",as.numeric(!(PON_FISHER <= bf.correction)))
#   ))





##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

# transplant.ch_pd2 <- read.table("/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB3.tsv",
#                                 header = T, sep = "\t", comment.char = "", quote="")
# 
# transplant.ch_pd2.2 <- dplyr::left_join(transplant.ch_pd2, vars %>% dplyr::select(CHROM,POS,REF,ALT,n),
#                  by = c("CHROM"="CHROM", "POS"="POS", 
#                         "REF"="REF", "ALT"="ALT"))
# write.table(transplant.ch_pd2.2, "/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB.2.tsv", row.names = F,
#            col.names = T, sep = "\t", quote=F)
# # gnomAD 
# gnomAD.col <- grep("^gnomAD_.*_VEP", colnames(transplant.ch_pd2))
# 

# 
# getVAFs <- function(x) {
#   mutect_VAF <- grep("Mutect2_gt_AF", colnames(x))
#   vardict_VAF <- grep("Vardict_gt_AF", colnames(x))
#   VAFS <- c(mutect_VAF, vardict_VAF)
#   VAFS
# }
# transplant.ch_pd2$max.under.0.35 <- apply(transplant.ch_pd2[,getVAFs(transplant.ch_pd2)],1,function(x) max(x, na.rm = T) < 0.35)



#####################################
bf.correction <- 0.05/38793002
# transplant.ch_pd2 <- test
# fillna <- function(x, r) {
#   x[is.na(x)] = r
#   return(x)
# }

transplant.ch_pd2$n.HGVSp<- fillna(transplant.ch_pd2$n.HGVSp, 0)
transplant.ch_pd2$n.HGVSc<- fillna(transplant.ch_pd2$n.HGVSc, 0)
transplant.ch_pd2$isBoltonBickHotspot <- transplant.ch_pd2$n.loci.vep >=5 | transplant.ch_pd2$n.HGVSp | transplant.ch_pd2$n.HGVSc
# colnames(transplant.ch_pd2)[colnames(transplant.ch_pd2)=="isHotSpot"] <- "isBoltonBickHotspot"
transplant.ch_pd2$Vardict_PASS <- ifelse(transplant.ch_pd2$Vardict_FILTER == "PASS" | transplant.ch_pd2$Vardict_FILTER == "BCBIO",1,0)
transplant.ch_pd2$Vardict_PASS <- fillna(transplant.ch_pd2$Vardict_PASS, 0)
transplant.ch_pd2$PON_FISHER <- fillna(transplant.ch_pd2$PON_FISHER, 0)
transplant.ch_pd2$failures2 <- apply(transplant.ch_pd2[,c("Vardict_PASS","Vardict_PON_2AT2_percent","Mutect2_PON_2AT2_percent",
                                                                    "alt_strand_counts_min_2_callers","max.under.0.35","max.over.0.02",
                                                                    "complexity_filters","passed","PON_FISHER","isHotSpot",
                                                                    "Mutect2_CALLER")], 
      1, function (x) {
        reasons <- c("Vardict_PON2:", "Mutect2_PON2:","Passed Both Callers:","SB both Callers:","max.under 35% & not HS:",
                     "max.over 2%:","Complexity filters:","PoN:")
        values <- c(x[["Vardict_PON_2AT2_percent"]], x[["Mutect2_PON_2AT2_percent"]], as.numeric(!x[["passed"]]),
                    as.numeric(!x[["alt_strand_counts_min_2_callers"]]), as.numeric(!x[["max.under.0.35"]]), 
                    as.numeric(!x[["max.over.0.02"]]), as.numeric(!x[["complexity_filters"]]), 
                    as.numeric(!(x[["PON_FISHER"]] <= bf.correction)))
        
        failures = list()
        for (i in 1:length(values)) {
          if (values[i]==1) {
            if (reasons[i]=="max.under 35% & not HS:") {
              if (!x["isHotSpot"]) {
                failures <- append(failures,paste(reasons[i],values[i]))
              }
            } else if (reasons[i]=="Passed Both Callers:") {
              if (!x["Vardict_PASS"]) {
                failures <- append(failures,paste(reasons[i],values[i]))
              }
            }
            else {
              failures <- append(failures,paste(reasons[i],values[i]))
            }
          }
        }
        failures.out <- paste(failures, collapse = ",  ")
        if (!x[["Mutect2_CALLER"]]) {
          failures.out <- paste0(failures.out, ", Vardict Complex (Not called by Mutect)")
        }
        return(failures.out)
      }
)

hotspots$n.HGVSp<- fillna(hotspots$n.HGVSp, 0)
hotspots$n.HGVSc<- fillna(hotspots$n.HGVSc, 0)
hotspots$isBoltonBickHotspot <- hotspots$n.loci.vep >=5 | hotspots$n.HGVSp | hotspots$n.HGVSc
hotspots$Vardict_PASS <- ifelse(hotspots$Vardict_FILTER == "PASS" | hotspots$Vardict_FILTER == "BCBIO",1,0)
hotspots$Vardict_PASS <- fillna(hotspots$Vardict_PASS, 0)
hotspots$PON_FISHER <- fillna(hotspots$PON_FISHER, 0)
ho
hotspots$failures2 <- apply(hotspots[,c("Vardict_PASS","Vardict_PON_2AT2_percent","Mutect2_PON_2AT2_percent",
                                                          "alt_strand_counts_min_2_callers","max.under.0.35","max.over.0.02",
                                                          "complexity_filters","passed","PON_FISHER","isBoltonBickHotspot",
                                                          "Mutect2_CALLER")], 
                                     1, function (x) {
                                       reasons <- c("Vardict_PON2:", "Mutect2_PON2:","Passed Both Callers:","SB both Callers:","max.under 35% & not BoltonBickHS:",
                                                    "max.over 2%:","Complexity filters:","PoN:")
                                       values <- c(x[["Vardict_PON_2AT2_percent"]], x[["Mutect2_PON_2AT2_percent"]], as.numeric(!x[["passed"]]),
                                                   as.numeric(!x[["alt_strand_counts_min_2_callers"]]), as.numeric(!x[["max.under.0.35"]]), 
                                                   as.numeric(!x[["max.over.0.02"]]), as.numeric(!x[["complexity_filters"]]), 
                                                   as.numeric(!(x[["PON_FISHER"]] <= bf.correction)))
                                       
                                       failures = list()
                                       for (i in 1:length(values)) {
                                         if (values[i]==1) {
                                           if (reasons[i]=="max.under 35% & not BoltonBickHS:") {
                                             if (!x["isBoltonBickHotspot"]) {
                                               failures <- append(failures,paste(reasons[i],values[i]))
                                             }
                                           } else if (reasons[i]=="Passed Both Callers:") {
                                             failures <- append(failures,paste(reasons[i],values[i]))
                                           }
                                           else {
                                             failures <- append(failures,paste(reasons[i],values[i]))
                                           }
                                         }
                                       }
                                       failures.out <- paste(failures, collapse = ",  ")
                                       if (!x[["Mutect2_CALLER"]]) {
                                         failures.out <- paste0(failures.out, ", Vardict Complex (Not called by Mutect)")
                                       }
                                       return(failures.out)
                                     }
)

# apply(transplant.ch_pd2[2,c("Vardict_FILTER","Vardict_PON_2AT2_percent",
#                                                           "Mutect2_PON_2AT2_percent","alt_strand_counts_min_2_callers",
#                                                           "max.under.0.35","max.over.0.02",
#                                                           "complexity_filters","passed","PON_FISHER","isHotSpot")], 1,
#                                      function(x) {
#                                        x[2:length(x)] <- as.numeric(x[2:length(x)])
#                                        print(sapply(x, class))
#                                        print(x)
#                                        print(x[["passed"]])
#                                        #print(as.numeric(!x[["passed"]]))
#                                      }
# )
                                       
# GnomAD Filter Less Stringent (All 3 Groups had a Max > .0007)
# GnomAD Filter More Stringent (Any Group had a Max > .0007)
# transplant.ch_pd2$failuresGnomAD <- apply(transplant.ch_pd2[,c("gnomAD_MAX.Stringent.0007","gnomAD_MAX.lessStringent.0007")], 
#                                     1, function (x) {
#                                       reasons <- c("GnomAD Filter More Stringent:", 
#                                                    "GnomAD Filter Less Stringent:")
#                                       values <- c(as.numeric(!x[["gnomAD_MAX.Stringent.0007"]]), as.numeric(!x[["gnomAD_MAX.lessStringent.0007"]]))
#                                       
#                                       failures = list()
#                                       for (i in 1:length(values)) {
#                                         if (values[i]==1) {
#                                           failures <- append(failures,paste(reasons[i],values[i]))
#                                         }
#                                       }
#                                       paste(failures, collapse = ",  ")
#                                     })

transplant.ch_pd2$failuresGnomAD2 <- apply(transplant.ch_pd2[,c("gnomAD_MAX.Stringent.005","gnomAD_MAX.lessStringent.005")], 
                                          1, function (x) {
                                            reasons <- c("GnomAD Filter More Stringent:", 
                                                         "GnomAD Filter Less Stringent:")
                                            values <- c(as.numeric(!x[["gnomAD_MAX.Stringent.005"]]), as.numeric(!x[["gnomAD_MAX.lessStringent.005"]]))
                                            
                                            failures = list()
                                            for (i in 1:length(values)) {
                                              if (values[i]==1) {
                                                failures <- append(failures,paste(reasons[i],values[i]))
                                              }
                                            }
                                            paste(failures, collapse = ",  ")
                                          })

hotspots$failuresGnomAD2 <- apply(hotspots[,c("gnomAD_MAX.Stringent.005","gnomAD_MAX.lessStringent.005")], 
                                           1, function (x) {
                                             reasons <- c("GnomAD Filter More Stringent:", 
                                                          "GnomAD Filter Less Stringent:")
                                             values <- c(as.numeric(!x[["gnomAD_MAX.Stringent.005"]]), as.numeric(!x[["gnomAD_MAX.lessStringent.005"]]))
                                             
                                             failures = list()
                                             for (i in 1:length(values)) {
                                               if (values[i]==1) {
                                                 failures <- append(failures,paste(reasons[i],values[i]))
                                               }
                                             }
                                             paste(failures, collapse = ",  ")
                                           })

t1 <- readxl::read_excel("~/Bolton/UKBB/ch_pd2_transplant_KB_BW_FINAL_reannotate.xlsx")
t1$key <- with(t1, paste(CHROM,POS,REF,ALT,SAMPLE,sep=":"))
transplant.ch_pd2$samplekey <- paste(transplant.ch_pd2$key, transplant.ch_pd2$SAMPLE, sep = ":")
transplant.ch_pd2.already.have <- transplant.ch_pd2$samplekey %in% t1$key
sum(as.numeric(t1$ch_pd_final[(t1$key %in% transplant.ch_pd2$samplekey)]), na.rm = T)
sum(!transplant.ch_pd2$samplekey %in% t1$key)
new.transplant <- transplant.ch_pd2[!transplant.ch_pd2$samplekey %in% t1$key,]
write.table(new.transplant, "/Users/brian/Bolton/UKBB/new.transplant.passed.mutect.tsv", sep = "\t", quote = F, row.names = F)
write.table(hotspots, "/Users/brian/Bolton/UKBB/hotspots.failed.mutect.tsv", sep = "\t", quote = F, row.names = F)

# COSMIChg38 <- read.table("~/Bolton/data/hg38_cosmic78_parsed.sorted.txt", sep = "\t", quote = "", header = T)
# 
# COSMIChg38.high.heme <- COSMIChg38[COSMIChg38$heme_cosmic_count>=10,]

# ct <- COSMIChg38 %>%
#   # filter(Chrom=="chr21", Start==34880643) %>%
#   group_by(Chrom,Start,End, Ref, Alt) %>%
#   mutate(n=sum(heme_cosmic_count))
# ct2 <- ct %>% filter(n >= 10)
# ct.heme <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.high.heme.tsv", 
#                       sep = "\t", quote = "", header = T)
# ct.myeloid <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.high.myeloid.tsv", 
#                          sep = "\t", quote = "", header = T)
# ct <- rbind(ct.heme,ct.myeloid) %>% group_by(GENOMIC_MUTATION_ID) %>% mutate(n=n()) %>% distinct()
# ct$gene <- gsub("_.*", "", ct$Gene_HGVSp_VEP)
# ct$aa.pos <- gsub(".*_", "", ct$Gene_HGVSp_VEP)
# ct$aa.pos <- as.numeric(str_extract(ct$aa.pos, "\\d+"))
# transplant.ch_pd2$aa.pos <- as.numeric(str_extract(gsub(".*_", "",transplant.ch_pd2$gene_loci_vep), "\\d+"))
# colnames(ct)[2] <- "key"
# ct$CHROM.POS <- unlist(lapply(ct$key, function(x) paste(str_split(x,":")[[1]][1],str_split(x,":")[[1]][2],sep = ":")))
# ct$GENE.AA.POS <- with(ct, paste(gene, aa.pos, sep=":"))
# 
# library(jsonlite)
# transplant.ch_pd2$near.cosmic.HS <- apply(transplant.ch_pd2[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange")], 1, function(x) {
#   p = c(-3:-1,1:9)
#   n = c(-9:-1,1:9)
#   prot = p + as.integer(x["aa.pos"])
#   vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
#   any.in.p <- vector.p %in% ct$GENE.AA.POS
#   nuc = n + as.integer(x["POS"])
#   vector.n = paste(x["CHROM"], nuc ,sep = ":")
#   any.in.n <- vector.n %in% ct$CHROM.POS
#   if (any(any.in.p)) {
#     # print(paste("Yes",x))
#     return(
#       c(x[["gene_aachange"]],
#         paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p], sep="|"))
#     )
#   } else if (any(any.in.n)) {
#     # print(paste("Yes",x))
#     return(
#       c(x[["gene_aachange"]],
#         paste(ct$CHROM.POS[ct$CHROM.POS %in% vector.n], ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n], sep="|"))
#       )
#   } else {
#     return()
#   }
# })
# transplant.ch_pd2$near.cosmic.HS <- sapply(transplant.ch_pd2$near.cosmic.HS, toJSON)
# transplant.ch_pd2$near.cosmic.HS <- ifelse(transplant.ch_pd2$near.cosmic.HS == "{}","",transplant.ch_pd2$near.cosmic.HS)
# write.table(transplant.ch_pd2, "/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB3.tsv", sep = "\t", quote = F, row.names = F)

# transplant.ch_pd2 <- transplant.ch_pd2[with(transplant.ch_pd2,order(CHROM,POS,eid)),]
# row.names(transplant.ch_pd2) <- with(transplant.ch_pd2, paste(CHROM.POS.REF.ALT,eid, sep = ":"))
# library(readxl)
# xl <- read.table("/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB_BW_FINAL2.txt", sep = "\t", header = T, quote = "", comment.char = "")
# with(xl, paste(CHROM,POS,REF,ALT,eid, sep = ":"))
# blah <- transplant.ch_pd2[with(xl, paste(CHROM,POS,REF,ALT,eid, sep = ":")),]
# 
# 
# 
# xl %>% group_by(SYMBOL_VEP, artifact) %>% summarise(n=n())
