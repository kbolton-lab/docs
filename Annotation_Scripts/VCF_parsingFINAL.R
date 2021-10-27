# VCF parsing
# 03/25/2021
# Yizhe Song
# Bolton lab
# Total ba
# VARIANT_CLASS_VEP

rm(list = ls())
if (!require('vcfR')) install.packages('vcfR')
if (!require('tidyverse')) install.packages('tidyverse')

library(vcfR)
library(tidyverse)

getVAFs <- function(x) {
  mutect_VAF <- grep("Mutect2_gt_AF", colnames(x))
  vardict_VAF <- grep("Vardict_gt_AF", colnames(x))
  varscan_VAF <- grep("Varscan2_gt_FREQ", colnames(x))
  VAFS <- c(mutect_VAF, vardict_VAF, varscan_VAF)
  VAFS
}

#mutect2_vcf_file <- "~/Bolton/final_vcfs/final.mutect2.vcf"
#mutect2_vcf_file <- "~/Bolton/final_vcfs/CORRECTED.final.mutect2.vcf"
mutect2_vcf_file <- "~/Bolton/final_vcfs/REALLY.CORRECTED.final.mutect2.vcf"
varscan2_vcf_file <- "~/Bolton/final_vcfs/final.varscan2.vcf"
vardict_vcf_file <- "~/Bolton/final_vcfs/final.vardict.vcf"


###################################################################
#####            extract all fields      ###
###################################################################
mutect2 <- vcfR2tidy(read.vcfR(mutect2_vcf_file, verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
mutect2.info <- mutect2$meta[mutect2$meta$Tag=="INFO",]$ID[1:28]
# mutect2.info <- c("AC","AF","AN","AS_FilterStatus","AS_SB_TABLE","AS_UNIQ_ALT_READ_COUNT","CONTQ","DP","ECNT","GERMQ",
#                   "MBQ","MFRL","MMQ","MPOS","NALOD","NCount","NLOD","OCM","PON","POPAF","ROQ","RPA","RU","SEQQ","STR",
#                   "STRANDQ","STRQ","TLOD" )
mutect.df <- mutect2[["dat"]]

varscan2 <- vcfR2tidy(read.vcfR(varscan2_vcf_file, verbose = FALSE),
                      single_frame = TRUE,
                      info_types = TRUE,
                      format_types = TRUE)
varscan2.info <- varscan2$meta[varscan2$meta$Tag=="INFO",]$ID[1:8]
# varscan2.info <- c("DP","SOMATIC","SS","SSC","GPV","SPV","AC","AN")
varscan.df <- varscan2[["dat"]]


vardict <- vcfR2tidy(read.vcfR(vardict_vcf_file, verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
vardict.info <- vardict$meta[vardict$meta$Tag=="INFO",]$ID[2:18]
# vardict.info <- c("TYPE","DP","END","VD","AF","SHIFT3","MSI","MSILEN","SSF","SOR",
#                   "LSEQ","RSEQ","STATUS","P0.01Likely","InDelLikely","AC","AN")
vardict.df <- vardict[["dat"]]

#####################################################################################
##### find duplicate
dup.mutect <- mutect.df[duplicated(mutect.df)==TRUE|duplicated(mutect.df, fromLast = TRUE),]
dup.varscan <- varscan.df[duplicated(varscan.df)==TRUE|duplicated(varscan.df, fromLast = TRUE),]
dup.vardict <- vardict.df[duplicated(vardict.df)==TRUE|duplicated(vardict.df, fromLast = TRUE),]
# remove duplicate
mutect.df <- mutect.df[!duplicated(mutect.df),]
varscan.df <- varscan.df[!duplicated(varscan.df),]
vardict.df <- vardict.df[!duplicated(vardict.df),]

#Split the FORMAT/AD for Varscan and Vardict only
# gt_AD = allelelic depths for ref then for alt sep by a comma
mutect.df <-mutect.df %>% separate(gt_AD, c("gt_AD_ref","gt_AD_alt"), 
                        sep=",",extra = "merge", fill = "right")
vardict.df <-vardict.df %>% separate(gt_AD, c("gt_AD_ref","gt_AD_alt"), 
                                               sep=",",extra = "merge", fill = "right")
#Rename it for Varscan
colnames(varscan.df)[colnames(varscan.df)=="gt_RD"] = "gt_AD_ref"
colnames(varscan.df)[colnames(varscan.df)=="gt_AD"] = "gt_AD_alt"

#covert AD from character to numeric
mutect.df[,c("gt_AD_ref","gt_AD_alt")] <- lapply(mutect.df[,c("gt_AD_ref","gt_AD_alt")], 
                                                       function(x) as.numeric(as.character(x)))
varscan.df[,c("gt_AD_ref","gt_AD_alt")] <- lapply(varscan.df[,c("gt_AD_ref","gt_AD_alt")], 
                                                        function(x) as.numeric(as.character(x)))
vardict.df[,c("gt_AD_ref","gt_AD_alt")] <- lapply(vardict.df[,c("gt_AD_ref","gt_AD_alt")], 
                                                        function(x) as.numeric(as.character(x)))

## ran R scipt already
all(varscan.df$PON_FISHER<.05/38663445)
all(mutect.df$PON_FISHER<.05/38663445)
all(vardict.df$PON_FISHER<.05/38663445)


# split FILTER here
mutect.df <- mutect.df %>% separate_rows(FILTER,sep =";", convert = TRUE)
# prefix filter with Mutect2 such as Mutect2_strand_bias
mutect.df$FILTER <- paste0("Mutect2_",mutect.df$FILTER)
mutect.df$value <- 1
mutect.df <- mutect.df %>% 
  spread(FILTER,value, fill = 0) # 1= passed the filter,0= didn't pass the filter

varscan.df <- varscan.df %>% separate_rows(FILTER,sep =";", convert = TRUE)
varscan.df$FILTER <- paste0("Varscan2_",varscan.df$FILTER)
varscan.df$value <- 1
varscan.df <- varscan.df %>% 
  spread(FILTER,value, fill = 0) # 1= passed the filter,0= didn't pass the filter

vardict.df <- vardict.df %>% separate_rows(FILTER,sep =";", convert = TRUE)
vardict.df$FILTER <- paste0("Vardict_",vardict.df$FILTER)
vardict.df$value <- 1
vardict.df <- vardict.df %>% 
  spread(FILTER,value, fill = 0) # 1= passed the filter,0= didn't pass the filter

# Prefix all FORMAT/ (gt_) columns with caller such as Vardict_gt_AD_ref & Vardict_gt_AD_alt
colnames(mutect.df)[which(grepl("^gt_", colnames(mutect.df)))] <-
  paste0("Mutect2_", colnames(mutect.df)[which(grepl("^gt_", colnames(mutect.df)))])
# then move all FILTER & FORMAT columns to front after ALT column
mutect2.pfx <- which(grepl("Mutect2_", colnames(mutect.df)))
mutect.new.col.ord <- c(1:5,mutect2.pfx,6:(min(mutect2.pfx)-1))
length(mutect.new.col.ord) == dim(mutect.df)[2]
mutect2.final <- mutect.df[,mutect.new.col.ord] ## 88 cols

colnames(varscan.df)[which(grepl("^gt_", colnames(varscan.df)))] <-
  paste0("Varscan2_", colnames(varscan.df)[which(grepl("^gt_", colnames(varscan.df)))])
varscan2.pfx <- which(grepl("Varscan2_", colnames(varscan.df)))
varscan2.new.col.ord <- c(1:5,varscan2.pfx,6:(min(varscan2.pfx)-1))
length(varscan2.new.col.ord) == dim(varscan.df)[2]
varscan2.final <- varscan.df[,varscan2.new.col.ord] ## 60 cols

colnames(vardict.df)[which(grepl("^gt_", colnames(vardict.df)))] <-
  paste0("Vardict_", colnames(vardict.df)[which(grepl("^gt_", colnames(vardict.df)))])
vardict.pfx <- which(grepl("Vardict_", colnames(vardict.df)))
vardict.new.col.ord <- c(1:5,vardict.pfx,6:(min(vardict.pfx)-1))
length(vardict.new.col.ord) == dim(vardict.df)[2]
vardict.final <- vardict.df[,vardict.new.col.ord] ## 79 cols

#################################################################################################
# VEP CSQ
# using tidyr::separate to tidy CSQ column even further,sep="\\|"
string <-"Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|VAR_SYNONYMS|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames <- str_split(string, "\\|")[[1]]
# CSQ has 91 columns, maybe we should suffix cols with VEP?
mutect2.final <- mutect2.final %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
varscan2.final <- varscan2.final %>% separate(CSQ, paste0(CSQnames, "_VEP"),  sep="\\|", extra = "merge", fill = "right")
vardict.final <- vardict.final %>% separate(CSQ, paste0(CSQnames, "_VEP"),  sep="\\|", extra = "merge", fill = "right")


## remove columns where all NA
na.mutect.cols <- names(mutect2.final[,apply(mutect2.final, 2, function(x) all(is.na(x)))])
na.varscan.cols <- names(varscan2.final[,apply(varscan2.final, 2, function(x) all(is.na(x)))])
na.vardict.cols <- names(vardict.final[,apply(vardict.final, 2, function(x) all(is.na(x)))])

mutect2.final <- mutect2.final[,-which(names(mutect2.final) %in% c(na.mutect.cols))]
varscan2.final <- varscan2.final[,-which(names(varscan2.final) %in% c(na.varscan.cols))]
vardict.final <- vardict.final[,-which(names(vardict.final) %in% c(na.vardict.cols))]

mutect2_vcf_file <- "~/Bolton/final_vcfs/final.mutect2.vcf"
varscan2_vcf_file <- "~/Bolton/final_vcfs/final.varscan2.vcf"
vardict_vcf_file <- "~/Bolton/final_vcfs/final.vardict.vcf"
write.table(mutect2.final, "~/Bolton/final_vcfs/REALLY.CORRECTED.final.mutect2.tsv",
            row.names = F, quote = F, sep = "\t")
write.table(varscan2.final, "~/Bolton/final_vcfs/final.varscan2.tsv",
            row.names = F, quote = F, sep = "\t")
write.table(vardict.final, "~/Bolton/final_vcfs/final.vardict.tsv",
            row.names = F, quote = F, sep = "\t")
########################################################################################
## ch_pd stuff
library(dplyr)
library(readxl)
library(stringr)
library(rtracklayer)
library(GenomicRanges)

mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- read.table("~/Bolton/data/hg38_mut_full_long_filtered_KB_deid_2.tsv",
                                                    quote = "", header = T, sep = "\t")
write.table(mut_full_long_filtered_KB_deid.ch_pd2, "~/Bolton/data/hg19_mut_full_long_filtered_KB_deid_2.tsv", sep = "\t",
          quote=F, row.names = F)
# mut_full_long_filtered_KB_deid.ch_pd2 <- readxl::read_xlsx("~/Bolton/data/mut_full_long_filtered_KB_deid_2.xlsx")
# 
# mut_full_long.gr <- GRanges(seqnames = paste0("chr",mut_full_long_filtered_KB_deid.ch_pd2$Chrom),
#                             ranges = IRanges(start=mut_full_long_filtered_KB_deid.ch_pd2$Start, 
#                                              end=mut_full_long_filtered_KB_deid.ch_pd2$Start), 
#                             mcols = mut_full_long_filtered_KB_deid.ch_pd2[,3:101])
# chain = import.chain("~/Bolton/data/hg19ToHg38.over.chain")
# mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- as.data.frame(liftOver(mut_full_long.gr,chain))
# mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- mut_full_long_filtered_KB_deid.ch_pd2.hg38[,c(3:4,8:106)]
# colnames(mut_full_long_filtered_KB_deid.ch_pd2.hg38) <- colnames(mut_full_long_filtered_KB_deid.ch_pd2) 
# write.table(mut_full_long_filtered_KB_deid.ch_pd2.hg38, "mut_full_long_filtered_KB_deid.ch_pd2.hg38.tsv", sep = "\t",
#           quote=F, row.names = F)


bick.topmed <- read.table("~/Bolton/data/bick_topmed_variants.txt",
                          quote = "", header = T, sep = "\t")
bick.topmed$CHROM <- paste0("chr", bick.topmed$CHROM)

num=2
## combine to count only
# https://community.rstudio.com/t/summarise-max-but-keep-all-columns/52449/3
topmed.loci.n <- rbind(bick.topmed %>%
                         filter(!grepl("c\\.", AAchange)) %>%
                         mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                         group_by(Gene, loci) %>%
                         mutate(n=dplyr::n()) %>%
                         ungroup() %>%
                         filter(n >= num) %>%
                         #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                         dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n),
                       bick.topmed %>%
                         filter(grepl("c\\.", AAchange)) %>%
                         mutate(loci = str_extract(AAchange, "c\\.\\d+[\\-\\+]\\d[TCGA-]")) %>% # c.769+2T>G and -> is for insertion
                         group_by(Gene, loci) %>%
                         mutate(n=dplyr::n()) %>%
                         ungroup() %>%
                         filter(n >= num) %>%
                         #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                         dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n)
                       )

topmed.loci.n <- topmed.loci.n %>% mutate(gene_loci = paste(Gene,loci,sep = "_"))

topmed.mutation.1 <- bick.topmed %>% 
  group_by(Gene, AAchange) %>% 
  dplyr::count() %>%  
  filter(n >= 1) %>% 
  mutate(gene_aachange = paste(Gene,AAchange,sep = "_"))

colnames(mut_full_long_filtered_KB_deid.ch_pd2.hg38)[1:4] <- colnames(bick.topmed)[1:4]
kelly.loci.n <- rbind(mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>%
                        filter(!grepl("c\\.", AAchange)) %>%
                        mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                        group_by(Gene, loci) %>%
                        mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
                        ungroup() %>%
                        filter(n >= num) %>%
                        #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                        dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,VariantClass),
                      mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>%
                        filter(grepl("c\\.", AAchange)) %>%
                        mutate(loci = str_extract(AAchange, "c\\.\\d+[\\-\\+]\\d[TCGA-]")) %>% # c.769+2T>G and -> is for insertion
                        group_by(Gene, loci) %>%
                        mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
                        ungroup() %>%
                        filter(n >= num) %>%
                        #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                        dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,VariantClass)
)
test = mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>%
  filter(!grepl("c\\.", AAchange)) %>%
  mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
  group_by(Gene, loci) %>%
  mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
  ungroup() %>%
  filter(Gene == "POLD1") %>%
  dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,VariantClass)



kelly.loci.n <- kelly.loci.n %>% mutate(gene_loci = paste(Gene,loci,sep = "_"))
kelly.mutation.1 <- mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>% 
  group_by(Gene, AAchange) %>% 
  dplyr::count() %>%  
  filter(n >= 1) %>% 
  mutate(gene_aachange = paste(Gene,AAchange,sep = "_"))

vars.loci.n <- rbind(mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>%
                dplyr::select(CHROM,POS,REF,ALT,Gene,AAchange,ch_my_pd,ch_pd2) %>%
                mutate(source="kelly"),
              bick.topmed %>%
                dplyr::select(CHROM,POS,REF,ALT,Gene,AAchange) %>%
                mutate(ch_my_pd = NA,
                       ch_pd2 = NA,
                       source = "bick"))

vars <- rbind(vars.loci.n %>%
                filter(!grepl("c\\.", AAchange)) %>%
                mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                group_by(Gene, loci) %>%
                mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
                ungroup() %>%
                dplyr::filter(n >= num) %>%
                # dplyr::filter(is.na(loci)) %>%
                #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,ch_my_pd,ch_pd2,source),
              vars.loci.n %>%
                filter(grepl("c\\.", AAchange)) %>%
                mutate(loci = str_extract(AAchange, "c\\.\\d+[\\-\\+]\\d[TCGA-]")) %>% # c.769+2T>G and -> is for insertion
                group_by(Gene, loci) %>%
                mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
                ungroup() %>%
                filter(n >= num) %>%
                #dplyr::filter(is.na(loci)) %>%
                #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,ch_my_pd,ch_pd2,source)) 

vars <- vars[!duplicated(vars[,c("CHROM","POS","REF","ALT")]),]
library(VariantAnnotation)
gr <- GRanges(seqnames = vars$CHROM,
              ranges = IRanges(start=vars$POS,
                               end = vars$POS),
              REF = vars$REF,
              ALT = vars$ALT)
gr <- sortSeqlevels(gr)
gr <- sort(gr)
final.vcf <- VariantAnnotation::VCF(rowRanges = gr, fixed = DataFrame(REF=DNAStringSet(gr$REF),
                                                                      ALT=DNAStringSetList(as.list(gr$ALT)),
                                                                      QUAL = 1))
if (num==5) {
  VariantAnnotation::writeVcf(final.vcf, "~/Bolton/vars.loci5.vcf")
  system(command = "gsed -i '1 i\\##fileformat=VCFv4.2' ~/Bolton/vars.loci5.vcf", intern = T)
  system(command = "bgzip -f ~/Bolton/vars.loci5.vcf && tabix ~/Bolton/vars.loci5.vcf.gz", intern = T)
} else if (num==2) {
  VariantAnnotation::writeVcf(final.vcf, "~/Bolton/vars.loci2.vcf")
  system(command = "gsed -i '1 i\\##fileformat=VCFv4.2' ~/Bolton/vars.loci2.vcf", intern = T)
  system(command = "bgzip -f ~/Bolton/vars.loci2.vcf && tabix ~/Bolton/vars.loci2.vcf.gz", intern = T)
}
########################################################################################
########################################################################################


mutect2.final <- read.table("~/Bolton/final_vcfs/REALLY.CORRECTED.final.mutect2.tsv", sep = "\t",
                            header = T, quote = "")
varscan2.final <- read.table("~/Bolton/final_vcfs/final.varscan2.tsv", sep = "\t",
                            header = T, quote = "")
vardict.final <- read.table("~/Bolton/final_vcfs/final.vardict.tsv", sep = "\t",
                             header = T, quote = "")


vardict.final.driver <- left_join(vardict.final, vars,
                                           by = c("CHROM"="CHROM", "POS"="POS", 
                                                  "REF"="REF", "ALT"="ALT"))
colnames(vardict.final.driver)[colnames(vardict.final.driver) %in% vardict.info] <-
  paste0("Vardict_", colnames(vardict.final.driver)[colnames(vardict.final.driver) %in% vardict.info])

mutect2.final.driver <- left_join(mutect2.final, vars,
                                           by = c("CHROM"="CHROM", "POS"="POS", 
                                                  "REF"="REF", "ALT"="ALT"))
colnames(mutect2.final.driver)[colnames(mutect2.final.driver) %in% mutect2.info] <-
  paste0("Mutect2_", colnames(mutect2.final.driver)[colnames(mutect2.final.driver) %in% mutect2.info])

varscan2.final.driver <- left_join(varscan2.final, vars,
                           by = c("CHROM"="CHROM", "POS"="POS", 
                                  "REF"="REF", "ALT"="ALT"))
colnames(varscan2.final.driver)[colnames(varscan2.final.driver) %in% varscan2.info] <-
  paste0("Varscan2_", colnames(varscan2.final.driver)[colnames(varscan2.final.driver) %in% varscan2.info])


library(plyr)
library(tidyverse)


## Way 3 works
#https://stackoverflow.com/questions/16042380/merge-data-frames-and-overwrite-values
library(dplyr)
colnames(mutect2.final.driver)[which(colnames(mutect2.final.driver)=="CALLER")] <- "Mutect2_CALLER"
colnames(varscan2.final.driver)[which(colnames(varscan2.final.driver)=="CALLER")] <- "Varscan2_CALLER"
colnames(vardict.final.driver)[which(colnames(vardict.final.driver)=="CALLER")] <- "Vardict_CALLER"
mutect2.final.driver$Mutect2_CALLER <- 1
varscan2.final.driver$Varscan2_CALLER <- 1
vardict.final.driver$Vardict_CALLER <- 1


intersection <- intersect(colnames(vardict.final.driver), colnames(varscan2.final.driver))
intersection <- intersection[6:length(intersection)] ## first 5 are the grouped by columns
intersection.cols.x <- paste0(intersection, ".x")
intersection.cols.y <- paste0(intersection, ".y")
final <- full_join(vardict.final.driver, varscan2.final.driver,
                   by = c("CHROM"="CHROM", "POS"="POS", 
                          "REF"="REF", "ALT"="ALT", "SAMPLE"="SAMPLE"))
final[,intersection.cols.x][is.na(final[,intersection.cols.x])] <- final[,intersection.cols.y][is.na(final[,intersection.cols.x])]

colnames(final)[colnames(final) %in% intersection.cols.x] <- intersection
final <- final[,!(colnames(final) %in% intersection.cols.y)]
final$PON_RefDepth <- as.numeric(final$PON_RefDepth)
final$PON_AltCounts <- as.numeric(final$PON_AltCounts)



intersection <- intersect(colnames(final), colnames(mutect2.final.driver))
intersection <- intersection[6:length(intersection)]
intersection.cols.x <- paste0(intersection, ".x")
intersection.cols.y <- paste0(intersection, ".y")
final <- full_join(final, mutect2.final.driver,
                   by = c("CHROM"="CHROM", "POS"="POS", 
                          "REF"="REF", "ALT"="ALT", "SAMPLE"="SAMPLE"))
# Error: Can't combine `Allele_VEP.x` <character> and `PON_RefDepth.x` <integer>.
final[,intersection.cols.x][is.na(final[,intersection.cols.x])] <- final[,intersection.cols.y][is.na(final[,intersection.cols.x])]
colnames(final)[colnames(final) %in% intersection.cols.x] <- intersection
final <- final[,!(colnames(final) %in% intersection.cols.y)] ## 227 columns
## variants from Kelly and Bick
sum(!is.na(final$n)) # 80 for n>=1, 28 for n>=2

final[,c("PON_AltCounts", "PON_FISHER", "PON_RefDepth")] <- sapply(final[,c("PON_AltCounts", "PON_FISHER", "PON_RefDepth")], as.numeric)


## fill NA with 0
final[, c("Mutect2_CALLER", "Varscan2_CALLER", "Vardict_CALLER")][is.na(final[, c("Mutect2_CALLER", 
                                                                                  "Varscan2_CALLER", 
                                                                                  "Vardict_CALLER")])] <- 0
## fill NA with 0
final[, c("Mutect2_PASS", "Varscan2_PASS", "Vardict_PASS")][is.na(final[, c("Mutect2_PASS", 
                                                                                  "Varscan2_PASS", 
                                                                                  "Vardict_PASS")])] <- 0
## coding 
final.coding <- final[!is.na(final$Consequence_VEP),]
sum(!is.na(final.coding$n)) # 27 for n>=2


write.table(final.coding, "~/Bolton/final_vcfs/final.coding.tsv", 
            row.names = F, quote = F, sep = "\t")
final.coding <- read.table("~/Bolton/final_vcfs/final.coding.tsv", sep = "\t", quote="")

VEP.cols <- which(grepl("VEP$", colnames(final)))
final.noncoding <-  final[is.na(final$Consequence_VEP),VEP.cols]
all(sapply(final.noncoding, function(x) all(is.na(x))))

final.coding$MAX_AF.Stringent.0007 <- final.coding$MAX_AF_VEP < .0007

## JUST COPYING for gnomad MAX AF stringent, MAX_AF is other databases
final.coding.gnomad = final.coding

final.coding.gnomad.sorted <- final.coding.gnomad[with(final.coding.gnomad, order(CHROM, POS)), ]
final.coding.gnomad.sorted$Varscan2_gt_FREQ <- as.numeric(sub("%", "", final.coding.gnomad.sorted$Varscan2_gt_FREQ))/100

gnomad.col <- grep("gnomAD.*AF", colnames(final.coding.gnomad.sorted))
final.coding.gnomad.sorted[,gnomad.col] <- apply(final.coding.gnomad.sorted[,gnomad.col], 2, as.numeric)
final.coding.gnomad.sorted[,gnomad.col][is.na(final.coding.gnomad.sorted[,gnomad.col])] <- 0
final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP <- as.numeric(apply(final.coding.gnomad.sorted[,gnomad.col], 1, function(x){
  max(x)
}))
final.coding.gnomad.sorted$gnomAD_MAX.Stringent.0007 <- final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.0007
AF <- grep("MAX_AF_VEP", colnames(final.coding.gnomad.sorted))
final.coding.gnomad.sorted <- final.coding.gnomad.sorted[,c(1:111,113:223,AF,224:ncol(final.coding.gnomad.sorted))]


#########################################################################################################
### gnomad age filter
X <- split(final.coding.gnomad.sorted, final.coding.gnomad.sorted$CHROM)

for (df in ls(X)) {
  write.table(X[[df]][,1:4], paste0("~/Bolton/gnomad_chrom/", df, ".tsv"),
              sep = "\t", quote = F, col.names = T, row.names = F)
}


mutect_VAF <- grep("Mutect2_gt_AF", colnames(final.coding.gnomad.sorted))
vardict_VAF <- grep("Vardict_gt_AF", colnames(final.coding.gnomad.sorted))
varscan_VAF <- grep("Varscan2_gt_FREQ", colnames(final.coding.gnomad.sorted))
VAFS <- c(mutect_VAF, vardict_VAF, varscan_VAF)

all_ages <- c(rep(32.5, 1428), rep(37.5, 1606), rep(42.5, 1905), 
              rep(47.5, 3228), rep(52.5, 4577), rep(57.5, 3924), 
              rep(62.5, 3656), rep(67.5, 3194), rep(72.5, 2153), rep(77.5, 1283))
ages <- c(32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,77.5)
parsed <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1)

gnomad.age.keep <- function(x) {
  parsed = as.integer(unlist(strsplit(gsub("\\[|\\]| ", "", x["agebins"]), split = ",")))
  if (all(parsed == 0)) {
    #print("true")
    if (max(x[VAFS], na.rm = T) > 0.35) {
      # remove
      return(F)
    }
    else {
      #keep
      return(T)
    }
  }
  else {
    variant.ages <- rep(cbind(ages,parsed)[,1], cbind(ages,parsed)[,2])
    wc <- wilcox.test(variant.ages, all_ages, alternative="greater")
    ## VAF > 35% and Gnomad variant mean age == Gnomad population mean age
    if (max(x[VAFS], na.rm = T) > 0.35 & wc$p.value > .05) {
      return(F)
    }
    else {
      return(T)
    }
  }
}

#x <- "chr9"
library(parallel)
X.1 <- parallel::mclapply(ls(X), function(x) {
  chr <- read.table(paste0("~/Bolton/gnomad_chrom2/", x, ".agebins.tsv"),
                    header = T,
                    sep = "\t")
  if (dim(X[[x]])[1] == dim(chr)[1]) {
    if (all(X[[x]][,1:4] == chr[,1:4])) {
      message(paste(x, "matches"))
      X[[x]]$agebins <- chr$Age_Bin
      X[[x]]$gnomad.age.keep <- apply(X[[x]], 1, gnomad.age.keep)
    }
  } else {
    chr <- chr[!duplicated(chr),]
    test2 <- left_join(X[[x]], chr,
                        by = c("CHROM"="CHROM", "POS"="POS", 
                               "REF"="REF", "ALT"="ALT"))
    if (all(X[[x]][,1:4] == test2[,1:4])) {
      message(paste(x, "matches after join"))
      X[[x]] = test2
      colnames(X[[x]])[length(colnames(X[[x]]))] <- "agebins"
      X[[x]]$gnomad.age.keep <- apply(X[[x]], 1, gnomad.age.keep)
    }
  }
  X[[x]]
})
names(X.1) <- ls(X)

##confirm #rows are same for all chromosomes
all(sapply(ls(X), function(x) {
  dim(X[[x]])[1] == dim(X.1[[x]])[1]
}))

## combine chrs for age bins
final.coding.gnomad.sorted.gnomadAge <- dplyr::bind_rows(X.1)

hotspot <- grep("^n$", colnames(final.coding.gnomad.sorted.gnomadAge))
ch_pd2 <- grep("^ch_pd2$", colnames(final.coding.gnomad.sorted.gnomadAge))
ch_my_pd <- grep("^ch_my_pd$", colnames(final.coding.gnomad.sorted.gnomadAge))
final.coding.gnomad.sorted.gnomadAge <- final.coding.gnomad.sorted.gnomadAge[,c(1:150,
                                                                                154:ncol(final.coding.gnomad.sorted.gnomadAge),
                                                                                hotspot,ch_my_pd,ch_pd2)]
sum(!is.na(final.coding.gnomad.sorted.gnomadAge$n)) # 27 for n>=2


write.table(final.coding.gnomad.sorted.gnomadAge,
            "~/Bolton/final_vcfs/final.coding.gnomad.sorted.gnomadAge.tsv", 
            sep = "\t",
            quote = F, row.names = F)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################

library(dplyr)
library(tidyverse)
## coding variants but has columns to filter out advanced gnomad filters 
#"gnomAD_MAX.Stringent.0007" "gnomad.age.keep"
final.coding.gnomad.sorted.gnomadAge <- read.table("~/Bolton/final_vcfs/final.coding.gnomad.sorted.gnomadAge.tsv", 
                                                   sep = "\t", header = T, quote = "")


#blacklist <- read.table("~/Bolton/HIV/HIV.blacklist.bed", header = F, sep = "\t")
blacklist <- read.table("/Users/brian/ENCFF356LFX.bed", header = F, sep = "\t")
colnames(blacklist) <- c("CHROM", "START", "END")
blacklist$blacklist <- T
blackist.gr <- GRanges(blacklist)
HIV.gr <- GRanges(seqnames = final.coding.gnomad.sorted.gnomadAge$CHROM,
                  ranges = IRanges(start=final.coding.gnomad.sorted.gnomadAge$POS, 
                                   end=final.coding.gnomad.sorted.gnomadAge$POS), 
                  mcols = final.coding.gnomad.sorted.gnomadAge[,3:ncol(final.coding.gnomad.sorted.gnomadAge)])
test <- GenomicRanges::findOverlaps(blackist.gr, HIV.gr)
test <- plyranges::join_overlap_left(HIV.gr, blackist.gr)
columns <- colnames(final.coding.gnomad.sorted.gnomadAge)
final.coding.gnomad.sorted.gnomadAge <- as.data.frame(test)[,c(1,2,6:ncol(as.data.frame(test)))]
colnames(final.coding.gnomad.sorted.gnomadAge) <- c(columns, "blacklist")
final.coding.gnomad.sorted.gnomadAge$blacklist[is.na(final.coding.gnomad.sorted.gnomadAge$blacklist)] <- F

############################################################
############################################################

## should be 188916 rows
final.gnomad.pass <- final.coding.gnomad.sorted.gnomadAge[(final.coding.gnomad.sorted.gnomadAge$gnomAD_MAX.Stringent.0007 &   #################################################
                                                final.coding.gnomad.sorted.gnomadAge$gnomad.age.keep),]                       #################################################


final.coding.gnomad.sorted.gnomadAge$called <- 
  apply(final.coding.gnomad.sorted.gnomadAge[,c("Mutect2_CALLER","Varscan2_CALLER","Vardict_CALLER")],
        1,
        function(x) sum(x)>=2)
sum((final.coding.gnomad.sorted.gnomadAge %>%
      filter(#gnomAD_MAX.Stringent.0007==1 & 
               Mutect2_CALLER &
               gnomad.age.keep==1) %>%
      dplyr::select(called))$called) # 26555 before mutect2 samples fix, after 26551

final.called <- final.coding.gnomad.sorted.gnomadAge

library(tidyverse)
final.called <- final.called %>% separate(Vardict_gt_ALD, c("Vardict_gt_ALD_forw","Vardict_gt_ALD_rev"), 
                                          sep=",",extra = "merge", fill = "right")
# vardict.final <- vardict.final %>% separate(Vardict_gt_ALD, c("Vardict_gt_ALD_forw","Vardict_gt_ALD_rev"), 
#                                            sep=",",extra = "merge", fill = "right")
# apply(vardict.final[,c("Vardict_gt_ALD_forw","Vardict_gt_ALD_rev")],2,as.numeric)
# vardict.final$altsum <- apply(apply(vardict.final[,c("Vardict_gt_ALD_forw","Vardict_gt_ALD_rev")],2,as.numeric),1,sum)
#                                                              "Mutect2_gt_alt_fwd","Mutect2_gt_alt_rev"), 
#                                           sep=",",extra = "merge", fill = "right")
final.called <- final.called %>% separate(Mutect2_gt_SB, c("Mutect2_gt_ref_fwd","Mutect2_gt_ref_rev",
                                                           "Mutect2_gt_alt_fwd","Mutect2_gt_alt_rev"), 
                                          sep=",",extra = "merge", fill = "right")
final.called <- final.called %>% separate(Varscan2_gt_DP4, c("Varscan2_gt_ref_fwd","Varscan2_gt_ref_rev",
                                                             "Varscan2_gt_alt_fwd","Varscan2_gt_alt_rev"), 
                                          sep=",",extra = "merge", fill = "right")

final.called$Mutect2_gt_alt_fwd <- as.numeric(final.called$Mutect2_gt_alt_fwd)
final.called$Mutect2_gt_alt_rev <- as.numeric(final.called$Mutect2_gt_alt_rev)
final.called$Vardict_gt_ALD_forw <- as.numeric(final.called$Vardict_gt_ALD_forw)
final.called$Vardict_gt_ALD_rev <- as.numeric(final.called$Vardict_gt_ALD_rev)
final.called$Varscan2_gt_alt_fwd <- as.numeric(final.called$Varscan2_gt_alt_fwd)
final.called$Varscan2_gt_alt_rev <- as.numeric(final.called$Varscan2_gt_alt_rev)

## strand bias + must have alt depth of >=6
## 1 means passed, 0 means failed
final.called <- final.called %>%
  mutate(Varscan2_SB = 
           ifelse(Varscan2_gt_alt_fwd/(Varscan2_gt_alt_fwd+Varscan2_gt_alt_rev)<.1 |
                    Varscan2_gt_alt_fwd/(Varscan2_gt_alt_fwd+Varscan2_gt_alt_rev)>.9,0,
                  ifelse((Varscan2_gt_alt_fwd+Varscan2_gt_alt_rev)<6,0,1)),
         Mutect2_SB = 
           ifelse(Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
                    Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,0,
                  ifelse((Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<6,0,1)),
         Vardict_SB = 
           ifelse(Vardict_gt_ALD_forw/(Vardict_gt_ALD_forw+Vardict_gt_ALD_rev)<.1 |
                    Vardict_gt_ALD_forw/(Vardict_gt_ALD_forw+Vardict_gt_ALD_rev)>.9,0,
                  ifelse((Vardict_gt_ALD_forw+Vardict_gt_ALD_rev)<6,0,1)),
  )

final.called$passed <-
  apply(final.called[,c("Mutect2_PASS","Varscan2_PASS","Vardict_PASS")],
        1,
        #function(x) (sum(x)>=2 & x["Mutect2_PASS"]))
        function(x) (sum(x)>=2))
sum(final.called$passed)


final.passed <- final.called # passed now less than 6394 to 6032


###########################################################################################
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(R453Plus1Toolbox)
#HG38 <- readDNAStringSet("~/Bolton/PON_PU/all_sequences.fa")

final.passed$context_5 = getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(seqnames = final.passed$CHROM, 
                                                                     ranges = IRanges(start=final.passed$POS-10, end=final.passed$POS-1)))
final.passed$context_3 = getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(seqnames = final.passed$CHROM, 
                                                                     ranges = IRanges(start=final.passed$POS+nchar(final.passed$REF), end=final.passed$POS+nchar(final.passed$REF)+9)))
final.passed$context_5_3 = getSeq(BSgenome.Hsapiens.UCSC.hg38, GRanges(seqnames = final.passed$CHROM, 
                                                                       ranges = IRanges(start=final.passed$POS-5, end=final.passed$POS+nchar(final.passed$REF)+4)))
final.passed$dust_score_5 <- R453Plus1Toolbox::complexity.dust(final.passed$context_5)
final.passed$dust_score_3 <- R453Plus1Toolbox::complexity.dust(final.passed$context_3)
final.passed$dust_score_5_3 <- R453Plus1Toolbox::complexity.dust(final.passed$context_5_3)
final.passed$dust_score <- apply(final.passed[,c("dust_score_5", "dust_score_3","dust_score_5_3")],1,max)

sameLetters <- function (x) {
  letter <- substr(x, 1, 1)
  for (pos in 1:nchar(x)){
    if (substr(x, pos, pos) != letter) { return(FALSE) }
  }
  return(TRUE)
}


library(stringr)
final.passed <- final.passed %>%
  dplyr::mutate(
    # Determines if the Downstream sequence has 3 of the same nucleotides in a row
    case_NXXX = str_detect(context_3, regex('^[A]{3,}|^[G]{3,}|^[T]{3,}|^[C]{3,}')) &
      # Checks if the last character of the mutation completes the homopolymer
      substr(ALT, nchar(ALT), nchar(ALT)) == substr(context_3, 1, 1),
    # Sandwich Mutation
    case_XNXX = ifelse(nchar(ALT) == 1, TRUE, FALSE) & # Can only occur with Single Nucleotides
      # Check if the last character of the upstream + the mutation + the first two characters of downstream
      # completes the homopolymer
      paste0(substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
             ALT,
             substr(context_3, start = 1, stop = 2)) %>% 
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Second Sandwich Case
    case_XXNX = ifelse(nchar(ALT) == 1, TRUE, FALSE) & 
      paste0(substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
             ALT,
             substr(context_3, start = 1, stop = 1)) %>% 
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Determines if the Upstream sequence ends with 3 of the same nucleotides in a row
    case_XXXN = sapply(substr(context_5, start = nchar(context_5)-2, stop = nchar(context_5)), sameLetters) &
      # Checks if the first character of the mutation completes the homopolymer
      substr(ALT, 1, 1) == substr(context_5, start = nchar(context_5), stop = nchar(context_5)),
    # Same as Case NXXX but for dinucleotides
    case_NNXX = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & str_detect(context_3, regex('^[A]{2,}|^[G]{2,}|^[T]{2,}|^[C]{2,}')) &
      substr(ALT, nchar(ALT)-1, nchar(ALT)) == substr(context_3, 1, 2),
    # Sandiwch Case
    case_XNNX = ifelse(nchar(ALT) == 2, TRUE, FALSE) & sapply(ALT, sameLetters) & 
      paste0(
        substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)),
        ALT,
        substr(context_3, start = 1, stop = 2)) %>% 
      str_detect(regex('[A]{4,}|[G]{4,}|[T]{4,}|[C]{4,}')),
    # Same as Case XXXN but for dinucleotides
    case_XXNN = ifelse(nchar(ALT) >= 2, TRUE, FALSE) & sapply(substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5)), sameLetters) &
      substr(ALT, 1, 2) == substr(context_5, start = nchar(context_5)-1, stop = nchar(context_5))
  )

final.passed <- final.passed %>%
  mutate(
    complexity_filters =
      !(case_NXXX) &
      !(case_XNXX) &
      !(case_XXNX) &
      !(case_XXXN) &
      !(case_NNXX) &
      !(case_XNNX) &
      !(case_XXNN) &
      dust_score < 7
  )

final.passed$context_5 <- as.character(final.passed$context_5)
final.passed$context_3 <- as.character(final.passed$context_3)
final.passed$context_5_3 <- as.character(final.passed$context_5_3)

## min VAFs > 35% remove
mutect_VAF <- grep("Mutect2_gt_AF", colnames(final.passed))
vardict_VAF <- grep("Vardict_gt_AF", colnames(final.passed))
varscan_VAF <- grep("Varscan2_gt_FREQ", colnames(final.passed))
VAFS <- c(mutect_VAF, vardict_VAF, varscan_VAF)
final.passed$min.under.0.35 <- apply(final.passed[,VAFS],1,function(x) min(x, na.rm = T) < 0.35)
final.passed$max.over.0.02 <- apply(final.passed[,VAFS],1,function(x) max(x, na.rm = T) >= 0.02)

#################################################################################
## test vardict from PON 2 @ >2%
library(VariantAnnotation)
testvcf <- VariantAnnotation::readVcf(VcfFile("~/Bolton/final_vcfs/all_coding_sorted.vcf.gz"))
# gr.called <- final %>% filter(called == 1)
######## gr <- GRanges(seqnames = gr.called$CHROM, 
#               ranges = IRanges(start=gr.called$POS, 
#                                end = gr.called$POS),
#               REF = gr.called$REF,
#               ALT = gr.called$ALT)
# gr <- GRanges(seqnames = final.passed$CHROM,
#               ranges = IRanges(start=final.passed$POS,
#                                end = final.passed$POS),
#               REF = final.passed$REF,
#               ALT = final.passed$ALT)
# gr <- sortSeqlevels(gr)
# gr <- sort(gr)
# final.vcf <- VariantAnnotation::VCF(rowRanges = gr, fixed = DataFrame(REF=DNAStringSet(gr$REF),
#                                                                       ALT=DNAStringSetList(as.list(gr$ALT)),
#                                                                       QUAL = 1))
# VariantAnnotation::writeVcf(final.vcf, "~/Bolton/final_vcfs/under35percent.vcf") ## really all coding rows 378,241


### need this file
new.PON <- read.table("~/Bolton/final_vcfs/filter.PON.vardict.txt", header = F, sep = "\t")
colnames(new.PON) <- colnames(final.passed)[1:4]
new.PON$PoN.2at2percent <- FALSE

# library(dplyr)
library(tidyverse)

final.passed <- left_join(final.passed, new.PON,
                          by = c("CHROM"="CHROM", "POS"="POS", 
                                 "REF"="REF", "ALT"="ALT"))
final.passed$PoN.2at2percent[is.na(final.passed$PoN.2at2percent)] <- T
#################################################################################
## Strand Bias 2 or more pass SB
final.passed$alt_strand_counts_min_2_callers <- apply(final.passed[,c("Mutect2_SB","Varscan2_SB","Vardict_SB")],
                                        1,
                                        function(x) sum(x, na.rm = T)>=2)
final.passed$alt_strand_counts_min_1_caller_only <- apply(final.passed[,c("Mutect2_SB","Varscan2_SB","Vardict_SB")],
                                                          1,
                                                          function(x) sum(x, na.rm = T)>=1)

mutect_VAF <- grep("Mutect2_gt_AF", colnames(final.passed))
vardict_VAF <- grep("Vardict_gt_AF", colnames(final.passed))
varscan_VAF <- grep("Varscan2_gt_FREQ", colnames(final.passed))
VAFS <- c(mutect_VAF, vardict_VAF, varscan_VAF)

final.passed$MIN_VAF <- apply(final.passed[,VAFS],1,function(x) min(x, na.rm = T))
final.passed$MAX_VAF <- apply(final.passed[,VAFS],1,function(x) max(x, na.rm = T))

final.passed$Rescuable <- apply(final.passed,1,function(x){
  if (as.numeric(x["MAX_VAF"])<.02 | as.numeric(x["MIN_VAF"])>.35) {
    return(F)
  } else {
    if (x["Vardict_CALLER"]==1 & (nchar(x["REF"])>10 | nchar(x["ALT"])>10) & x["Existing_variation_VEP"]=="") {
      return(F)
    } else {
      return(T)
    }
  }
})

nrow(final.passed %>%
       filter(
         gnomad.age.keep &
           #called &
           passed &
           min.under.0.35 & #1
           max.over.0.02 & #2
           alt_strand_counts_min_2_callers & #3
           PoN.2at2percent & #4
           #Mutect2_CALLER &
           #Mutect2_PASS &
           complexity_filters #5
       ))


write.table(final.passed, "final_all_rows_all_filters.tsv", row.names = F, quote = F, sep= "\t")
#################################################################################
final.passed <- read.table("~/Bolton/HIV/final_all_rows_all_filters.tsv", sep = "\t", header = T, quote = "")
# ##########
# 
# sum(final.passed$Mutect2_CALLER)
# sum(final.passed$Mutect2_CALLER & final.passed$Varscan2_CALLER)
# sum(final.passed$Mutect2_CALLER & final.passed$Vardict_CALLER)
# sum(final.passed$Mutect2_CALLER & final.passed$Varscan2_CALLER & final.passed$Vardict_CALLER)
# 
# sum(final.passed$Mutect2_PASS)
# sum(final.passed$Mutect2_PASS & final.passed$Varscan2_PASS)
# sum(final.passed$Mutect2_PASS & final.passed$Vardict_PASS)
# sum(final.passed$Mutect2_PASS & final.passed$Varscan2_PASS & final.passed$Vardict_PASS)
# 
# ###########
final.passed.rescuable.didpass <- final.passed %>%
  filter(
    gnomad.age.keep &
      #passed &
      #dust_score < 7 &
      alt_strand_counts_min_2_callers &
      min.under.0.35 & 
      PoN.2at2percent &
      Rescuable
  )

sum(final.passed.rescuable.didpass$Mutect2_PASS)
sum()
##########
# send to annotate_PD.R script
source("~/Bolton/R_stuff/annotate.PD.R")
final.passed.rescuable.notpassed <- final.passed %>%
  filter(
      #gnomAD_MAX.Stringent.0007 &
      gnomad.age.keep &
      !passed &
      called &
      alt_strand_counts_min_1_caller_only &
      #dust_score < 7 &
      min.under.0.35 & 
      PoN.2at2percent &
      Rescuable & 
      Mutect2_CALLER
  )
final.passed.rescuable.notpassed$source_rescue <- "not_passed_rescued"
M.notpassed <- annotate.PD(final.passed.rescuable.notpassed)
final.passed.rescuable.notpassed <- left_join(final.passed.rescuable.notpassed,
                                              M.notpassed %>%
                                                # dplyr::select("CHROM","POS","REF","ALT","SAMPLE","ch_my_pd","ch_pd2",
                                                #               "VariantClass","AAchange"),
                                                dplyr::select("CHROM","POS","REF","ALT","SAMPLE","CosmicCount","heme_cosmic_count","MDS",
                                                            "AML","MPN","ch_my_pd","ch_pd","ch_pd2","VariantClass","AAchange","Gene"),
                                              by=c("CHROM","POS","REF","ALT","SAMPLE"))


final.passed.rescuable.notpassed.onlyMutect2 <- final.passed %>%
  filter(
    (gnomad.age.keep &
    !called &
    !passed &
    alt_strand_counts_min_1_caller_only &
    #dust_score < 7 &
    min.under.0.35 & 
    PoN.2at2percent &
    Rescuable & 
    final.passed$Mutect2_CALLER) | (CHROM=="chr20"&POS==32434638))
final.passed.rescuable.notpassed.onlyMutect2$source_rescue <- "mutect2Only_not_passed_rescued"
M.notpassed.onlyMutect2 <- annotate.PD(final.passed.rescuable.notpassed.onlyMutect2)
final.passed.rescuable.notpassed.onlyMutect2 <- left_join(final.passed.rescuable.notpassed.onlyMutect2,
                                                          M.notpassed.onlyMutect2 %>%
                                                            dplyr::select("CHROM","POS","REF","ALT","SAMPLE","CosmicCount","heme_cosmic_count","MDS",
                                                                          "AML","MPN","ch_my_pd","ch_pd","ch_pd2","VariantClass","AAchange","Gene"),
                                            by=c("CHROM","POS","REF","ALT","SAMPLE"))


## passed and "rescuable" for annotating, passed everything but complexity
final.passed.rescuable.didpass <- final.passed %>%
  filter(
    gnomad.age.keep &
    passed &
    #dust_score < 7 &
    alt_strand_counts_min_2_callers &
    min.under.0.35 & 
    PoN.2at2percent &
    Rescuable
  )
sum(final.passed.rescuable.didpass$Mutect2_PASS)

final.passed.rescuable.didpass$source_rescue <- "passed_rescued"
M.passed <- annotate.PD(final.passed.rescuable.didpass)
final.passed.rescuable.didpass <- left_join(final.passed.rescuable.didpass,
                                            M.passed %>%
                                              dplyr::select("CHROM","POS","REF","ALT","SAMPLE","CosmicCount","heme_cosmic_count","MDS",
                                                            "AML","MPN","ch_my_pd","ch_pd","ch_pd2","VariantClass","AAchange","Gene"),
                                            by=c("CHROM","POS","REF","ALT","SAMPLE"))

###################################################################################################

pd <- rbind(final.passed.rescuable.notpassed,
            final.passed.rescuable.notpassed.onlyMutect2,
            final.passed.rescuable.didpass) %>%
  dplyr::filter(ch_pd2.y>=1 & MAX_VAF>=.02 & 
           #!(SYMBOL_VEP %like% "HLA") &
           #VariantClass == "Missense" &
           !(POS %in% c(152265208, 152273811, 44873619,  ##  red by Kelly
                        3850457, 43338229, # 3850457 failed low base quality in M2, 43338229 snpish in gnomad
                          111598951, 111598965, 111598967, 111598973, 111598976, ## mutect only
                          132062544, 132062545, 132062547, 132062548, 132062553, ## mutect only
                          86510608, 86510611, 33301851, 137871488, 156778288, 41749340, ## mutect only 41749340,
                          137871488, 156778288, 156778291, 2039792, 2039778, ## mutect only
                        31226550, 42472359, 5126715 # missense review
                        ))
         ) 
pd <- pd %>%
  mutate(change=paste(CHROM,POS,REF,ALT,SAMPLE, sep = ":"))
# write.table(pd, "~/Bolton/HIV/pd.tsv", sep = "\t", row.names = F, quote = F)
pd.orig <- read.table("~/Bolton/HIV/pd.tsv", sep = "\t", quote = "", header = T)
##########

# sum(pd$Mutect2_CALLER)
# sum(pd$Mutect2_CALLER & pd$Varscan2_CALLER)
# sum(pd$Mutect2_CALLER & pd$Vardict_CALLER)
# sum(pd$Mutect2_CALLER & pd$Varscan2_CALLER & pd$Vardict_CALLER)
# 
# sum(pd$Mutect2_PASS)
# sum(pd$Mutect2_PASS & pd$Varscan2_PASS)
# sum(pd$Mutect2_PASS & pd$Vardict_PASS)
# sum(pd$Mutect2_PASS & pd$Varscan2_PASS & pd$Vardict_PASS)
# 
# sum(pd$Mutect2_PASS) / sum(pd$Mutect2_CALLER)
# sum(pd$Mutect2_PASS & pd$Varscan2_PASS) / sum(pd$Mutect2_CALLER & pd$Varscan2_CALLER)
# sum(pd$Mutect2_PASS & pd$Vardict_PASS) / sum(pd$Mutect2_CALLER & pd$Vardict_CALLER)
# sum(pd$Mutect2_PASS & pd$Varscan2_PASS & pd$Vardict_PASS) / sum(pd$Mutect2_CALLER & pd$Varscan2_CALLER & pd$Vardict_CALLER)

###########
pd2 <- pd2 %>%
  mutate(change=paste(CHROM,POS,REF,ALT,SAMPLE, sep = ":"))

new.pd <- pd[!(pd$change %in% pd2$change),]
## PD
ages.HIV.case.control <- readxl::read_xlsx("~/Bolton/HIV/id2.xlsx")
ages.HIV.case.control[!duplicated(ages.HIV.case.control$id),]


ages.HIV <- read.csv("~/Bolton/HIV/id2.csv")
ages.HIV <- ages.HIV[!duplicated(ages.HIV$id),]
graph.pd <- pd %>% dplyr::select(SAMPLE) %>% dplyr::group_by(SAMPLE) %>% 
  dplyr::mutate(n = ifelse(n()>=1,1,0)) %>% dplyr::distinct()# %>% ungroup()
ages.HIV <- left_join(ages.HIV[,c("id","Age")], graph.pd, by = c("id"="SAMPLE"))
# ages.HIV.age.order <- ages.HIV[order(ages.HIV$Age),]
# test5 <- ages.HIV.age.order %>%
#   mutate(label = ifelse(is.na(n),"",Age))
# x <- barplot(ages.HIV.age.order$n)
# labs <- test5$label
# text(cex=.85, x=x-.25, y=-0.25, labs, xpd=TRUE, srt=45)
#graph.pd <- left_join(graph, ages.HIV[,c("id","Age")], by = c("SAMPLE"="id"))
#write.table(graph.pd, "graph.pd.tsv", row.names = F, sep = "\t")
#write.table(ages.HIV, "ages.HIV.tsv", row.names = F, sep = "\t")
#ages.HIV <- left_join(ages.HIV[,c("id","Age")], graph.pd, by = c("id"="SAMPLE"))
#graph.pd <- left_join(graph, ages.HIV[,c("id","Age")], by = c("SAMPLE"="id"))
#ages.HIV <- read.table("~/Bolton/data/ages.HIV.tsv", header=T, sep="\t")
ages.HIV$bins <- cut(x=ages.HIV$Age, breaks=c(40,45,55,65,Inf), include.lowest = T)
Bygroup.count = data.frame(n=tapply(ages.HIV$n, ages.HIV$bins, function(x) sum(x, na.rm = T)))
Bygroup.total = data.frame(ages.HIV %>% group_by(bins) %>% dplyr::summarise(total=as.numeric(n())))
row.names(Bygroup.total) <- Bygroup.total$bins
df <- cbind(Bygroup.count, Bygroup.total[,2])
colnames(df)[2] <- "total"
df$prop <- df$n/df$total
df <- df %>% rowwise() %>%
  dplyr::mutate(lower = prop.test(n, total, conf.level = .95)$conf.int[1], upper = prop.test(n, total, conf.level = .95)$conf.int[2])
df <- as.data.frame(df)
rownames(df) <- Bygroup.total$bins
df$bin <- rownames(df)
df$bin <- factor(df$bin, levels = df$bin)
df$lab <- paste0(df$n,"(",round(df$prop*100,2),"%",")")

ggplot(df, aes(x = bin, y = prop, group = 1)) +
  geom_bar(stat='identity', fill = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.4, fill="green") +
  geom_text(aes(x=bin,y=prop,label=lab),vjust=-1) +
  ggtitle("All Passed") +
  #ylim(c(0, max(df$prop)*1.1)) +
  xlab("Age") +
  ylab("Variants per sample with Any Mutation") +
  theme(plot.title = element_text(hjust = 0.5))

write.table(pd %>% filter(POS %in% c(41749340, 25239199, 2039778),
                          source_rescue=="mutect2Only_not_passed_rescued"),
            "weak_evidence.tsv", row.names=F, sep="\t")
## pd by gene
ggplot(pd %>% dplyr::select(SYMBOL_VEP) %>%
         group_by(SYMBOL_VEP) %>%
         dplyr::summarise(n=n()), aes(reorder(x = SYMBOL_VEP,-n), y = n, fill = blues9)) +
  geom_bar(stat = 'identity', position = "dodge", color = 'black', fill = pal_nejm('default')(1), size = 0.1) + 
  ylab("Number of Mutations") +
  xlab(' ')


## up stream and downstream of mutations for ddPCR (if needed)
pd <- pd %>%
  mutate(CHROM2=CHROM, flank5=POS-100,flank3=nchar(REF)+POS+100-1)
forbam <- pd %>%
  mutate(CHROM2=CHROM, 
         flank5=paste0(POS-100,"-",POS-1),
         ALT_SEQ_ddPCR=ALT,
         flank3=ifelse(
           nchar(ALT)>nchar(REF), 
           paste0(POS+1,"-",POS+100),
           paste0(POS+(nchar(REF)-nchar(ALT)+1),"-",POS+(nchar(REF)-nchar(ALT)+100))
           )
         #,
         #AAchange= ifelse(is.na(AAchange.x),AAchange.y,AAchange.x)
  ) %>%
  dplyr::select(CHROM,POS,REF,ALT,MAX_VAF,SYMBOL_VEP,SAMPLE,VARIANT_CLASS_VEP,
                AAchange, 
                #AAchange.y,
                HGVSc_VEP,
                HGVSp_VEP,
                CHROM2,flank5,ALT_SEQ_ddPCR,flank3)
chrOrder<-c(paste("chr",1:22,sep=""),"chrX","chrY","chrM")
forbam$CHROM <- factor(forbam$CHROM, levels=chrOrder)
forbam <- forbam[order(forbam$CHROM,forbam$POS),]
write.table(forbam[,-12], "~/Bolton/HIV/PD_for_Grant.tsv", sep="\t", row.names = FALSE, quote=F)
nwcs489 <- readxl::read_xlsx("~/Bolton/HIV/nwcs489.xlsx")
colnames(nwcs489)
as.numeric(gsub("-","",stringr::str_extract(ages.HIV$id,"-.*-")))
all(nwcs489$`Random\r\nparticipant\r\nnumber` %in% as.numeric(gsub("-","",stringr::str_extract(ages.HIV$id,"-.*-"))))
ages.HIV.case.control <- readxl::read_xlsx("~/Bolton/HIV/id2.xlsx")
test <- ages.HIV.case.control %>%
  group_by(gender, `Race/Ethnicity`) %>%
  arrange(gender,`Race/Ethnicity`, Age) %>%
  dplyr::select(id,name,gender,`Race/Ethnicity`, Age, `Case/Control`)

pd$avgVAF <- apply(pd[,getVAFs(pd)],1,function(x) mean(x, na.rm=T))

pd.removed <- rbind(final.passed.rescuable.notpassed,
                       final.passed.rescuable.notpassed.onlyMutect2,
                       final.passed.rescuable.didpass) %>%
  filter(ch_pd2.test>=1 & MAX_VAF>=.02 & 
           !(SYMBOL_VEP %in% c("HLA-A", "HLA-B", "HLA-C")) & 
  #VariantClass == "Missense" &
  (POS %in% c(152265208, 152273811, 44873619,  ##  red by Kelly
               3850457, 43338229, # 3850457 failed low base quality in M2, 43338229 snpish in gnomad
               111598951, 111598965, 111598967, 111598973, 111598976, ## mutect only
               132062544, 132062545, 132062547, 132062548, 132062553, ## mutect only
               86510608, 86510611, 33301851, 137871488, 156778288, 41749340, ## mutect only  
               137871488, 156778288, 156778291, 2039792, 2039778, ## mutect only
               31226550, 42472359, 5126715 # missense review
  ))) 
table(pd.removed$source_rescue)
###################################################################################################
# some passed are not "rescuable"?? these passed everything but some are not rescuable, i.e. not > 2% or complex vardict
#final.passed$blacklist <- final.coding.gnomad.sorted.gnomadAge$blacklist
# final.passed$passed2 <-
#   apply(final.passed[,c("Mutect2_PASS","Varscan2_PASS","Vardict_PASS")],
#         1,
#         function(x) (sum(x)>=2))
# final.passed.everything <- final.passed %>%
#   filter(
#     (#gnomAD_MAX.Stringent.0007 &
#     gnomad.age.keep &
#     #passed &
#     passed2 &
#     Mutect2_PASS &
#     alt_strand_counts_min_2_callers &
#     min.under.0.35 &
#     max.over.0.02 &
#     complexity_filters & # this is not included above in final.passed.rescuable.didpass rescue
#     PoN.2at2percent 
#     & !(blacklist))
#   )
# write.table(final.passed.everything, "~/Bolton/data/final.passed.everything.tsv", sep="\t", row.names = F)
# final.passed.everything <- read.table("~/Bolton/data/final.passed.everything.tsv", sep="\t", header=T)
# library(dplyr)
# final.passed.everything$source_rescue <- "passed_all"
# M.allpass <- annotate.PD(final.passed.everything)
# final.passed.everything <- left_join(final.passed.everything,
#                                      M.allpass %>%
#                                               dplyr::select("CHROM","POS","REF","ALT","SAMPLE","ch_my_pd","ch_pd2","VariantClass"),
#                                             by=c("CHROM","POS","REF","ALT","SAMPLE"))
#############################
## using refflat to look at down and upstream regions of exons to get average sequencing depth over the exons
x <- barplot(sort(table(final.passed.everything$SYMBOL_VEP), decreasing = T)[1:50], xaxt="n")
labs <- names(sort(table(final.passed.everything$SYMBOL_VEP), decreasing = T)[1:50])
text(cex=.75, x=x-.25, y=-2.25, labs, xpd=TRUE, srt=45)
# genes <- rtracklayer::import("~/Bolton/Mocha_files/refFlat.bed")
genes <- rtracklayer::import.bed("~/Bolton/Mocha_files/refFlat.exons.bed")
final.passed.everything.gr <- GRanges(seqnames = final.passed.everything$CHROM,
                                      ranges = IRanges(start=final.passed.everything$POS, 
                                                       end=final.passed.everything$POS), 
                                      mcols = final.passed.everything[,3:ncol(final.passed.everything)])

##
genes@elementMetadata <-  DataFrame(genes@elementMetadata$name, start(genes@ranges), end(genes@ranges))
genes.ovr <- plyranges::join_overlap_left(final.passed.everything.gr, genes)
genes.df <-  as.data.frame(genes.ovr)[,c(1,2,6:ncol(as.data.frame(genes.ovr)))]
columns <- colnames(final.passed.everything)
colnames(genes.df) <- c(columns, c("gene.name","start","end"))
genes.df <- genes.df[!duplicated(genes.df[,c("CHROM","POS","REF","ALT","SAMPLE")]),]
genes.df[,c("SYMBOL_VEP")] == genes.df[,c("gene.name")] 
genes.df[!(genes.df[,c("SYMBOL_VEP")] == gsub("_.*","",genes.df[,c("gene.name")])),c("CHROM","POS","SYMBOL_VEP","gene.name")]
hist((genes.df$end-genes.df$start)[(genes.df$end-genes.df$start)<2000])
genes.df <- genes.df %>%
  mutate(depth.start = pmax(POS-135,start, na.rm=T),
         depth.end = pmin(POS+135,end,na.rm=T))

#library(Rsamtools)
write.table(genes.df[,c("CHROM","depth.start", "depth.end", "POS","REF","ALT","SAMPLE","SYMBOL_VEP","gene.name")],
            "~/Bolton/HIV/depth.bed", row.names = F, quote = F, col.names = F, sep = "\t")

ave.depth <- sapply(with(genes.df,paste0("samtools depth /Volumes/bolton/Active/projects/HIV/input/",
                     SAMPLE,".bam", " -r ", CHROM, ":",
                     depth.start,"-",depth.end, " | awk '{sum+=$3;} END {print sum/NR}'")),
       function(x) {
         as.numeric(system(command = x, intern = T))
       }
)
genes.df$ave.depth <- ave.depth
x <- barplot(sort(table(genes.df$SYMBOL_VEP[genes.df$ave.depth <= 1000]), decreasing = T)[1:30], xaxt="n")
labs <- names(sort(table(genes.df$SYMBOL_VEP[genes.df$ave.depth <= 1000]), decreasing = T)[1:30])
text(cex=.75, x=x-.25, y=-2.25, labs, xpd=TRUE, srt=45)
  
# system(command = "samtools depth /Volumes/bolton/Active/projects/HIV/input/CHAJ-814863-171.bam -r chr1:952362-952601 | awk '{sum+=$3;} END {print sum/NR}'",
#      intern = T)

low.depth <- genes.df[genes.df$ave.depth <= 1000,]

#######################################################################################
## this is for exon coordinates from UCSC table browser: refFlat
library(dplyr )
final.passed.everything <- read.table("~/Bolton/data/final.passed.everything.tsv", sep="\t", header=T)
genes <- rtracklayer::import.bed("~/Bolton/Mocha_files/refFlat.exons.bed")
UKBB <- rtracklayer::import("/Users/brian/Bolton/CWL_TESTS/output.fa")
library(seqinr)
library("Biostrings")

s = readDNAStringSet("nm.fasta")
UKBB <- readDNAStringSet("/Users/brian/Bolton/CWL_TESTS/output.fa")
C = complexity.dust(UKBB)
final.passed.everything.gr <- GRanges(seqnames = final.passed.everything$CHROM,
                                      ranges = IRanges(start=final.passed.everything$POS, 
                                                       end=final.passed.everything$POS), 
                                      mcols = final.passed.everything[,3:ncol(final.passed.everything)])
genes@elementMetadata <- DataFrame(genes@elementMetadata$name, start(genes@ranges), end(genes@ranges))
genes.ovr <- plyranges::join_overlap_left(final.passed.everything.gr, genes)
columns <- colnames(final.passed.everything)
final.passed.everything1 <-  as.data.frame(genes.ovr)[,c(1,2,6:ncol(as.data.frame(genes.ovr)))]
colnames(final.passed.everything1) <- c(columns, c("gene.name","start","end"))
final.passed.everything1 <- final.passed.everything1[!duplicated(final.passed.everything1[,c("CHROM","POS","REF","ALT","SAMPLE")]),]
sum(final.passed.everything1[,c("SYMBOL_VEP")] == gsub("_.*","",final.passed.everything1[,c("gene.name")]), na.rm = T)
final.passed.everything1[!(final.passed.everything1[,c("SYMBOL_VEP")] == gsub("_.*","",final.passed.everything1[,c("gene.name")])),c("CHROM","POS","SYMBOL_VEP","gene.name")]
#hist((final.passed.everything1$end-final.passed.everything1$start)[(final.passed.everything1$end-final.passed.everything1$start)<2000])
final.passed.everything1 <- final.passed.everything1 %>%
  mutate(depth.start = pmax(POS-135,start, na.rm=T),
         depth.end = pmin(POS+135,end,na.rm=T))



dups <- rtracklayer::import("/Users/brian/Bolton/final_vcfs/dup.grch38.bed.gz")
dups$dups <- T
repetitive <- rtracklayer::import("/Users/brian/Bolton/HIV/simpleRepeat.bed")
repetitive$repetitive <- T
#repeatMasker <- rtracklayer::import("/Users/brian/Bolton/HIV/repeatMasker.bed")
repeatMasker <- rtracklayer::import("/Users/brian/Bolton/HIV/repeatMaskerJoinedCurrent.bed")
repeatMaskerRegion <- read.table("/Users/brian/Bolton/HIV/repeatMaskerJoinedCurrent2.bed", sep="\t", quote="", header = F)
repeatMasker$region <- repeatMaskerRegion$V1
repeatMasker$repeatMasker <- T

final.passed.everything1.gr <- GRanges(seqnames = final.passed.everything1$CHROM,
                  ranges = IRanges(start=final.passed.everything1$POS, 
                                   end=final.passed.everything1$POS), 
                  mcols = final.passed.everything1[,3:ncol(final.passed.everything1)])
test.dups <- plyranges::join_overlap_left(final.passed.everything1.gr, dups)
test.dups.rep <- plyranges::join_overlap_left(test.dups, repetitive)
test.dups.rep.repMasker <- plyranges::join_overlap_left(test.dups.rep, repeatMasker)

columns <- colnames(final.passed.everything1)
final.passed.everything2 <- as.data.frame(test.dups.rep.repMasker)[,c(1,2,6:ncol(as.data.frame(test.dups.rep.repMasker)))]
colnames(final.passed.everything2) <- c(columns, c("dup_score","dups","trf","trf_score","simpleTandemRepeat","RMName",
                                                   "RMScore","RMRegion","repeatMasker"))
final.passed.everything2$dups[is.na(final.passed.everything2$dups)] <- F
final.passed.everything2$simpleTandemRepeat[is.na(final.passed.everything2$simpleTandemRepeat)] <- F
final.passed.everything2$repeatMasker[is.na(final.passed.everything2$repeatMasker)] <- F
# test.dup <- final.passed.everything2[duplicated(final.passed.everything2[,c("CHROM","POS","REF","ALT","SAMPLE")])==TRUE|duplicated(final.passed.everything2[,c("CHROM","POS","REF","ALT","SAMPLE")], fromLast = TRUE),]
# sum(duplicated(final.passed.everything2[,c("CHROM","POS","REF","ALT","SAMPLE")]))
final.passed.everything2 <- final.passed.everything2[!duplicated(final.passed.everything2[,c("CHROM","POS","REF","ALT","SAMPLE")]),]
write.table(final.passed.everything2, "final.passed.everything2.tsv", row.names = F, quote = F, sep="\t")


ave.depth <- sapply(with(final.passed.everything2,paste0("samtools depth /Volumes/bolton/Active/projects/HIV/input/",
                                         SAMPLE,".bam", " -r ", CHROM, ":",
                                         depth.start,"-",depth.end, " | awk '{sum+=$3;} END {print sum/NR}'")),
                    function(x) {
                      as.numeric(system(command = x, intern = T))
                    }
)
final.passed.everything2$ave.depth

final.passed.everything3 <- final.passed.everything2 %>% dplyr::filter(!dups & !simpleTandemRepeat)
#final.passed.everything3 <- final.passed.everything2 %>% dplyr::filter(!dups & !simpleTandemRepeat & !repeatMasker)
## All Depths
x <- barplot(sort(table(final.passed.everything3$SYMBOL_VEP), decreasing = T)[1:40], xaxt="n")
labs <- names(sort(table(final.passed.everything3$SYMBOL_VEP), decreasing = T)[1:40])
text(cex=.75, x=x-.25, y=-2.25, labs, xpd=TRUE, srt=45)



x <- barplot(sort(table(final.passed.everything3$SAMPLE), decreasing = T)[1:30], xaxt="n")
labs <- names(sort(table(final.passed.everything3$SAMPLE), decreasing = T)[1:30])
text(cex=.75, x=x-.25, y=-2.25, labs, xpd=TRUE, srt=45)
sort(table(final.passed.everything3$SYMBOL_VEP[final.passed.everything3$SAMPLE==labs[1]]), decreasing = T)[1:5]



# last.test <- left_join(final.passed.everything3,
#                       low.depth %>%
#                         dplyr::select("CHROM","POS","REF","ALT","SAMPLE","ave.depth"),
#                       by=c("CHROM","POS","REF","ALT","SAMPLE"))
# last.test <- last.test[!is.na(last.test$ave.depth),]
# x <- barplot(sort(table(last.test$SYMBOL_VEP), decreasing = T)[1:30], xaxt="n")
# labs <- names(sort(table(last.test$SYMBOL_VEP), decreasing = T)[1:30])
# text(cex=.75, x=x-.25, y=-2.25, labs, xpd=TRUE, srt=45)



library(dplyr)

# symbols
x <- barplot(sort(table(final.passed.everything$SYMBOL_VEP), decreasing = T)[1:30], xaxt="n")
labs <- names(sort(table(final.passed.everything$SYMBOL_VEP), decreasing = T)[1:30])
text(cex=.75, x=x-.25, y=-8.25, labs, xpd=TRUE, srt=45)

mutect_VAF <- grep("Mutect2_gt_AF", colnames(final.passed.everything))
vardict_VAF <- grep("Vardict_gt_AF", colnames(final.passed.everything))
varscan_VAF <- grep("Varscan2_gt_FREQ", colnames(final.passed.everything))
VAFS <- c(mutect_VAF, vardict_VAF, varscan_VAF)
final.passed.everything$avg.VAF <- apply(final.passed.everything[,VAFS], 1, function(x) mean(x, na.rm=T))

final.passed.everything <- final.passed.everything %>% mutate(
  VAF_bin = case_when(
    avg.VAF < 0.05 ~ '<0.05',
    avg.VAF >= 0.05 & avg.VAF < 0.1 ~ '0.05-0.10',
    avg.VAF >= 0.1 & avg.VAF < 0.25 ~ '0.10-0.25',
    avg.VAF >= 0.25 & avg.VAF ~ '>=.25'
  ))


top30 <- names(sort(table(final.passed.everything3$SYMBOL_VEP), decreasing = T)[1:30])
final.passed.everything.df3 <- final.passed.everything3 %>%
  dplyr::select(SYMBOL_VEP, VAF_bin, gnomAD_MAX.Stringent.0007) %>% 
  dplyr::filter(SYMBOL_VEP %in% top30) %>%
  dplyr::count(SYMBOL_VEP, VAF_bin, gnomAD_MAX.Stringent.0007)
  #dplyr::count(SYMBOL_VEP, sort = TRUE)
final.passed.everything.df3$VAF_bin <- factor(final.passed.everything.df3$VAF_bin,
                                              levels = c("<0.05","0.05-0.10","0.10-0.25",">=.25"))
final.passed.everything.df3$SYMBOL_VEP <- factor(final.passed.everything.df3$SYMBOL_VEP,
                                                 levels=top30)

ggplot(final.passed.everything.df3, aes(x = SYMBOL_VEP, y = n, fill = gnomAD_MAX.Stringent.0007)) +
  geom_bar(stat = 'identity', position = "stack", color = 'black', size = 0.1) + 
  ylab("Number of Mutations") +
  facet_wrap(~VAF_bin, nrow=4) +
  xlab(' ') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


library(ggsci)
ggplot(pd %>% dplyr::select(SYMBOL_VEP) %>%
         group_by(SYMBOL_VEP) %>%
         dplyr::summarise(n=n()), aes(reorder(x = SYMBOL_VEP,-n), y = n, fill = blues9)) +
  geom_bar(stat = 'identity', position = "dodge", color = 'black', fill = pal_nejm('default')(1), size = 0.1) + 
  ylab("Number of Mutations") +
  xlab(' ')
# symbols passed vardict
x <- barplot(sort(table(final.passed.everything$SYMBOL_VEP[final.passed.everything$Vardict_PASS==1]), decreasing = T)[1:30], xaxt="n")
labs <- names(sort(table(final.passed.everything$SYMBOL_VEP), decreasing = T)[1:30])
text(cex=.75, x=x-.25, y=-8.25, labs, xpd=TRUE, srt=45)
##############################################################################################################
## vardict PON



##############################################################################################################
# SAMPLE all passsdd stuff
x <- barplot(sort(table(final.passed.everything3$SAMPLE), decreasing = T), xaxt="n")

labs <- names(sort(table(final.passed.everything3$SAMPLE), decreasing = T))
text(cex=.55, x=x-.25, y=-15.25, labs, xpd=TRUE, srt=90)
sort(table(final.passed.everything3$SAMPLE), decreasing = T)[1:12]
high.count <- names(sort(table(final.passed.everything3$SAMPLE), decreasing = T)[1:12])
ages.HIV <- readxl::read_xlsx("~/Bolton/HIV/metadata.xlsx")
test.HIV <- ages.HIV %>% 
  group_by(id) %>%
  dplyr::count(id, sort = TRUE)
  #%>%filter(n > 1)
test.HIV$n <- as.factor(test.HIV$n)
ggplot(test.HIV, aes(reorder(x = id,-Freq), y=Freq, fill=n)) +
  geom_bar(stat = "identity") +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        panel.grid.major = element_blank())  +
  guides(fill=guide_legend(title="Times Sequenced"))
  
test.HIV <- left_join(test.HIV, data.frame(sort(table(final.passed.everything3$SAMPLE), decreasing = T)), by = c("id"="Var1"))
test.HIV <- test.HIV[order(test.HIV$Freq, decreasing=T),]
high.count %in% test.HIV$id
  summarise(n=n()) %>%
  arrange(n, descen)
ages.HIV <- ages.HIV[!duplicated(ages.HIV$id),]



pd.in.passed <- final.passed.everything %>% 
  dplyr::filter(SYMBOL_VEP %in% unique(pd$SYMBOL_VEP))
table(pd$source_rescue)
test3 <- final.passed %>%
  filter(
    (#gnomAD_MAX.Stringent.0007 &
      gnomad.age.keep &
        #passed &
        passed2 &
        Mutect2_PASS &
        alt_strand_counts_min_2_callers &
        min.under.0.35 &
        max.over.0.02 &
        #complexity_filters & # this is not included above in final.passed.rescuable.didpass rescue
        PoN.2at2percent 
      & !(blacklist)) &
      SYMBOL_VEP %in% unique(pd$SYMBOL_VEP) &
      !(POS %in% c(152265208, 152273811, 44873619,  ##  red by Kelly
                   3850457, 43338229, # 3850457 failed low base quality in M2, 43338229 snpish in gnomad
                   111598951, 111598965, 111598967, 111598973, 111598976, ## mutect only
                   132062544, 132062545, 132062547, 132062548, 132062553, ## mutect only
                   86510608, 86510611, 33301851, 137871488, 156778288, 41749340, ## mutect only 41749340,
                   137871488, 156778288, 156778291, 2039792, 2039778, ## mutect only
                   31226550, 42472359, 5126715 # missense review
      ))
  )

## ALL PASSED
graph <- final.passed.everything3
graph <- graph %>% dplyr::select(SAMPLE) %>% dplyr::group_by(SAMPLE) %>% mutate(n = n()) %>% distinct()# %>% ungroup()
ages.HIV <- left_join(ages.HIV[,c("id","Age")], graph, by = c("id"="SAMPLE"))
graph <- left_join(graph, ages.HIV[,c("id","Age")], by = c("SAMPLE"="id"))
ages.HIV$bins <- cut(x=ages.HIV$Age, breaks=c(40,45,55,65,Inf), include.lowest = T)
Bygroup.count = data.frame(n=tapply(ages.HIV$n, ages.HIV$bins, function(x) sum(x, na.rm = T)))
Bygroup.total = data.frame(ages.HIV %>% group_by(bins) %>% dplyr::summarise(total=as.numeric(n())))
row.names(Bygroup.total) <- Bygroup.total$bins
df <- cbind(Bygroup.count, Bygroup.total[,2])
colnames(df)[2] <- "total"
df$prop <- df$n/df$total
df <- as.data.frame(df)
rownames(df) <- Bygroup.total$bins
df$bin <- rownames(df)
df$bin <- factor(df$bin, levels = df$bin)
df$lab <- paste0(df$n,"(",round(df$prop,2),")")

ggplot(df, aes(x = bin, y = prop, group = 1)) +
  geom_bar(stat='identity', fill = "red") +
  geom_text(aes(x=bin,y=prop,label=lab),vjust=-1) +
  ggtitle("All Passed") +
  ylim(c(0, max(df$prop)*1.1)) +
  xlab("Age") +
  ylab("Variants per sample with Any Mutation") +
  theme(plot.title = element_text(hjust = 0.5))










######


