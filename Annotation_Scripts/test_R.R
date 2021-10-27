# test100 <- read.table("/Users/brian/Bolton/UKBB/results/12/update/transplant_match/match/1232933_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "")
test <- read.table("/Users/brian/Bolton/UKBB/results/10/1000294_23153_0_0.final.tsv", 
                   header = T, sep = "\t", comment.char = "", quote="")
test1 <- read.table("/Users/brian/Bolton/UKBB/results/lung/update/5993871_23153_0_0.final.tsv", 
                   header = T, sep = "\t", comment.char = "", quote="")
test2 <- read.table("/Users/brian/Bolton/UKBB/results/lung/update/1218265_23153_0_0.final.tsv", 
                    header = T, sep = "\t", comment.char = "", quote="")

# test <- read.table("~/Bolton/UKBB/results/12/combined.UKBB.12.tsv", header = T, sep = "\t", comment.char = "", quote="")
# test <- read.table("~/Bolton/UKBB/results/10/1000043_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "", quote="")
# test <- read.table("~/Bolton/UKBB/results/10/1007298_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "", quote="")
# test <- read.table("~/Bolton/UKBB/results/10/1006140_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "", quote="")
# test=test[-775,]
# test2 <- read.table("~/Bolton/data/final_all_rows_all_filters.tsv", sep = "\t", header = T, quote = "")
# test <- read.table("~/Bolton/UKBB/docs/DNAnexus/apps/test/1227284_23153_0_0.final.tsv", sep = "\t", header = T, quote = "", comment.char = "")
# 
# test <- read.table("/Users/brian/Bolton/UKBB/results/12/update/transplant_match/match/1230132_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "")

library(tidyverse)
mutect <- read.table("~/Bolton/UKBB/docs/DNAnexus/apps/test/mutect.tsv", header = F, sep = "\t", comment.char = "", quote="") 
vardict <- read.table("~/Bolton/UKBB/docs/DNAnexus/apps/test/vardict.tsv", header = F, sep = "\t", comment.char = "", quote="")
test <- read.table("~/Bolton/UKBB/results/59/5977335_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "", quote="")
colnames(mutect) <- c("CHROM", "POS", "REF", "ALT", "FILTER", "AD")
colnames(vardict) <- c("CHROM", "POS", "REF", "ALT", "FILTER", "AD")
mutect <- mutect %>% separate(AD, c("gt_AD_ref","gt_AD_alt"), 
                              sep=",",extra = "merge", fill = "right")
vardict <- vardict %>% separate(AD, c("gt_AD_ref","gt_AD_alt"), 
                                sep=",",extra = "merge", fill = "right")


((FMT/AF * FMT/DP < 6) && ((FMT/MQ < 55.0 && FMT/NM > 1.0) || (FMT/MQ < 60.0 && FMT/NM > 2.0) || (FMT/DP < 10) || (FMT/QUAL < 45)))

B <- read.table("/Users/brian/Bolton/data/ASXL1.txt")
test <- read.table("/Users/brian/Downloads/5977335_23153_0_0.final.tsv", header = T, sep = "\t", comment.char = "", quote="")






test$Mutect2_gt_alt_fwd[1]/(test$Mutect2_gt_alt_fwd+test$Mutect2_gt_alt_rev)[1]<.1


final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions %>%
  mutate(Mutect2_SB2 = 
           ifelse(Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
                    Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,0,
                  ifelse((Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1,0,1)),
  )


test <- test %>%
  mutate(Mutect2_SB_test1 = Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev),
         Mutect2_SB_test2 = (Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev))


test <- test %>%
  mutate(Mutect2_SB_test3 = Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
           Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,
         Mutect2_SB_test4 = (Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1)

test <- test %>%
  mutate(Mutect2_gt_alt_fwd2 = Mutect2_gt_alt_fwd,
         Mutect2_gt_alt_rev2 = Mutect2_gt_alt_rev)


test <- test %>%
  mutate(Mutect2_SB2 = 
           ifelse(Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
                    Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,0,
                  ifelse((Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1,0,1))
  )


test2=final.coding.gnomad.sorted.regions
test2 <- test2 %>%
  mutate(Mutect2_SB_test1 = Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev),
         Mutect2_SB_test2 = (Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev))


test2 <- test2 %>%
  mutate(Mutect2_SB_test3 = Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
           Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,
         Mutect2_SB_test4 = (Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1)

test2 <- test2 %>%
  mutate(Mutect2_gt_alt_fwd2 = Mutect2_gt_alt_fwd,
         Mutect2_gt_alt_rev2 = Mutect2_gt_alt_rev)


test2 <- test2 %>%
  mutate(Mutect2_SB2 = 
           ifelse(Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
                    Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,0,
                  ifelse((Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1,0,1))
  )

normalized.read.counts <- readRDS(file = "/Users/brian/Bolton/CCOC/CCOC.RNAseq.normalized.rds")


library(vcfR)
library(tidyverse)
fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}
fillblank <- function(x, r) {
  x[nchar(x)==0] = r
  return(x)
}
mutect2 <- vcfR2tidy(read.vcfR("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/1230132_23153_0_0.vep.annotated.vcf.gz", verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
vardict.df = vcfR2tidy(read.vcfR("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/complex.only.variant.query.vep.annotated.vcf.gz", verbose = FALSE),
                 single_frame = TRUE,
                 info_types = TRUE,
                 format_types = TRUE)$dat
mutect2.info <- mutect2$meta[mutect2$meta$Tag=="INFO",]$ID
mutect2.info <- mutect2.info[!(mutect2.info %in% c("SAMPLE","PON_RefDepth", "PON_AltDepth", "PON_FISHER", "CSQ"))]
mutect.df <- mutect2[["dat"]]

## new VEP has 96 fields
string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames <- str_split(string, "\\|")[[1]]
# CSQ has 91 columns, maybe we should suffix cols with VEP?
mutect.df <- mutect.df %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
mutect.df$gnomAD_AF_VEP <- fillblank(mutect.df$gnomAD_AF_VEP, NA)
vardict.df <- vardict.df %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")

# mutect2_vep3 <- vcfR2tidy(read.vcfR("~/Bolton/UKBB/results/mutect_vep.test.vcf.gz", verbose = FALSE),
#                      single_frame = TRUE,
#                      info_types = TRUE,
#                      format_types = TRUE)
# mutect2.info <- mutect2$meta[mutect2$meta$Tag=="INFO",]$ID
# mutect2.info <- mutect2.info[!(mutect2.info %in% c("SAMPLE","PON_RefDepth", "PON_AltDepth", "PON_FISHER", "CSQ"))]
# mutect.df <- mutect2[["dat"]]
# 
# ## new VEP has 96 fields
# string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
# CSQnames <- str_split(string, "\\|")[[1]]
# # CSQ has 91 columns, maybe we should suffix cols with VEP?
# mutect.df <- mutect.df %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
# mutect.df$gnomAD_AF_VEP <- fillblank(mutect.df$gnomAD_AF_VEP, NA)










# apply(mutect.df.test)
# 
# fillblank(mutect.df$gnomAD_AF_VEP, NA)[1:20]
# mutect.df$gnomAD_AF_VEP[1:20]
# mutect.df.test = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
#   filter(is.na(gnomAD_AF_VEP) & nchar(Existing_variation_VEP)>0) %>% 
#   select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),starts_with("gnomAD"), SAMPLE)
# 
# mutect.df.test = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
#   filter(is.na(gnomAD_AF_VEP) & nchar(Existing_variation_VEP)>0) %>% 
#   select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),starts_with("gnomAD"), SAMPLE)
# mutect.df.test = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
#   filter(is.na(gnomAD_AF_VEP) & nchar(Existing_variation_VEP)>0 & nchar(REF)==1 & nchar(ALT==1)) %>% 
#   select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),starts_with("gnomAD"), SAMPLE)
# > mutect.df.test = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
#   +   filter(is.na(gnomAD_AF_VEP) & nchar(Existing_variation_VEP)>0 & nchar(REF)==1 & nchar(ALT==1) & grepl("^rs",Existing_variation_VEP)) %>% 
#   +   select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),starts_with("gnomAD"), SAMPLE)
mutect.df.test2 = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
   filter(is.na(gnomAD_AF_VEP) & nchar(Existing_variation_VEP)>0 & nchar(REF)==1 & nchar(ALT==1) & grepl("^rs",Existing_variation_VEP)) %>% 
   select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),starts_with("gnomAD"), SAMPLE)

mutect.df.test3 = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
  filter(is.na(gnomADe_AF_VEP) & nchar(Existing_variation_VEP)>0 & nchar(REF)==1 & nchar(ALT==1) & grepl("^rs",Existing_variation_VEP)) %>% 
  select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),(starts_with("gnomADe") & ends_with("VEP")), SAMPLE)

mutect.df.test4 = as.data.frame(apply(mutect.df, 2, fillblank, NA)) %>% 
  filter(is.na(gnomAD_AF_VEP) & is.na(gnomADe_AF_VEP)) %>% 
  select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),(starts_with("gnomADe") & ends_with("VEP")), SAMPLE)

         
         
         
m2 <- vcfR2tidy(read.vcfR("~/Bolton/UKBB/results/1000144_23153_0_0.vep.annotated.vcf.gz", verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
m2.info <- m2$meta[m2$meta$Tag=="INFO",]$ID
m2.info <- m2.info[!(m2.info %in% c("SAMPLE","PON_RefDepth", "PON_AltDepth", "PON_FISHER", "CSQ"))]
m2.df <- m2[["dat"]]

## new VEP has 96 fields
string2 <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames2 <- str_split(string2, "\\|")[[1]]
# CSQ has 91 columns, maybe we should suffix cols with VEP?
m2.df <- m2.df %>% separate(CSQ, paste0(CSQnames2, "_VEP"), sep="\\|", extra = "merge", fill = "right")
m2.df$gnomAD_AF_VEP <- fillblank(mutect.df$gnomAD_AF_VEP, NA)

m2.df.test2 = as.data.frame(apply(m2.df, 2, fillblank, NA)) %>% 
  filter(is.na(gnomAD_AF_VEP) & is.na(gnomADe_AF_VEP) & is.na(gnomADg_AF_VEP)) %>% 
  select(CHROM,POS,REF,ALT,starts_with("Existing_variation"),starts_with("gnomAD"), SAMPLE)

##################
my.list=list("Frame_Shift_Del","Frame_Shift_Ins","Nonsense_Mutation","Splice_Site")

A.filter1 <- A %>% filter((grepl("truncat", source) | grepl("Frame_Shift", VariantClass) |
                             grepl("Nonsense_Mutation", VariantClass) | grepl("Splice_Site", VariantClass)) &
                            Exon=="*")

A.filter2 <- A %>% filter((grepl("truncat", source) | grepl("Frame_Shift", VariantClass) |
                            grepl("Nonsense_Mutation", VariantClass) | grepl("Splice_Site", VariantClass)) &
                           Exon=="*" & (
                             VariantClass!="Frame_Shift_Del" &
                               VariantClass!="Frame_Shift_Ins" &
                               VariantClass!="Nonsense_Mutation" &
                               VariantClass!="Splice_Site"
                           ))
# 292 rows
# should be 516 rows
A.filter3 <- A %>% filter((grepl("truncat", source) | grepl("Frame_Shift", VariantClass) |
                             grepl("Nonsense_Mutation", VariantClass) | grepl("Splice_Site", VariantClass)) &
                            Exon=="*" & (
                              VariantClass=="Frame_Shift_Del" |
                                VariantClass=="Frame_Shift_Ins" |
                                VariantClass=="Nonsense_Mutation" |
                                VariantClass=="Splice_Site"
                            ))
A.filter4 <- dplyr::setdiff(A, A.filter3)

dim(A)
1079-292
## need 787 not truncating + 516 truncating for 129 genes




final.A <- data.frame()
Genes=unique(A.filter3$Gene)

for (gene in Genes) {
  temp <- A.filter3 %>% filter(Gene==gene)
  for (vc in my.list) {
    vc
    if (!(vc %in% temp$VariantClass)) {
      row <- temp[1,]
      row$VariantClass <- vc
      temp <- rbind(temp, row)
    }
  }
  final.A <- rbind(final.A, temp)
}
table(final.A$VariantClass)

A.done <- rbind(A.filter4, final.A)
A.done <- A.done[with(A.done, order(Gene, VariantClass)), ]

n1=length(unique(A$Gene))
n2=length(unique(A.done$Gene))
length(unique(A.not.A.filter3$Gene))

write.table(A.done[,1:9], "~/Bolton/data/pd_table_kbreview_bick_trunc.tsv", sep = "\t", quote = F, row.names = F)
xlsx::write.xlsx(A.done[,1:9], "~/Bolton/data/pd_table_kbreview_bick_trunc.xlsx", row.names = F)








library(vcfR)
library(tidyverse)
pm <- vcfR2tidy(read.vcfR("/Users/brian/Bolton/data/panmyeloid_variant_counts.vep.annotated.vcf.gz", verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
                                           
pm <- read.table("/Users/brian/Bolton/data/panmyeloid_variant_counts.vep.tsv", header = T, sep = "\t", comment.char = "", quote="")
colnames(pm)[1] <- "CHROM"
pm <- pm[,-3]
panmyeloid_variant_counts = read_delim(args$pan_myeloid, delim = "\t", guess_max = 5e6)

row.names(pm) <- with(pm[,c(1:4)],paste(gsub("chr","",CHROM),POS,REF,ALT,sep="_"))
pm$ID_VARIANT = with(pm[,c(1:4)],paste(gsub("chr","",CHROM),POS,REF,ALT,sep="_"))
pm2 <- pm[panmyeloid_variant_counts$ID_VARIANT,]
colnames(pm2)[7] <- "CSQ"


string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINE"
CSQnames <- str_split(string, "\\|")[[1]]
# CSQ has 91 columns, maybe we should suffix cols with VEP?
pm2 <- pm2 %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
pm3 <- cbind(pm2[,c(1:4)], panmyeloid_variant_counts,pm2[,c(8:114)])
write.table(pm3, "/Users/brian/Bolton/data/panmyeloid_variant_counts.vep.annotated.vcf.tsv", sep = "\t", row.names = F, col.names = T, quote = F)


