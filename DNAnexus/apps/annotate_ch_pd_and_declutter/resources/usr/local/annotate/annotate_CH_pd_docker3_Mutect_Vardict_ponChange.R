#!/usr/local/bin/Rscript

library(argparser)
library(vcfR)
library(tidyverse)
library(stringr)
# library(readxl)
library(VariantAnnotation)
library(ggsci)
library(data.table)
library(httr)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(R453Plus1Toolbox)
library(sqldf)
library(jsonlite)

'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}


cosmic.run <- function(df, sample.column) {

  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  X.1 <- parallel::mclapply(ls(X1), function(x) {
  # X.1 <- lapply(ls(X1)[7:8], function(x) {
    df.chr <- X1[[x]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP)
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    cosmic <- read.table(paste0("/home/dnanexus/CosmicMutantExport.final.more.minimal.test.",x,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
    colnames(cosmic) <- c("COSMIC_ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")

    cosmic.test <- sqldf("SELECT l.*, r.COSMIC_ID, r.var_key as var_key_cosmic, HGVSP, r.CosmicCount,
              r.heme_cosmic_count, r.myeloid_cosmic_count
              FROM `df.chr` as l
              LEFT JOIN `cosmic` as r
              on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
    cosmic.test <- cosmic.test[,c(colnames(df.chr),"COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count")]
    cosmic.test$CosmicCount <- fillna(cosmic.test$CosmicCount,0)
    cosmic.test$heme_cosmic_count <- fillna(cosmic.test$heme_cosmic_count,0)
    cosmic.test$myeloid_cosmic_count <- fillna(cosmic.test$myeloid_cosmic_count,0)
    cosmic.test <- cosmic.test %>%
      group_by(var_key,cosmic.test[[sample.column]]) %>%
      slice_max(order_by = heme_cosmic_count, n = 1, with_ties = F)
    cosmic.test <- data.frame(cosmic.test)

    row.names(cosmic.test) <- cosmic.test$var_key
    cosmic.test <- cosmic.test[df.chr$var_key,]

    if (all(cosmic.test[,c("CHROM","POS","REF","ALT",sample.column)]==df.chr[,c("CHROM","POS","REF","ALT",sample.column)])) {
      print(paste("good inner", x))
      return(cosmic.test)
    } else {
      print(paste("bad inner", x))
    }
  })
  proc.time() - ptm

  if (length(X1) == length(X.1)) {
    names(X.1) <- ls(X1)
  } else {
    stop("WRONG")
  }

  all.match <- all(sapply(ls(X1), function(x) {
    dim(X1[[x]])[1] == dim(X.1[[x]])[1]
  }))
  if (all.match) {
    cosmic.final2 <- bind_rows(X.1)
  }
  if (all(cosmic.final2[,c("CHROM","POS","REF","ALT",sample.column)]==df[,c("CHROM","POS","REF","ALT",sample.column)])) {
    print("good")
    MUTS <- cbind(df, cosmic.final2 %>% dplyr::select(COSMIC_ID,CosmicCount,heme_cosmic_count,myeloid_cosmic_count))
    return(MUTS)
  } else {
    print("bad")
  }
}

parser <- arg_parser("Process NGS variant calls")
parser <- add_argument(parser, "--impact-file", type="character", help="kelly's file")
parser <- add_argument(parser, "--bick-file", type="character", help="bick's file")
parser <- add_argument(parser, "--cosmic-file", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--variant-calls-files", type="character", help="variant caller files")
parser <- add_argument(parser, "--TSG-file", type="character", help="tumor suppressor genes")
parser <- add_argument(parser, "--eid", type="character", help="eid for output")
parser <- add_argument(parser, "--oncoKB-curated", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--pd-annotation-file", type="character", help="PD annoation files from papaemmanuil lab") # updated
parser <- add_argument(parser, "--pan-myeloid", type="character", help="panmyeloid_variant_counts.tsv") # updated
parser <- add_argument(parser, "--number-samples-annot", type="integer", help="number samples in loci for variant", default=5)
parser <- add_argument(parser, "--blacklist", type="character", help="ENCODE blacklist hg38")
parser <- add_argument(parser, "--segemental-duplications", type="character", help="segmental duplications")
parser <- add_argument(parser, "--simple-repeats", type="character", help="Simple Tandem Repeats by TRF")
parser <- add_argument(parser, "--repeat-masker", type="character", help="repeatMaskerJoinedCurrent")
parser <- add_argument(parser, "--target-length", type="integer", help="total sum of target intervals", default=39000000)
parser <- add_argument(parser, "--bolton-bick-vars", type="character", help="help") # updated
parser <- add_argument(parser, "--mut2-bick", type="character", help="") # updated
parser <- add_argument(parser, "--mut2-kelly", type="character", help="") # updated
parser <- add_argument(parser, "--col-order", type="character", help="") # updated
parser <- add_argument(parser, "--matches2", type="character", help="") # updated MORE

args <- parse_args(parser)


if (is.na(args$impact_file) || is.na(args$bick_file) || is.na(args$cosmic_file) || 
    is.na(args$variant_calls_files) || is.na(args$TSG_file) || is.na(args$eid) || 
    is.na(args$oncoKB_curated) || is.na(args$pd_annotation_file) || 
    is.na(args$pan_myeloid) || is.na(args$target_length) ||
    is.na(args$bolton_bick_vars) || is.na(args$mut2_bick) || is.na(args$mut2_kelly) || is.na(args$matches2)) {
  stop("Missing one of the annotation files or output")
}

variant.files <- str_split(args$variant_calls_files, ",")[[1]]
bf.correction <- 0.05/args$target_length


mutect2_vcf_file <- variant.files[which(grepl("mutect", variant.files))]
vardict_vcf_file <- variant.files[which(grepl("vardict", variant.files))]
       
### read in                        
mutect2 <- vcfR2tidy(read.vcfR(mutect2_vcf_file, verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
mutect2.info <- mutect2$meta[mutect2$meta$Tag=="INFO",]$ID
mutect2.info <- mutect2.info[!(mutect2.info %in% c("SAMPLE","PON_RefDepth", "PON_AltDepth", "PON_FISHER", "CSQ"))]
mutect.df <- mutect2[["dat"]]

vardict <- vcfR2tidy(read.vcfR(vardict_vcf_file, verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
vardict.info <- vardict$meta[vardict$meta$Tag=="INFO",]$ID
vardict.info <- vardict.info[!(vardict.info %in% c("OLD_MULTIALLELIC", "SAMPLE", "CSQ"))]
vardict.df <- vardict[["dat"]]
vardict.df <- vardict.df[,!colnames(vardict.df) %in% c("OLD_MULTIALLELIC")]


mutect.df <- mutect.df[!duplicated(mutect.df),]
vardict.df <- vardict.df[!duplicated(vardict.df),]

# test <-  mutect.df[is.na(mutect.df$CSQ),]
## remove blank CSQ, OTA non-coding now doing this after merge since vardict depends on mutect
#mutect.df <- mutect.df[!is.na(mutect.df$CSQ),]
#vardict.df <- vardict.df[!is.na(vardict.df$CSQ),]

mutect.df <- mutect.df %>% separate(gt_AD, c("gt_AD_ref","gt_AD_alt"), 
                                    sep=",",extra = "merge", fill = "right")
vardict.df <- vardict.df %>% separate(gt_AD, c("gt_AD_ref","gt_AD_alt"), 
                                      sep=",",extra = "merge", fill = "right")


#covert AD from character to numeric
mutect.df[,c("gt_AD_ref","gt_AD_alt")] <- lapply(mutect.df[,c("gt_AD_ref","gt_AD_alt")], 
                                                 function(x) as.numeric(as.character(x)))
vardict.df[,c("gt_AD_ref","gt_AD_alt")] <- lapply(vardict.df[,c("gt_AD_ref","gt_AD_alt")], 
                                                  function(x) as.numeric(as.character(x)))


# Prefix all FORMAT/ (gt_) columns with caller such as Vardict_gt_AD_ref & Vardict_gt_AD_alt
colnames(mutect.df)[which(grepl("^gt_", colnames(mutect.df)))] <-
  paste0("Mutect2_", colnames(mutect.df)[which(grepl("^gt_", colnames(mutect.df)))])
# then move all FILTER & FORMAT columns to front after ALT/FILTER column
mutect2.pfx <- which(grepl("Mutect2_", colnames(mutect.df)))
mutect.new.col.ord <- c(1:7,mutect2.pfx,8:(min(mutect2.pfx)-1))
length(mutect.new.col.ord) == dim(mutect.df)[2]
mutect2.final <- mutect.df[,mutect.new.col.ord] 
colnames(mutect2.final)[colnames(mutect2.final)=="FILTER"] = "Mutect2_FILTER"
mutect2.final$Mutect2_PASS <- ifelse(grepl("PASS", mutect2.final$Mutect2_FILTER), 1, 0) ## 56 cols

## prefix FMT data with Caller
colnames(vardict.df)[which(grepl("^gt_", colnames(vardict.df)))] <-
  paste0("Vardict_", colnames(vardict.df)[which(grepl("^gt_", colnames(vardict.df)))])
vardict.pfx <- which(grepl("Vardict_", colnames(vardict.df)))
vardict.new.col.ord <- c(1:7,vardict.pfx,8:(min(vardict.pfx)-1))
length(vardict.new.col.ord) == dim(vardict.df)[2]
vardict.final <- vardict.df[,vardict.new.col.ord] 
colnames(vardict.final)[colnames(vardict.final)=="FILTER"] = "Vardict_FILTER"
vardict.final$Vardict_PASS <- ifelse(grepl("PASS", vardict.final$Vardict_FILTER), 1, 0) ## 52 cols



## first fill PON_2AT2 percent column
mutect2.final$PON_2AT2_percent <- fillna(mutect2.final$PON_2AT2_percent, 0)
colnames(mutect2.final)[colnames(mutect2.final)=="PON_2AT2_percent"] = "Mutect2_PON_2AT2_percent"

vardict.final$PON_2AT2_percent <- fillna(vardict.final$PON_2AT2_percent, 0)
colnames(vardict.final)[colnames(vardict.final)=="PON_2AT2_percent"] = "Vardict_PON_2AT2_percent"

## remove columns where all NA
## only do this for columns that are not from Alex's script
na.mutect.cols <- names(mutect2.final[,apply(mutect2.final, 2, function(x) all(is.na(x)))])
na.vardict.cols <- names(vardict.final[,apply(vardict.final, 2, function(x) all(is.na(x)))])

length(na.mutect.cols)
length(na.vardict.cols)
## only do this for columns that are not from Alex's script
na.vardict.cols <- setdiff(na.vardict.cols, c("calpos","calref","calalt","calpos_best","calref_best","calalt_best"))
length(na.vardict.cols)
mutect2.final <- mutect2.final[,-which(names(mutect2.final) %in% c(na.mutect.cols))]
vardict.final <- vardict.final[,-which(names(vardict.final) %in% c(na.vardict.cols))]




########################################################################################

# mutect2.final.driver <- dplyr::left_join(mutect2.final, vars,
#                                          by = c("CHROM"="CHROM", "POS"="POS", 
#                                                 "REF"="REF", "ALT"="ALT"))
## prefix INFO data with Caller
colnames(mutect2.final)[colnames(mutect2.final) %in% mutect2.info] <-
  paste0("Mutect2_", colnames(mutect2.final)[colnames(mutect2.final) %in% mutect2.info])

# vardict.final.driver <- left_join(vardict.final, vars,
#                                   by = c("CHROM"="CHROM", "POS"="POS", 
#                                          "REF"="REF", "ALT"="ALT"))
colnames(vardict.final)[colnames(vardict.final) %in% vardict.info] <-
  paste0("Vardict_", colnames(vardict.final)[colnames(vardict.final) %in% vardict.info])



## Way 3 works
#https://stackoverflow.com/questions/16042380/merge-data-frames-and-overwrite-values
colnames(mutect2.final)[which(colnames(mutect2.final)=="CALLER")] <- "Mutect2_CALLER"
colnames(vardict.final)[which(colnames(vardict.final)=="CALLER")] <- "Vardict_CALLER"
mutect2.final$Mutect2_CALLER <- 1
vardict.final$Vardict_CALLER <- 1

write.table(mutect2.final, "~/mutect.tmp", row.names = F, quote = F, sep = "\t")
mutect2.final <- read.table("~/mutect.tmp", sep = "\t", header = T, quote = "")
write.table(vardict.final, "~/vardict.tmp", row.names = F, quote = F, sep = "\t")
vardict.final <- read.table("~/vardict.tmp", sep = "\t", header = T, quote = "")


## vardict with mutect2
intersection <- intersect(colnames(vardict.final), colnames(mutect2.final))
intersection <- intersection[6:length(intersection)] ## first 5 are the grouped by columns
intersection.cols.x <- paste0(intersection, ".x")
intersection.cols.y <- paste0(intersection, ".y")
final <- full_join(vardict.final, mutect2.final,
                   by = c("CHROM"="CHROM", "POS"="POS", 
                          "REF"="REF", "ALT"="ALT", "SAMPLE"="SAMPLE"))
final[,intersection.cols.x][is.na(final[,intersection.cols.x])] <- final[,intersection.cols.y][is.na(final[,intersection.cols.x])]
#final <- final[!is.na(final$CSQ),]

## remove *.x from column names
colnames(final)[colnames(final) %in% intersection.cols.x] <- intersection
## remove *.y columns completely
final <- final[,!(colnames(final) %in% intersection.cols.y)]

## remove blank CSQ columns because they are non-coding and separate
final <- final[!is.na(final$CSQ),]
# VEP CSQ
## new VEP has 96 fields
## new VEP with gnomADg has 106 fields
string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames <- str_split(string, "\\|")[[1]]
# CSQ has 91 columns, maybe we should suffix cols with VEP?
final <- final %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
############################################################
##
############################################################
final$PON_RefDepth <- as.numeric(final$PON_RefDepth)
final$PON_AltDepth <- as.numeric(final$PON_AltDepth)

final[,c("PON_AltDepth", "PON_FISHER", "PON_RefDepth")] <- sapply(final[,c("PON_AltDepth", "PON_FISHER", "PON_RefDepth")], as.numeric)
## fill NA with 0
final[, c("Mutect2_CALLER","Vardict_CALLER")][is.na(final[, c("Mutect2_CALLER","Vardict_CALLER")])] <- 0
final[, c("Mutect2_PASS","Vardict_PASS")][is.na(final[, c("Mutect2_PASS","Vardict_PASS")])] <- 0


final.coding.gnomad.sorted <- final[with(final, order(CHROM, POS)), ]

# gnomAD 
gnomAD.col <- grep("gnomAD_.*_VEP", colnames(final.coding.gnomad.sorted))
final.coding.gnomad.sorted[,gnomAD.col] <- apply(final.coding.gnomad.sorted[,gnomAD.col], 2, as.numeric)
## we want to know which are NA
#final.coding.gnomad.sorted[,gnomAD.col][is.na(final.coding.gnomad.sorted[,gnomAD.col])] <- 0
final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP <- as.numeric(apply(final.coding.gnomad.sorted[,gnomAD.col], 1, function(x){
  max(x, na.rm = T)
}))
# gnomADe /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz
gnomADe.col <- grep("gnomADe_.*_VEP", colnames(final.coding.gnomad.sorted))
final.coding.gnomad.sorted[,gnomADe.col] <- apply(final.coding.gnomad.sorted[,gnomADe.col], 2, as.numeric)
## we want to know which are NA
#final.coding.gnomad.sorted[,gnomADe.col][is.na(final.coding.gnomad.sorted[,gnomADe.col])] <- 0
final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP <- as.numeric(apply(final.coding.gnomad.sorted[,gnomADe.col], 1, function(x){
  max(x, na.rm = T)
}))
# gnomADg http://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/r3.0/
gnomADg.col <- grep("gnomADg_.*_VEP", colnames(final.coding.gnomad.sorted))
final.coding.gnomad.sorted[,gnomADg.col] <- apply(final.coding.gnomad.sorted[,gnomADg.col], 2, as.numeric)
## we want to know which are NA
#final.coding.gnomad.sorted[,gnomADg.col][is.na(final.coding.gnomad.sorted[,gnomADg.col])] <- 0
final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP <- as.numeric(apply(final.coding.gnomad.sorted[,gnomADg.col], 1, function(x){
  max(x, na.rm = T)
}))


getVAFs <- function(x) {
  mutect_VAF <- grep("Mutect2_gt_AF", colnames(x))
  vardict_VAF <- grep("Vardict_gt_AF", colnames(x))
  VAFS <- c(mutect_VAF, vardict_VAF)
  VAFS
}
#test <- final.coding.gnomad.sorted[final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP > 0.0007,]
# final.coding.gnomad.sorted$gnomAD_MAX.Stringent.0007 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.0007 &
#                                                            final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP < 0.0007 &
#                                                            final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP < 0.0007)
# final.coding.gnomad.sorted$gnomAD_MAX.lessStringent.0007 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.0007 |
#                                                            final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP < 0.0007 |
#                                                            final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP < 0.0007)
## stringentness is w.r.t keeping variant, any can be > .005 and get filtered out
## if 0.003,0.002,0.001 keep but if 0.006,0.003,0.001 filter out
final.coding.gnomad.sorted$gnomAD_MAX.Stringent.005 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.005 &
                                                          final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP < 0.005 &
                                                          final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP < 0.005)
## if False means all 3 groups has a Max VAF > .005 so this is less stringent for filtering out because 
## all 3 needs to have a Max above .005
## if 0.003,0.002,0.001 keep, if 0.006,0.003,0.001 keep, but if 0.006,0.006,0.006 filter out
## all have to be > .005 to get filtered out
final.coding.gnomad.sorted$gnomAD_MAX.lessStringent.005 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.005 |
                                                              final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP < 0.005 |
                                                              final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP < 0.005)

#########################################################################################################
# blacklist, segmental duplications, simple tandem repeatts, and masked repeat regions
blacklist <- read.table(args$blacklist, header = F, sep = "\t")
dups <- rtracklayer::import(args$segemental_duplications)
repeats <- rtracklayer::import(args$simple_repeats)
repeatMasker <- rtracklayer::import(args$repeat_masker) 
colnames(blacklist) <- c("CHROM", "START", "END")
blacklist$blacklist <- T
dups$dups <- T
repeats$repeats <- T
#repeatMasker$region <- with(as.data.frame(repeatMasker), paste0(seqnames, ":", start, "-", end))
repeatMasker$repeatMasker <- T

blackist.gr <- GRanges(blacklist)
sample.gr <- GRanges(seqnames = final.coding.gnomad.sorted$CHROM,
                  ranges = IRanges(start=final.coding.gnomad.sorted$POS, 
                                   end=final.coding.gnomad.sorted$POS), 
                  mcols = final.coding.gnomad.sorted[,3:ncol(final.coding.gnomad.sorted)])
join <- plyranges::join_overlap_left(sample.gr, blackist.gr)
join2 <- plyranges::join_overlap_left(join, dups)
join3 <- plyranges::join_overlap_left(join2, repeats)
join4 <- plyranges::join_overlap_left(join3, repeatMasker)

columns <- colnames(final.coding.gnomad.sorted)
final.coding.gnomad.sorted.regions <- as.data.frame(join4)[,c(1,2,6:ncol(as.data.frame(join4)))]
colnames(final.coding.gnomad.sorted.regions) <- c(columns, "blacklist", 
                                                           "dup_score","dups",
                                                           "trf","trf_score","simpleTandemRepeat",
                                                           "RMName","RMScore","repeatMasker")
final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions[!duplicated(final.coding.gnomad.sorted.regions[,c("CHROM","POS","REF","ALT","SAMPLE")]),]
final.coding.gnomad.sorted.regions$blacklist[is.na(final.coding.gnomad.sorted.regions$blacklist)] <- F
final.coding.gnomad.sorted.regions$dups[is.na(final.coding.gnomad.sorted.regions$dups)] <- F
final.coding.gnomad.sorted.regions$simpleTandemRepeat[is.na(final.coding.gnomad.sorted.regions$simpleTandemRepeat)] <- F
final.coding.gnomad.sorted.regions$repeatMasker[is.na(final.coding.gnomad.sorted.regions$repeatMasker)] <- F

final.coding.gnomad.sorted.regions$called <- 
  apply(final.coding.gnomad.sorted.regions[,c("Mutect2_CALLER","Vardict_CALLER")],
        1,
        function(x) sum(x)>=2)

final.coding.gnomad.sorted.regions <-
  final.coding.gnomad.sorted.regions %>% separate(Vardict_gt_ALD, c("Vardict_gt_ALD_forw","Vardict_gt_ALD_rev"), 
                                                  sep=",",extra = "merge", fill = "right")

final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions %>% 
  separate(Mutect2_gt_SB, c("Mutect2_gt_ref_fwd","Mutect2_gt_ref_rev",
                            "Mutect2_gt_alt_fwd","Mutect2_gt_alt_rev"), 
           sep=",",extra = "merge", fill = "right")

final.coding.gnomad.sorted.regions$Mutect2_gt_alt_fwd <- as.numeric(final.coding.gnomad.sorted.regions$Mutect2_gt_alt_fwd)
final.coding.gnomad.sorted.regions$Mutect2_gt_alt_rev <- as.numeric(final.coding.gnomad.sorted.regions$Mutect2_gt_alt_rev)
final.coding.gnomad.sorted.regions$Vardict_gt_ALD_forw <- as.numeric(final.coding.gnomad.sorted.regions$Vardict_gt_ALD_forw)
final.coding.gnomad.sorted.regions$Vardict_gt_ALD_rev <- as.numeric(final.coding.gnomad.sorted.regions$Vardict_gt_ALD_rev)

## 1 passes SB and 0 fails SB
final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions %>%
  mutate(Mutect2_SB = 
           ifelse(Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
                    Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,0,
                  ifelse((Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1,0,1)),
         Vardict_SB = 
           ifelse(Vardict_gt_ALD_forw/(Vardict_gt_ALD_forw+Vardict_gt_ALD_rev)<.1 |
                    Vardict_gt_ALD_forw/(Vardict_gt_ALD_forw+Vardict_gt_ALD_rev)>.9,0,
                  ifelse((Vardict_gt_ALD_forw+Vardict_gt_ALD_rev)<1,0,1)),
  )

## Passed by both Mutect2 and vardict
## if add in varscan, passed by Mutect2 and (Vardict or Varscan)
final.coding.gnomad.sorted.regions$passed <-
  apply(final.coding.gnomad.sorted.regions[,c("Mutect2_PASS","Vardict_PASS")],
        1,
        #function(x) (sum(x)>=2 & x["Mutect2_PASS"]))
        function(x) (sum(x)>=2))



################################################
## rename for ease
final.passed <- final.coding.gnomad.sorted.regions
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
VAFS <- getVAFs(final.passed)
final.passed$max.under.0.35 <- apply(final.passed[,VAFS],1,function(x) max(x, na.rm = T) < 0.35)
final.passed$max.over.0.02 <- apply(final.passed[,VAFS],1,function(x) max(x, na.rm = T) >= 0.02)

## Strand Bias 2 or more pass SB
final.passed$alt_strand_counts_min_2_callers <- apply(final.passed[,c("Mutect2_SB","Vardict_SB")],
                                                      1,
                                                      function(x) sum(x, na.rm = T)>=2)
## Strand Bias 2 or more pass SB
final.passed$alt_strand_counts_min_1_caller_only <- apply(final.passed[,c("Mutect2_SB","Vardict_SB")],
                                                          1,
                                                          function(x) sum(x, na.rm = T)>=1)
################################################
### Annotate PD

AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               '%3D'='=', '='='=')

## ch_pd stuff
vars <- read.table(args$bolton_bick_vars, sep = "\t", header = T, comment.char = "")
vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")
topmed.mutation.2 <- read.table(args$mut2_bick, sep = "\t", header = T, comment.char = "")
kelly.mutation.2 <- read.table(args$mut2_kelly, sep = "\t", header = T, comment.char = "")
##both Bick/Bolton
matches.2.c.p <- read.table(args$matches2, sep = "\t", header = T, comment.char = "")

final.passed$AAchange <- gsub("(.*p\\.)(.*)", "\\2", final.passed$HGVSp)
for (i in 1:length(AminoAcids)) {
  final.passed$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], final.passed$AAchange)
}
final.passed$gene_loci_p <- paste(final.passed$SYMBOL_VEP,
                                  paste0(sapply(final.passed$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1]),
                                         as.numeric(str_extract(final.passed$AAchange, "\\d+"))),
                                  sep = "_")
final.passed$gene_loci_c <- paste(final.passed$SYMBOL_VEP,
                                  gsub(".*:", "", final.passed$HGVSc),
                                  sep = "_")
final.passed$gene_loci_vep <- ifelse(is.na(final.passed$gene_loci_p),final.passed$gene_loci_c,final.passed$gene_loci_p)
final.passed$key <- with(final.passed, paste(CHROM,POS,REF,ALT,sep = ":"))
final.passed$gene_aachange <- with(final.passed, paste(SYMBOL_VEP, AAchange, sep = "_"))
final.passed$gene_cDNAchange <- paste(final.passed$SYMBOL_VEP, gsub(".*:","",final.passed$HGVSc_VEP), sep="_")

dims <- dim(final.passed)[[1]]
final.passed <- sqldf("SELECT l.*, r.`n.loci.vep`, r.`source.totals.loci`
            FROM `final.passed` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
final.passed <- final.passed[!duplicated(final.passed),]
## make sure aachange exists as in doesn't end with an '_'; example: DNMT3A_ for splice 
final.passed <- sqldf("SELECT l.*, r.`n.HGVSp`, r.`source.totals.p`
            FROM `final.passed` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.gene_aachange = r.gene_aachange) AND r.gene_aachange NOT LIKE '%_'")
final.passed <- final.passed[!duplicated(final.passed),]
final.passed <- sqldf("SELECT l.*, r.`n.HGVSc`, r.`source.totals.c`
            FROM `final.passed` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange) AND r.gene_cDNAchange NOT LIKE '%_'")
final.passed <- final.passed[!duplicated(final.passed),]
dim(final.passed)[[1]] == dims
########################################################################################
annotate.PD <- function(x) {
  
  MUTS <- x
  translate_consequence = c(
    "frameshift_variant" = "Frame_Shift_", 
    "stop_gained" = "Nonsense_Mutation", 
    "splice_region_variant" = "Splice_Site",
    "inframe_deletion" = "In_Frame_Del", 
    "inframe_insertion" = "In_Frame_Ins",
    "synonymous_variant" = "synonymous_variant",
    "missense_variant" = "Missense_Mutation",
    "intron_variant" = "intron_variant",
    "splice_donor_variant" = "Splice_Site",
    "coding_sequence_variant" = "coding_sequence_variant",
    "protein_altering_variant" = "protein_altering_variant",
    "splice_acceptor_variant" = "Splice_Site",
    "3_prime_UTR_variant" = "3_prime_UTR_variant",
    "5_prime_UTR_variant" = "5'Flank",
    "stop_lost" = "_Nonsense_Mutation",
    "start_lost" = "Nonsense_Mutation",
    "stop_retained_variant" = "Nonsense_Mutation")
  
  translate_type = c(
    'SNV' = 'SNV',
    'deletion' = 'Del',
    'insertion' = 'Ins')
  MUTS <- MUTS %>%
    mutate(VariantClass = translate_consequence[str_split(Consequence_VEP,"&",simplify = TRUE)[,1]]) %>%
    mutate(VariantClass = ifelse(VariantClass == 'Frame_Shift_', paste0('Frame_Shift_', translate_type[VARIANT_CLASS_VEP]), VariantClass))
  
  
  MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP","HGVSp_VEP","HGVSc_VEP","VariantClass","SAMPLE","EXON_VEP","AAchange")]
  MUTS$aa_ref <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1])
  MUTS$aa_alt <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][2])
  MUTS$aa_pos <- as.numeric(str_extract(MUTS$AAchange, "\\d+"))
  MUTS$var_key = paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  
  ###create list of tumor suppressor genes####
  ## file 1
  gene_census = read_delim(args$TSG_file, delim = "\t", guess_max = 5e6)
  ## file 2
  oncoKB_curated = read_delim(args$oncoKB_curated, delim = "\t", guess_max = 5e6)
  oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
  
  TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
  
  
  ###########annotate with COSMIC########
  #source("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/cosmic.parallel.compute1.R")
  ## need to add a parallel
  MUTS <- cosmic.run(MUTS, "SAMPLE")
  MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  
  
  #annotate oncokb
  get_oncokb = function(mut) {
    
    apiKey = "a83627cd-47e4-4be0-82dd-8f4cc9e4d6d0"
    request_url = paste0("https://oncokb.org//api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
                         mut['SYMBOL_VEP'], "&alteration=", mut['AAchange'], "&consequence=", mut['VariantClass'])
    
    res=httr::content(httr::GET(request_url, httr::add_headers(Authorization = paste("Bearer", apiKey))))
    return(res$oncogenic)
  }
  
  message("annotating variants with oncoKB...")
  h = curl::new_handle()
  curl::handle_setopt(h, http_version = 2)
  httr::set_config(httr::config(http_version = 0))
  cl = parallel::makeCluster(8)
  MUTS$oncoKB = parallel::parApply(cl, MUTS, 1, get_oncokb)
  parallel::stopCluster(cl)
  
  
  ##annotate myeloid-PD
  
  # Variant Annotation
  # Variants were annotated according to evidence for functional relevance in cancer (putative driver or CH-PD)
  # and for relevance to myeloid neoplasms specifically (CH-Myeloid-PD).
  # We annotated variants as oncogenic in myeloid disease (CH-Myeloid-PD) if they fulfilled any of the following criteria:
  #
  #   1. Truncating variants in NF1, DNMT3A, TET2, IKZF1, RAD21, WT1, KMT2D, SH2B3, TP53, CEBPA, ASXL1, RUNX1, BCOR, KDM6A, STAG2, PHF6, KMT2C, PPM1D, ATM, ARID1A, ARID2, ASXL2, CHEK2, CREBBP, ETV6, EZH2, FBXW7, MGA, MPL, RB1,SETD2, SUZ12, ZRSR2 or in CALR exon 9
  #   2. Translation start site mutations in SH2B3
  #   3. TERT promoter mutations
  #   4. FLT3-ITDs
  #   5. In-frame indels in CALR, CEBPA, CHEK2, ETV6, EZH2
  
  ##read in PD annoation files from papaemmanuil lab##
  ## file 4
  A <- read.table(args$pd_annotation_file, sep = "\t", header = T, quote = "")[,1:9]
  A$keys = apply(A, 1, function(x) {
    names(x)[x != '*'] %>% .[!. %in% c('source', 'ch_my_pd', 'ch_pancan_pd')]
  })
  pd_dfs = split(A, paste(A$keys))
  
  ch_my_variants = NULL
  pd_dfs = split(A, paste(A$keys))
  MUTS$aa_pos <- as.character(MUTS$aa_pos)
  MUTS$Exon <- sapply(MUTS$EXON_VEP, function(x) str_split(x,"/")[[1]][1])
  colnames(MUTS)[5] <- c("Gene")
  MUTS.temp <- MUTS
  matched = list()
  
  for (i in 1:length(pd_dfs)){
    temp = as.data.frame(pd_dfs[i])
    colnames(temp) = colnames(A)
    keys = unlist(unique(temp$keys))
    if ('Exon' %in% keys) {
      MUTS.temp <- MUTS.temp %>% filter(!(var_key %in% matched))
      temp <- temp %>% filter(Exon != '*')
    } else {
      MUTS.temp <- MUTS.temp %>% filter(!(var_key %in% matched))
      temp <- temp %>% filter(Exon == '*')
    }
    # are all keys is both tables?
    print(all(keys %in% colnames(temp), keys %in% colnames(MUTS.temp)))
    pds = left_join(MUTS.temp, temp[,c(keys,'ch_my_pd','source')], by=keys)
    pds = pds %>% filter(ch_my_pd>0)
    matched = append(matched, pds$var_key)
    ch_my_variants = rbind(ch_my_variants, pds)
  }
  
  ch_my_variants = ch_my_variants %>% dplyr::select(source, var_key, ch_my_pd) %>% unique()
  tryCatch(
    expr = {
      ch_my_variants <- aggregate(source ~ var_key + ch_my_pd, data = ch_my_variants, FUN = paste, collapse = ",")
    },
    error = function(e){ 
      message(e)
    }
  )
  ## easier to just not return anything from TC and call this below it
  MUTS = left_join(MUTS, ch_my_variants, by="var_key")
  MUTS <- MUTS %>% 
    mutate(WHY_CH_ch_my_pd = ifelse(MUTS$ch_my_pd>0, "PD_table;", ""))
  
  
  # 6. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times
  
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, paste0(WHY_CH_ch_my_pd, " Cosmic_heme;"), WHY_CH_ch_my_pd))
  
  # 7. Any variant noted as potentially oncogenic in an in-house dataset of 7,000 individuals with myeloid neoplasm greater than or equal to 5 times
  panmyeloid_variant_counts = read_delim(args$pan_myeloid, delim = "\t", guess_max = 5e6)
  panmyeloid_variant_counts$var_key = with(panmyeloid_variant_counts, paste(CHROM,POS,REF,ALT, sep = ":"))
  panmyeloid_variant_counts <- panmyeloid_variant_counts[panmyeloid_variant_counts$Annotation == "ONCOGENIC" | 
                                                           panmyeloid_variant_counts$Annotation == "SNP",]
  
  MUTS$Gene_HGVSp_VEP <- paste(MUTS$Gene, gsub(".*p.","", MUTS$HGVSp_VEP), sep = "_")
  panmyeloid_variant_counts$Gene_HGVSp_VEP <- paste(panmyeloid_variant_counts$SYMBOL_VEP, gsub(".*p.","", panmyeloid_variant_counts$HGVSp_VEP), sep = "_")
  MUTS <- sqldf("SELECT l.*, r.MDS, r.AML, r.MPN, r.Annotation
            FROM `MUTS` as l
            LEFT JOIN `panmyeloid_variant_counts` as r
            on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
  MUTS <- MUTS %>%
    mutate(
      MDS = fillna(MDS, 0),
      AML = fillna(AML, 0),
      MPN = fillna(MPN, 0),
      Annotation = fillna(Annotation, ''),
      n_panmyeloid = ifelse(Annotation == 'ONCOGENIC', MDS + AML + MPN, 0))
  
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(n_panmyeloid>=5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(MUTS$n_panmyeloid>=5, paste0(WHY_CH_ch_my_pd, " Pan_Myeloid;"), WHY_CH_ch_my_pd))
  
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", 1, 0),
           WHY_CH_ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", "OncoKB;", ""))
  MUTS$isOncogenic <- MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic"
  MUTS$isTSG <- MUTS$Gene %in% TSG
  
  # 2. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  # Genes not listed in the cancer census or OncoKB were reviewed in the literature to determine if they were potentially tumor suppressor genes.# Annotate PD based on prevelance in cancer databases (COSMIC, Oncokb)
  
  
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(isTSG & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"), 1, ch_pd),
           WHY_CH_ch_pd = ifelse(isTSG & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"), 
                                 paste0(WHY_CH_ch_pd, " TSG & Truncating;"), WHY_CH_ch_pd))
  
  #3. Any variant reported as somatic at least 20 times in COSMIC
  MUTS$CosmicCount <- as.numeric(MUTS$CosmicCount)
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(CosmicCount>=20, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(CosmicCount>=20, paste0(WHY_CH_ch_pd, " Cosmic;"), WHY_CH_ch_pd))
  
  #4. Any variant meeting criteria for CH-Myeloid-PD as above.
  
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, paste0(WHY_CH_ch_pd, " ch_my_pd>=1"), WHY_CH_ch_pd))
  
  
  annotate_PD.R = MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  
  sample=grep("SAMPLE",colnames(x))
  
  if(all(annotate_PD.R[,c(1:4,9)]==x[,c(1:4,sample)])) {
    message("good same dim as orig")
    annotate_PD.R <- annotate_PD.R %>% # annotate_PD.R is the MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
      dplyr::mutate(gene_loci_p = paste(Gene, paste0(aa_ref, aa_pos), sep = "_"),
                    cDNAchange = gsub(".*:","", HGVSc_VEP),
                    gene_loci_c = paste(Gene, cDNAchange, sep = "_"),
                    gene_aachange = paste(Gene, AAchange, sep = "_"),
                    gene_loci_vep = ifelse(is.na(AAchange) | AAchange == "", gene_loci_c, gene_loci_p),
                    ch_pd2 = case_when(ch_pd==1 ~ 1,
                                       (gene_loci_vep %in% vars$gene_loci_vep | 
                                          gene_aachange %in% matches.2.c.p |
                                          gene_loci_c %in% matches.2.c.p) ~ 1,
                                       TRUE ~ 0),
                    WHY_CH_ch_pd2 = ifelse(ch_pd==1, "CH_pd=1;", ""),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_vep %in% vars$gene_loci_vep, paste0(WHY_CH_ch_pd2, " Gene_loci>=5;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_aachange %in% matches.2.c.p$exact_match, paste0(WHY_CH_ch_pd2, " gene_aachange>=2 BB;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_c %in% matches.2.c.p$exact_match, paste0(WHY_CH_ch_pd2, " gene_loci_c>=2 BB;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_aachange %in% topmed.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_aachange>=2 Bick;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_c %in% topmed.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_loci_c>=2 Bick;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_aachange %in% kelly.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_aachange>=2 Bolton;"), WHY_CH_ch_pd2),
                    WHY_CH_ch_pd2 = ifelse(gene_loci_c %in% kelly.mutation.2$exact_match, paste0(WHY_CH_ch_pd2, " gene_loci_c>=2 Bolton;"), WHY_CH_ch_pd2))
    
    
    annotate_PD.R$WHY_CH <- apply(annotate_PD.R, 1, function(x) {
      return(toJSON(list('ch_my_pd'=x["WHY_CH_ch_my_pd"],
                         'ch_pd2'=x["WHY_CH_ch_pd2"],
                         'ch_pd'=x["WHY_CH_ch_pd"])))
    })
    
    annotate_PD.R$ch_my_pd <- fillna(annotate_PD.R$ch_my_pd, 0)
    annotate_PD.R$ch_pd <- fillna(annotate_PD.R$ch_pd, 0)
    
    return(
      annotate_PD.R %>% ## add n_panmyeloid
        dplyr::select("CHROM","POS","REF","ALT","SAMPLE",
                      "COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count",
                      "MDS", "AML","MPN","n_panmyeloid",
                      "isOncogenic","isTSG",
                      "ch_my_pd","ch_pd","ch_pd2",
                      "VariantClass","AAchange","Gene","WHY_CH") 
    )
  } else {
    message("something wrong with dims")
  }
}


final.passed$Mutect2_PON_2AT2_percent <- fillna(final.passed$Mutect2_PON_2AT2_percent, 0)
final.passed$Vardict_PON_2AT2_percent <- fillna(final.passed$Vardict_PON_2AT2_percent, 0)

## let's look at 
    # gnomAD_MAX.Stringent.0007, gnomAD_MAX.lessStringent.0007, 
    # MAX_gnomAD_AF_VEP, MAX_gnomADe_AF_VEP, MAX_gnomADg_AF_VEP
# final.passed$Vardict_PON_2AT2_percent = 1 = fail
final.passed$passed_everything <- (!final.passed$Vardict_PON_2AT2_percent &
                                     !final.passed$Mutect2_PON_2AT2_percent &
                                     final.passed$alt_strand_counts_min_2_callers &
                                     final.passed$max.under.0.35 &
                                     final.passed$max.over.0.02 &
                                     final.passed$passed &
                                     final.passed$complexity_filters &
                                     final.passed$PON_FISHER <= bf.correction)
                                 
final.passed$passed_everything_mutect <- final.passed$passed_everything & final.passed$Mutect2_PASS

annotation <- annotate.PD(final.passed)

final.df <-  left_join(final.passed,
                       annotation,
                       by=c("CHROM","POS","REF","ALT","SAMPLE"))

final.df$eid <- args$eid
# col.order <- read.table(args$col_order)$V1
# final.df <- final.df[,col.order]
write.table(final.df, paste0(args$eid,".final.tsv"), row.names = F, sep = "\t", quote = F)
write.table(colnames(final.df), paste0(args$eid,".columns.txt"), row.names = F, col.names = F)
