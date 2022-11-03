#!/usr/local/bin/Rscript

library(argparser)
library(vcfR)
library(tidyverse)
library(stringr)
library(VariantAnnotation)
library(ggsci)
library(data.table)
library(httr)
library(rtracklayer)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
library(R453Plus1Toolbox)
options(gsubfn.engine = "R")
library(sqldf)
library(jsonlite)

'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}



parser <- arg_parser("Process NGS variant calls")
parser <- add_argument(parser, "--mutect2-vcf-file", type="character", help="mutect vcf")
parser <- add_argument(parser, "--merged", type="character", help="merged tsv")
parser <- add_argument(parser, "--truncating", type="character", help="BB truncating")
parser <- add_argument(parser, "--TSG-file", type="character", help="tumor suppressor genes")
parser <- add_argument(parser, "--sample_id", type="character", help="sample_id/file id for output")  
parser <- add_argument(parser, "--out_folder", type="character", help="out_folder")
parser <- add_argument(parser, "--oncoKB-curated", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--pd-annotation-file", type="character", help="PD annoation files from papaemmanuil lab") # updated
parser <- add_argument(parser, "--blacklist", type="character", help="ENCODE blacklist hg38")
parser <- add_argument(parser, "--segemental-duplications", type="character", help="segmental duplications")
parser <- add_argument(parser, "--simple-repeats", type="character", help="Simple Tandem Repeats by TRF")
parser <- add_argument(parser, "--repeat-masker", type="character", help="repeatMaskerJoinedCurrent")
parser <- add_argument(parser, "--bolton-bick-vars", type="character", help="help") # updated
parser <- add_argument(parser, "--mut2-bick", type="character", help="") # updated
parser <- add_argument(parser, "--mut2-kelly", type="character", help="") # updated
parser <- add_argument(parser, "--matches2", type="character", help="") # updated MORE

args <- parse_args(parser)
print(args$sample_id)
print(args$merged)


### read in       
mutect2_vcf_file <- args$mutect2_vcf_file
mutect <- vcfR2tidy(read.vcfR(mutect2_vcf_file, verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
mutect.info <- mutect$meta[mutect$meta$Tag=="INFO",]$ID
mutect.info <- mutect.info[!(mutect.info %in% c("SAMPLE","PON_RefDepth", "PON_AltDepth", "PON_FISHER", "CSQ"))]
mutect.df <- mutect[["dat"]]


# test <- mutect2.df %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
mutect.df <- mutect.df[!duplicated(mutect.df),]

print(args$sample_id)


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




## first fill PON_2AT2 percent column
mutect2.final$PON_2AT2_percent <- fillna(mutect2.final$PON_2AT2_percent, 0)
colnames(mutect2.final)[colnames(mutect2.final)=="PON_2AT2_percent"] = "Mutect2_PON_2AT2_percent"


## remove columns where all NA
## only do this for columns that are not from Alex's script
na.mutect.cols <- names(mutect2.final[,apply(mutect2.final, 2, function(x) all(is.na(x)))])
na.mutect.cols <- setdiff(na.mutect.cols, c("Mutect2_gt_PGT","Mutect2_gt_PID","Mutect2_gt_PS","Mutect2_RPA","Mutect2_RU","Mutect2_STRQ"))
length(na.mutect.cols)
mutect2.final <- mutect2.final[,-which(names(mutect2.final) %in% c(na.mutect.cols))]




########################################################################################

# mutect2.final.driver <- dplyr::left_join(mutect2.final, vars,
#                                          by = c("CHROM"="CHROM", "POS"="POS", 
#                                                 "REF"="REF", "ALT"="ALT"))
## prefix INFO data with Caller
colnames(mutect2.final)[colnames(mutect2.final) %in% mutect.info] <-
  paste0("Mutect2_", colnames(mutect2.final)[colnames(mutect2.final) %in% mutect.info])





## Way 3 works
#https://stackoverflow.com/questions/16042380/merge-data-frames-and-overwrite-values
colnames(mutect2.final)[which(colnames(mutect2.final)=="CALLER")] <- "Mutect2_CALLER"
mutect2.final$Mutect2_CALLER <- 1


final <- mutect2.final


## remove blank CSQ columns because they are non-coding and separate
final <- final[!is.na(final$CSQ),]
# VEP CSQ
## new VEP has 96 fields
## new VEP with gnomADg has 106 fields
string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|VAR_SYNONYMS|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames <- str_split(string, "\\|")[[1]]
# mutect2.final <- mutect2.final %>% separate(CSQ, paste0(CSQnamesM, "_VEP"), sep="\\|", extra = "merge", fill = "right")
# if (!is.null(vardict.final$CSvardict.final$CSQ)) {
#   vardict.final <- vardict.final %>% separate(CSQ, paste0(CSQnamesV, "_VEP"), sep="\\|", extra = "merge", fill = "right")
# }
# CSQ has 91 columns, maybe we should suffix cols with VEP?
final <- final %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
############################################################
##
############################################################






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
  VAFS <- c(mutect_VAF)
  VAFS
}


## stringentness is w.r.t keeping variant, any can be > .005 and get filtered out and we don't look at it
## if 0.003,0.002,0.001 keep but if 0.006,0.003,0.001 filter out
final.coding.gnomad.sorted$gnomAD_MAX.Stringent.005 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.005 &
                                                          final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP < 0.005 &
                                                          final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP < 0.005)
#this is same as ...
# final.coding.gnomad.sorted$gnomAD_MAX.Stringent.005 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP >= 0.005 |
#                                                           final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP >= 0.005 |
#                                                           final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP > 0.005)


## if False means all 3 groups has a Max VAF > .005 so this is less stringent for filtering out because 
## all 3 needs to have a Max above .005 which will always be equal to or less than max stringent :)
## if 0.003,0.002,0.001 keep, if 0.006,0.003,0.001 keep, but if 0.006,0.006,0.006 filter out
## all have to be > .005 to get filtered out
final.coding.gnomad.sorted$gnomAD_MAX.lessStringent.005 <- (final.coding.gnomad.sorted$MAX_gnomAD_AF_VEP < 0.005 |
                                                              final.coding.gnomad.sorted$MAX_gnomADe_AF_VEP < 0.005 |
                                                              final.coding.gnomad.sorted$MAX_gnomADg_AF_VEP < 0.005)
final.coding.gnomad.sorted$key <- with(final.coding.gnomad.sorted, paste(CHROM,POS,REF,ALT, sep = ":"))
#########################################################################################################
# blacklist, segmental duplications, simple tandem repeatts, and masked repeat regions
blacklist <- rtracklayer::import(args$blacklist)
print(blacklist)

dups <- rtracklayer::import(args$segemental_duplications)
repeats <- rtracklayer::import(args$simple_repeats)
repeatMasker <- rtracklayer::import(args$repeat_masker) 
#colnames(blacklist) <- c("CHROM", "START", "END")
blacklist$blacklist <- T
print(dups)
dups$dups <- T
repeats$repeats <- T
#repeatMasker$region <- with(as.data.frame(repeatMasker), paste0(seqnames, ":", start, "-", end))
repeatMasker$repeatMasker <- T
#blackist.gr <- GRanges(blacklist)
sample.gr <- GRanges(seqnames = final.coding.gnomad.sorted$CHROM,
                  ranges = IRanges(start=final.coding.gnomad.sorted$POS, 
                                   end=final.coding.gnomad.sorted$POS), 
                  mcols = final.coding.gnomad.sorted[,3:ncol(final.coding.gnomad.sorted)])
join <- plyranges::join_overlap_left(sample.gr, blacklist)
join2 <- plyranges::join_overlap_left(join, dups)
join3 <- plyranges::join_overlap_left(join2, repeats)
join4 <- plyranges::join_overlap_left(join3, repeatMasker)
columns <- colnames(final.coding.gnomad.sorted)
final.coding.gnomad.sorted.regions <- as.data.frame(join4)[,c(1,2,6:ncol(as.data.frame(join4)))]
colnames(final.coding.gnomad.sorted.regions) <- c(columns, "blacklist", 
                                                           "dup_score","dups",
                                                           "trf","trf_score","simpleTandemRepeat",
                                                           "RMName","RMScore","repeatMasker")
final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions[!duplicated(final.coding.gnomad.sorted.regions[,c("CHROM","POS","REF","ALT")]),]
final.coding.gnomad.sorted.regions$key <- with(final.coding.gnomad.sorted.regions, paste(CHROM,POS,REF,ALT, sep = ":"))
final.coding.gnomad.sorted.regions$blacklist[is.na(final.coding.gnomad.sorted.regions$blacklist)] <- F
final.coding.gnomad.sorted.regions$dups[is.na(final.coding.gnomad.sorted.regions$dups)] <- F
final.coding.gnomad.sorted.regions$simpleTandemRepeat[is.na(final.coding.gnomad.sorted.regions$simpleTandemRepeat)] <- F
final.coding.gnomad.sorted.regions$repeatMasker[is.na(final.coding.gnomad.sorted.regions$repeatMasker)] <- F
final.coding.gnomad.sorted.regions$called <- 1

final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions %>% 
  separate(Mutect2_gt_SB, c("Mutect2_gt_ref_fwd","Mutect2_gt_ref_rev",
                            "Mutect2_gt_alt_fwd","Mutect2_gt_alt_rev"), 
           sep=",",extra = "merge", fill = "right")
final.coding.gnomad.sorted.regions$Mutect2_gt_alt_fwd <- as.numeric(final.coding.gnomad.sorted.regions$Mutect2_gt_alt_fwd)
final.coding.gnomad.sorted.regions$Mutect2_gt_alt_rev <- as.numeric(final.coding.gnomad.sorted.regions$Mutect2_gt_alt_rev)


## 1 passes SB and 0 fails SB
final.coding.gnomad.sorted.regions <- final.coding.gnomad.sorted.regions %>%
  mutate(Mutect2_SB = 
           ifelse(Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<.1 |
                    Mutect2_gt_alt_fwd/(Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)>.9,0,
                  ifelse((Mutect2_gt_alt_fwd+Mutect2_gt_alt_rev)<1,0,1))
  )




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
final.passed$max.under.0.35 <- sapply(final.passed[,VAFS],function(x) max(x, na.rm = T) < 0.35)
final.passed$max.over.0.02 <- sapply(final.passed[,VAFS],function(x) max(x, na.rm = T) >= 0.02)


## Strand Bias 2 or more pass SB
final.passed$alt_strand_counts_min_1_caller_only <- sapply(final.passed[,c("Mutect2_SB")],
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

final.passed$AAchange <- gsub("(.*p\\.)(.*)", "\\2", final.passed$HGVSp_VEP)
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
paste0("dims match after sqldf: ",dim(final.passed)[[1]] == dims)

final.passed <- final.passed[!duplicated(final.passed[,c("CHROM","POS","REF","ALT")]),]
########################################################################################
annotate.PD <- function(x) {
  MUTS <- x
  
  ## has 9
  ## missing 6
  # additional we can have "upstream_gene_variant", "non_coding_transcript_variant", "non_coding_transcript_exon_variant"
  # total 18 unique(unlist(sapply(MUTS$Consequence_VEP, function(x) { str_split(x,"&",simplify = TRUE)})))
  translate_consequence = c(
    "frameshift_variant" = "frameshift_variant", #1
    
    "inframe_deletion" = "inframe_deletion", #2
    "inframe_insertion" = "inframe_insertion", #3
    "synonymous_variant" = "synonymous_variant", #4
    "missense_variant" = "missense_variant", #5
    "intron_variant" = "intron_variant", #6
    
    "splice_region_variant" = "splice_region_variant", #7
    "splice_donor_variant" = "splice_region_variant", # missing
    "splice_acceptor_variant" = "splice_region_variant", # missing
    "feature_truncation" = "feature_truncation",
    "3_prime_UTR_variant" = "3_prime_UTR_variant", # missing
    "5_prime_UTR_variant" = "5_prime_UTR_variant", # missing
    
    "stop_gained" = "stop_gained", #8
    "start_lost" = "start_lost", #9
    "stop_lost" = "missense_variant", # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
    "stop_retained_variant" = "synonymous_variant" # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
  )
  
  MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP",
                  "HGVSp_VEP","n.HGVSp","HGVSc_VEP","n.HGVSc",
                  "Consequence_VEP","SAMPLE","EXON_VEP","AAchange")]
  
  MUTS <- MUTS %>%
    mutate(VariantClass = str_split(Consequence_VEP,"&",simplify = TRUE)[,1]) %>%
    mutate(VariantClass2 = case_when((VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) > 0 ~ "inframe_deletion",
                                    (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) < 0 ~ "inframe_insertion",
                                    (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 != 0 ~ "frameshift_variant",
                                    (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF) == nchar(ALT) ~ "missense_variant",
                                    TRUE ~ VariantClass)) %>%
    mutate(VariantClass2 = translate_consequence[VariantClass])
  
  MUTS$aa_ref <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1])
  MUTS$aa_alt <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][2])
  MUTS$aa_pos <- as.numeric(str_extract(MUTS$AAchange, "\\d+"))
  MUTS$var_key = paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  
  ###create list of tumor suppressor genes####
  ## file 1
  gene_census = read.table(args$TSG_file, comment.char = "", sep = "\t", quote = "", header = T)
  #gene_census.og <- gene_census %>% filter(grepl("oncogene", Role_in_Cancer))
  ## file 2
  oncoKB_curated = read.table(args$oncoKB_curated, comment.char = "", sep = "\t", quote = "", header = T)
  oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
  # oncoKB_curated.og = oncoKB_curated %>% filter(Is_Oncogene=="Yes")
  TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
  
  
  
  message("annotating variants with oncoKB...")

  merged <- read.table(args$merged, comment.char = "", sep = "\t", quote = "", header = T)
  COSMIC <- merged %>% dplyr::select(key, COSMIC_ID, CosmicCount, heme_cosmic_count, myeloid_cosmic_count)
  oncoKB <- merged %>% dplyr::select(key, oncoKB)
  MUTS <- MUTS %>% left_join(oncoKB, by=c("var_key"="key"))
  print(table(MUTS$oncoKB))
  message("oncoKB done\n")
  
  ###########annotate with COSMIC########
  ## MUTS MUST BE ORDERED by CHROM (and maybe POS? just to be safe?? Basically whatever column is paralleled needs to be sorted by...)
  message("annotating variants with COSMIC...")
  # MUTS <- cosmic.run(MUTS, "SAMPLE", "var_key")
  MUTS <- MUTS %>% left_join(COSMIC, by=c("var_key"="key"))
  MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  print(table(MUTS$myeloid_cosmic_count))
  message("COSMIC done\n")
  ##
  
  
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
  
  ## about ~68 genes
  # ch.my.pd
  ch.my.pd.genes <- unique(A$Gene[A$ch_my_pd==1])
  
  A$keys = apply(A, 1, function(x) {
    names(x)[x != '*'] %>% .[!. %in% c('source', 'ch_my_pd', 'ch_pancan_pd')]
  })
  pd_dfs = split(A, paste(A$keys))
  
  # MUTS <- MUTS2
  
  MUTS$aa_pos <- as.character(MUTS$aa_pos)
  MUTS$Exon <- sapply(MUTS$EXON_VEP, function(x) str_split(x,"/")[[1]][1])
  colnames(MUTS)[colnames(MUTS)=="SYMBOL_VEP"] <- c("Gene")
  MUTS.temp <- MUTS
  matched = list()
  ch_my_variants = NULL
  
  message("annotating variants with PD Table...")
  
  ## anything truncating in PD_table
  truncating <- c("frameshift_variant", "splice_acceptor_variant","splice_donor_variant", "stop_gained")
  for (i in 1:length(pd_dfs)){
    print(i)
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
    pds = left_join(MUTS.temp, temp[,c(keys,'ch_my_pd','ch_pancan_pd','source')], by=keys)
    pds_ch_my_pd = pds %>% filter(ch_my_pd > 0)
    pds_ch_pd = pds %>% filter(ch_pancan_pd >=1) # & VariantClass %in% truncating
    # matched = append(matched, pds$var_key)
    matched = append(matched, pds_ch_my_pd$var_key)
    matched = append(matched, pds_ch_pd$var_key)
    ch_my_variants = rbind(ch_my_variants, pds_ch_my_pd, pds_ch_pd)
  }
  
  # just so we have the gene for reviewing
  ch_my_variants = ch_my_variants %>%
    dplyr::rename(gene = Gene) %>%
    dplyr::select(source, var_key, ch_my_pd, ch_pancan_pd, gene) %>% unique()
  
  tryCatch(
    expr = {
      ch_my_variants <- aggregate(source ~ var_key + ch_my_pd + ch_pancan_pd, data = ch_my_variants, FUN = paste, collapse = ",")
    },
    error = function(e){
      message(e)
    }
  )
  message("PD Table done\n")
  
  ## easier to just not return anything from TC and call this below it
  MUTS = left_join(MUTS, ch_my_variants, by="var_key")
  ## for the ifelse statement below, ch_my_pd CANNOT be NA!!!!!!!!!!!
  MUTS$ch_my_pd <- fillna(MUTS$ch_my_pd, 0)
  ## update ch_my_pd for those ~71 genes listed in Slack if B/B hotspot (12/22/2021 9:41 PM)
  MUTS$n.HGVSp <- fillna(MUTS$n.HGVSp, 0)
  MUTS$n.HGVSc <- fillna(MUTS$n.HGVSc, 0)
  MUTS$ch_my_pd <- ifelse((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes, 1, MUTS$ch_my_pd)
  ## initial ch_pd just get pancan and truncating, myeloid stuff comes below
  MUTS$ch_pd <- ifelse(MUTS$VariantClass %in% truncating & MUTS$ch_pancan_pd>=1, 1, 0)
  
  MUTS <- MUTS %>%
    mutate(WHY_CH_ch_my_pd = ifelse(MUTS$ch_my_pd>0 | MUTS$ch_pd>0, "PD_table;", ""))
  
  
  # 6. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, paste0(WHY_CH_ch_my_pd, " Cosmic_heme;"), WHY_CH_ch_my_pd))
  
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", 1, ch_pd),
           WHY_CH_ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", "OncoKB;", ""))
  
  MUTS$isOncogenic <- MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic" | MUTS$oncoKB=="Predicted Oncogenic"
  MUTS$isTSG <- MUTS$Gene %in% TSG
  
  # 2. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  # Genes not listed in the cancer census or OncoKB were reviewed in the literature to determine if they were potentially tumor suppressor genes.# Annotate PD based on prevelance in cancer databases (COSMIC, Oncokb)
  
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(isTSG & VariantClass %in% truncating, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(isTSG & VariantClass %in% truncating,
                                 paste0(WHY_CH_ch_pd, " TSG & Truncating;"), WHY_CH_ch_pd))
  
  ## is gene a truncating hotpot gene and is the variant in truncating class
  # vars.truncating <- vars %>%
  #   filter(truncating=="truncating") %>%
  #   group_by(Gene) %>%
  #   mutate(n=n()) %>%
  #   filter(n >=2)
  vars.truncating <- read.table(args$truncating, sep = "\t", header = T)
  
  truncating_genes <- unique(vars.truncating[vars.truncating$n>=10,"SYMBOL_VEP"])
  MUTS$isTruncatingHotSpot <- ifelse(MUTS$Gene %in% truncating_genes & MUTS$VariantClass %in% truncating, 1, 0)
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(isTruncatingHotSpot, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(isTruncatingHotSpot, paste0(WHY_CH_ch_pd, " BB truncating-HS;"), WHY_CH_ch_pd),
           ch_my_pd = ifelse(isTruncatingHotSpot & Gene %in% ch.my.pd.genes, 1, ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(isTruncatingHotSpot & Gene %in% ch.my.pd.genes, paste0(WHY_CH_ch_my_pd, " BB truncating-HS;"), WHY_CH_ch_my_pd))
  
  #3. Any variant reported as somatic at least 20 times in COSMIC
  MUTS$CosmicCount <- as.numeric(MUTS$CosmicCount)
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(CosmicCount>=20, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(CosmicCount>=20, paste0(WHY_CH_ch_pd, " Cosmic;"), WHY_CH_ch_pd))
  
  
  #4. Any variant meeting criteria for CH-Myeloid-PD as above.
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, 1, ch_pd),
           WHY_CH_ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, paste0(WHY_CH_ch_pd, " ch_my_pd>=1"), WHY_CH_ch_pd))
  
  MUTS$ch_pd <- fillna(MUTS$ch_pd, 0)
  annotate_PD.R = MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  
  # sample.x=grep("SAMPLE",colnames(x))
  # sample.annorate=grep("SAMPLE",colnames(annotate_PD.R))  # sample.annorate=grep("SAMPLE",colnames(annotate_PD.R))
  
  if(all(annotate_PD.R[,c("CHROM","POS","REF","ALT","SAMPLE")]==x[,c("CHROM","POS","REF","ALT","SAMPLE")])) {
    message("good same dim as orig")
    annotate_PD.R <- annotate_PD.R %>% # annotate_PD.R is the MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
      dplyr::mutate(gene_loci_p = paste(Gene, paste0(aa_ref, aa_pos), sep = "_"),
                    cDNAchange = gsub(".*:","", HGVSc_VEP),
                    gene_loci_c = paste(Gene, cDNAchange, sep = "_"),
                    gene_aachange = paste(Gene, AAchange, sep = "_"),
                    gene_loci_vep = ifelse(is.na(AAchange) | AAchange == "", gene_loci_c, gene_loci_p),
                    ch_pd2 = case_when(ch_pd==1 ~ 1,
                                       (gene_loci_vep %in% vars$gene_loci_vep |
                                          gene_aachange %in% matches.2.c.p$exact_match |
                                          gene_loci_c %in% matches.2.c.p$exact_match) ~ 1,
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
    return(
      annotate_PD.R %>% #
        dplyr::select("CHROM","POS","REF","ALT","SAMPLE",
                      #"n.HGVSp","n.HGVSc", ## ugh .x and .y oh well
                      "COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count",
                      "oncoKB","isOncogenic","isTSG", "isTruncatingHotSpot",
                      "ch_my_pd","ch_pd","ch_pd2",
                      "VariantClass","AAchange","Gene","WHY_CH")
    )
  } else {
    message("Annotate PD failed: something wrong with dims")
  }
}
# annotation <- annotate.PD(final.passed)


final.passed$Mutect2_PON_2AT2_percent <- fillna(final.passed$Mutect2_PON_2AT2_percent, 0)
final.passed$passed <- final.passed$Mutect2_FILTER=="PASS"

## let's look at 
    # gnomAD_MAX.Stringent.0007, gnomAD_MAX.lessStringent.0007, 
    # MAX_gnomAD_AF_VEP, MAX_gnomADe_AF_VEP, MAX_gnomADg_AF_VEP
# final.passed$Vardict_PON_2AT2_percent = 1 = fail
final.passed$passed_everything <- (!final.passed$Mutect2_PON_2AT2_percent &
                                     final.passed$alt_strand_counts_min_1_caller_only &
                                     final.passed$max.under.0.35 &
                                     final.passed$max.over.0.02 &
                                     final.passed$passed &
                                     final.passed$complexity_filters)
                                 
final.passed$passed_everything_mutect <- final.passed$passed_everything & final.passed$Mutect2_PASS
final.passed <- 
annotation <- annotate.PD(final.passed)
final.df <-  left_join(final.passed,
                       annotation,
                       by=c("CHROM","POS","REF","ALT","SAMPLE"))

final.df$sample_id <- args$sample_id
colnames(final.df)[(which(sapply(final.df, class)=="list"))]
write.table(final.df, paste0(args$out_folder,"/",args$sample_id,".final.tsv"), row.names = F, sep = "\t", quote = F)


