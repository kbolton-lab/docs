#!/usr/local/bin/Rscript

library(argparser)
library(tidyverse)
library(stringr)
library(ggsci)
library(data.table)
library(httr)
library(sqldf)


'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

#'@returns dataframe with equal number of rows as input with binded counts columns
#'@param df dataframe: this is the MUTS dataframe; needs to have `HGVSp_VEP` column
#'@param sample.column string: column name for sample name/id
cosmic.run <- function(df, sample.column) {
  
  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  X.1 <- parallel::mclapply(ls(X1), function(x) {
    # X.1 <- lapply(ls(X1)[7:8], function(x) {
    df.chr <- X1[[x]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP)
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    cosmic <- read.table(paste0("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.final.more.minimal.test.",x,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
    colnames(cosmic) <- c("COSMIC_ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")
    
    cosmic.test <- sqldf("SELECT l.*, r.COSMIC_ID, r.var_key as var_key_cosmic, HGVSP, r.CosmicCount, 
              r.heme_cosmic_count, r.myeloid_cosmic_count
              FROM `df.chr` as l
              LEFT JOIN `cosmic` as r
              on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
    cosmic.test <- cosmic.test[,c(colnames(df.chr),"CosmicCount","heme_cosmic_count","myeloid_cosmic_count")]
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
parser <- add_argument(parser, "--input", type="character", help="in file")
parser <- add_argument(parser, "--output", type="character", help="out file")
parser <- add_argument(parser, "--TSG-file", type="character", help="tumor suppressor genes")
parser <- add_argument(parser, "--oncoKB-curated", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--pd-annotation-file", type="character", help="PD annoation files from papaemmanuil lab") # updated
parser <- add_argument(parser, "--pan-myeloid", type="character", help="panmyeloid_variant_counts.tsv") # updated
parser <- add_argument(parser, "--vars", type="character",help="help") # updated
parser <- add_argument(parser, "--mut1-bick", type="character", help="") # updated
parser <- add_argument(parser, "--mut1-kelly", type="character",help="") # updated


args <- parse_args(parser)
args <- parse_args(parser, list(c("-i","/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/12/update/transplant/4006226_23153_0_0.final.tsv"),
                                c("-o","/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/4006226_23153_0_0.final.reannotate.tsv"),
                                c("-T","~/Bolton/data/gene_census_TSG.txt"),
                                c("--oncoKB-curated","~/Bolton/data/all_curated_genes_v2.0.tsv"),
                                c("-p","~/Bolton/data/pd_table_kbreview_bick_trunc2.txt"),
                                c("--pan-myeloid","~/Bolton/data/panmyeloid_variant_counts.vep.annotated.vcf.tsv"),
                                c("-v","~/Bolton/data/bick.bolton.vars2.txt"),
                                c("--mut1-bick","~/Bolton/data/topmed.mutation.1.c.p.txt"),
                                c("--mut1-kelly","~/Bolton/data/kelly.mutation.1.c.p.txt")))
#####
## REDO Annotate_PD
# source("hotspots_final.R")
getVAFs <- function(x) {
  mutect_VAF <- grep("Mutect2_gt_AF", colnames(x))
  vardict_VAF <- grep("Vardict_gt_AF", colnames(x))
  VAFS <- c(mutect_VAF, vardict_VAF)
  VAFS
}
AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               '%3D'='=', '='='=')


final.test <- read.table(args$input, sep = "\t", header = T, quote = "", comment.char = "")
final.test <- final.test[,c(1:249,263)]
final.test$max.under.0.35 <- apply(final.test[,getVAFs(final.test)],1,function(x) max(x, na.rm = T) < 0.35)
vardict.pon2 <- read.table("/Users/brian/apps/vardict_pon2at2percent_final/resources/usr/local/test/vardict.2N.tsv", 
                           sep = "\t", header = F, quote = "", comment.char = "")
colnames(vardict.pon2) <- c("CHROM","POS","REF","ALT")
vardict.pon2$Vardict_PON_2AT2_percent <- 1
mutect.pon2 <- read.table("/Users/brian/Bolton/CWL_TESTS/mutect2.merged.SB.2N.tsv", 
                           sep = "\t", header = F, quote = "", comment.char = "")
colnames(mutect.pon2) <- c("CHROM","POS","REF","ALT")
mutect.pon2$Mutect2_PON_2AT2_percent <- 1

final.test <- final.test %>%
  dplyr::select(-c(Mutect2_PON_2AT2_percent,Vardict_PON_2AT2_percent)) %>% 
  left_join(mutect.pon2 %>% 
              dplyr::select(CHROM, POS, REF, ALT, Mutect2_PON_2AT2_percent),
            by = c("CHROM"="CHROM", "POS"="POS", "REF"="REF", "ALT"="ALT")) %>%
  left_join(vardict.pon2 %>% 
              dplyr::select(CHROM, POS, REF, ALT, Vardict_PON_2AT2_percent),
            by = c("CHROM"="CHROM", "POS"="POS", "REF"="REF", "ALT"="ALT"))

final.test$Mutect2_PON_2AT2_percent <- fillna(final.test$Mutect2_PON_2AT2_percent, 0)
final.test$Vardict_PON_2AT2_percent <- fillna(final.test$Vardict_PON_2AT2_percent, 0)

final.test$AAchange <- gsub("(.*p\\.)(.*)", "\\2", final.test$HGVSp)
for (i in 1:length(AminoAcids)) {
  final.test$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], final.test$AAchange)
}
final.test$gene_loci_p <- paste(final.test$SYMBOL_VEP,
                                paste0(sapply(final.test$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1]),
                                       as.numeric(str_extract(final.test$AAchange, "\\d+"))),
                                sep = "_")
final.test$gene_loci_c <- paste(final.test$SYMBOL_VEP,
                                gsub(".*:", "", final.test$HGVSc),
                                sep = "_")
final.test$gene_loci_vep <- ifelse(is.na(final.test$gene_loci_p),final.test$gene_loci_c,final.test$gene_loci_p)
final.test$key <- with(final.test, paste(CHROM,POS,REF,ALT,sep = ":"))
final.test$gene_aachange <- with(final.test, paste(SYMBOL_VEP, AAchange, sep = "_"))

vars <- read.table(args$vars,sep = "\t", header = T, comment.char = "")
vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
topmed.mutation.1 <- read.table(args$mut1_bick, sep = "\t", header = T, comment.char = "")
kelly.mutation.1 <- read.table(args$mut1_kelly, sep = "\t", header = T, comment.char = "")


dims <- dim(final.test)[[1]]
# final.test <- sqldf("SELECT l.*, r.`n.loci.vep`, r.`source.totals.loci`
#             FROM `final.test` as l
#             LEFT JOIN `vars` as r
#             on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
final.test <- sqldf("SELECT l.*, r.`n.loci.vep`
            FROM `final.test` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
final.test <- final.test[!duplicated(final.test),]
## make sure aachange exists as in doesn't end with an '_'; example: DNMT3A_ for splice 
# final.test <- sqldf("SELECT l.*, r.`n.HGVSp`, r.`source.totals.p`
#             FROM `final.test` as l
#             LEFT JOIN `vars` as r
#             on l.key = r.key OR (l.gene_aachange = r.gene_aachange) AND r.gene_aachange NOT LIKE '%_'")
# final.test <- final.test[!duplicated(final.test),]
# final.test <- sqldf("SELECT l.*, r.`n.HGVSc`, r.`source.totals.c`
#             FROM `final.test` as l
#             LEFT JOIN `vars` as r
#             on l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange) AND r.gene_cDNAchange NOT LIKE '%_'")
# final.test <- final.test[!duplicated(final.test),]
print(paste0("dims match: ", dim(final.test)[[1]] == dims))

################################################
### Annotate PD

AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               'Ter'= 'X', '='='=', 'del'='DEL')

annotate.PD <- function(x) {
  x <- final.test
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
    "splice_accepter_variant" = "Splice_Site",
    "3_prime_UTR_variant" = "3_prime_UTR_variant",
    "5_prime_UTR_variant" = "5'Flank",
    "stop_lost" = "Nonsense_Mutation",
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
  source("~/Bolton/data/trying.parallel.cosmic.test.R")
  ## need to add a parallel
  MUTS.t <- cosmic.run(MUTS,"SAMPLE")
  MUTS <- MUTS.t
  MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  MUTS$heme_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  
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
  unique(A$keys)
  
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
  MUTS$ch_my_pd = fillna(MUTS$ch_my_pd, 0)
  
  
  # 6. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times
  
  MUTS$ch_my_pd = ifelse(MUTS$heme_cosmic_count>=10,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  
  # 7. Any variant noted as potentially oncogenic in an in-house dataset of 7,000 individuals with myeloid neoplasm greater than or equal to 5 times
  panmyeloid_variant_counts = read_delim(args$pan_myeloid, delim = "\t", guess_max = 5e6)
  panmyeloid_variant_counts$var_key = with(panmyeloid_variant_counts, paste(CHROM,POS,REF,ALT, sep = ":"))
  panmyeloid_variant_counts <- panmyeloid_variant_counts[panmyeloid_variant_counts$Annotation == "ONCOGENIC" | 
                                                           panmyeloid_variant_counts$Annotation == "SNP",]
  # panmyeloid_variant_counts$MDS <- fillna(panmyeloid_variant_counts$MDS, 0)
  # panmyeloid_variant_counts$AML <- fillna(panmyeloid_variant_counts$AML, 0)
  # panmyeloid_variant_counts$MPN <- fillna(panmyeloid_variant_counts$MPN, 0)
  # MUTS1 = left_join(MUTS, panmyeloid_variant_counts %>% dplyr::select(Annotation,MDS,AML,MPN,var_key), by = 'var_key') %>%
  #   mutate(
  #     MDS = fillna(MDS, 0),
  #     AML = fillna(AML, 0),
  #     MPN = fillna(MPN, 0),
  #     Annotation = fillna(Annotation, '')
  #   ) %>%
  #   mutate(n_panmyeloid = ifelse(Annotation == 'ONCOGENIC', MDS + AML + MPN, 0))
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
  
  
  MUTS$ch_my_pd = ifelse(MUTS$n_panmyeloid>=5,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  
  MUTS$ch_pd = ifelse(MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic", 1, 0)
  MUTS$isOncogenic <- MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic"
  MUTS$isTSG <- MUTS$Gene %in% TSG
  
  # 2. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  # Genes not listed in the cancer census or OncoKB were reviewed in the literature to determine if they were potentially tumor suppressor genes.# Annotate PD based on prevelance in cancer databases (COSMIC, Oncokb)
  
  
  MUTS$ch_pd = ifelse(MUTS$isTSG &
                        MUTS$VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"), 1, MUTS$ch_pd)
  
  #3. Any variant reported as somatic at least 20 times in COSMIC
  MUTS$CosmicCount <- as.numeric(MUTS$CosmicCount)
  MUTS$ch_pd = ifelse(MUTS$CosmicCount>=20, 1, MUTS$ch_pd)
  
  #4. Any variant meeting criteria for CH-Myeloid-PD as above.
  
  MUTS = MUTS %>% mutate(ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, 1, ch_pd))
  
  
  annotate_PD.R = MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  
  sample=grep("SAMPLE",colnames(x))
  
  if(all(annotate_PD.R[,c(1:4,9)]==x[,c(1:4,sample)])) {
    message("good same dim as orig")
    return(
      annotate_PD.R %>% # annotate_PD.R is the MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
        dplyr::mutate(gene_loci_p = paste(Gene, paste0(aa_ref, aa_pos), sep = "_"),
                      cDNAchange = gsub(".*:","", HGVSc_VEP),
                      gene_loci_c = paste(Gene, cDNAchange, sep = "_"),
                      gene_aachange = paste(Gene, AAchange, sep = "_"),
                      gene_loci_vep = ifelse(is.na(AAchange) | AAchange == "", gene_loci_c, gene_loci_p),
                      ch_pd2 = case_when(ch_pd==1 ~ 1,
                                         (gene_loci_vep %in% vars$gene_loci_vep | 
                                            gene_aachange %in% topmed.mutation.1$exact_match | gene_loci_c %in% topmed.mutation.1$exact_match |
                                            gene_aachange %in% kelly.mutation.1$exact_match | gene_loci_c %in% kelly.mutation.1$exact_match
                                         ) & VariantClass=="Missense_Mutation" ~ 1,
                                         gene_loci_vep %in% vars$gene_loci_vep & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
                                         TRUE ~ 0))
    )
  } else {
    message("something wrong with dims")
  }
  
}

bf.correction <- 0.05/38793002
final.test$passed_everything <- (!final.test$Vardict_PON_2AT2_percent &
                                     !final.test$Mutect2_PON_2AT2_percent &
                                     final.test$alt_strand_counts_min_2_callers &
                                     final.test$max.under.0.35 &
                                     final.test$max.over.0.02 &
                                     final.test$passed &
                                     final.test$complexity_filters &
                                     final.test$PON_FISHER <= bf.correction)
                                 
final.test$passed_everything_mutect <- final.test$passed_everything & final.test$Mutect2_PASS

annotation <- annotate.PD(final.test)
final.df <-  left_join(final.test,
                       annotation %>%
                         dplyr::select("CHROM","POS","REF","ALT","SAMPLE","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","MDS",
                                       "AML","MPN","ch_my_pd","ch_pd","ch_pd2","VariantClass","Gene"),
                       by=c("CHROM","POS","REF","ALT","SAMPLE"))
final.df$AAchange2 <- final.df$AAchange
write.table(final.df, args$output, row.names = F, sep = "\t", quote = F)

