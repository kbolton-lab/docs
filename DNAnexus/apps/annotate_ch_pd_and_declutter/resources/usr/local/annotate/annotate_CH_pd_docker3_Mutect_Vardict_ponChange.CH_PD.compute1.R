#!/usr/local/bin/Rscript

library(argparser)
library(tidyverse)
library(stringr)
# library(readxl)
library(data.table)
library(httr)
library(sqldf)
library(jsonlite)

# source("~/Bolton/UKBB/args.R")
'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}




parser <- arg_parser("Process NGS variant calls")
parser <- add_argument(parser, "--input", type="character", help="in")
parser <- add_argument(parser, "--output", type="character", help="out")

parser <- add_argument(parser, "--TSG-file", type="character", help="tumor suppressor genes")
parser <- add_argument(parser, "--oncoKB-curated", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--pd-annotation-file", type="character", help="PD annoation files from papaemmanuil lab") # updated
parser <- add_argument(parser, "--pan-myeloid", type="character", help="panmyeloid_variant_counts.tsv") # updated
parser <- add_argument(parser, "--bolton-bick-vars", type="character", help="help") # updated
parser <- add_argument(parser, "--mut2-bick", type="character", help="") # updated
parser <- add_argument(parser, "--mut2-kelly", type="character", help="") # updated
parser <- add_argument(parser, "--matches2", type="character", help="") # updated MORE


args <- parse_args(parser)
# args$input <- "/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/lung/update/pass/pon_pass_myeloid/1007589_23153_0_0.pon.pass.myeloid.tsv"
# args$input <- "/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/lung/update/pass/pon_pass_myeloid/update_ch_pd/1574340_23153_0_0.pon.pass.myeloid.ch_pd.tsv"
# args$input2 <- "/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/lung/update/pass/pon_pass_myeloid/1574340_23153_0_0.pon.pass.myeloid.tsv"
input <- read.table(args$input, header = T, quote = "", sep = "\t", comment.char = "")
input2 <- read.table(args$input2, header = T, quote = "", sep = "\t", comment.char = "")
input <- input[with(input, order(CHROM, POS, ALT, REF)),]
input_clean <- input
input_clean <- input_clean[,-which(colnames(input_clean)=="AAchange.y")]

colnames(input_clean)[colnames(input_clean) == "AAchange.x"] <- "AAchange"
which(input_clean$ch_pd2==1)
# test.clean <- input_clean[which(input_clean$ch_pd2==1),]
 ################################################
### Annotate PD

AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               '%3D'='=', '='='=')

input_clean$AAchange <- gsub("(.*p\\.)(.*)", "\\2", input_clean$HGVSp)
for (i in 1:length(AminoAcids)) {
  input_clean$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], input_clean$AAchange)
}

input_clean$gene_loci_p <- paste(input_clean$SYMBOL_VEP,
                                  paste0(sapply(input_clean$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1]),
                                         as.numeric(str_extract(input_clean$AAchange, "\\d+"))),
                                  sep = "_")
input_clean$gene_loci_c <- paste(input_clean$SYMBOL_VEP,
                                  gsub(".*:", "", input_clean$HGVSc),
                                  sep = "_")
input_clean$gene_loci_vep <- ifelse(is.na(input_clean$gene_loci_p),input_clean$gene_loci_c,input_clean$gene_loci_p)
input_clean$key <- with(input_clean, paste(CHROM,POS,REF,ALT,sep = ":"))
input_clean$gene_aachange <- with(input_clean, paste(SYMBOL_VEP, AAchange, sep = "_"))
input_clean$gene_cDNAchange <- paste(input_clean$SYMBOL_VEP, gsub(".*:","",input_clean$HGVSc_VEP), sep="_")

## ch_pd stuff
vars <- read.table(args$bolton_bick_vars, sep = "\t", header = T, comment.char = "")
vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")
topmed.mutation.2 <- read.table(args$mut2_bick, sep = "\t", header = T, comment.char = "")
kelly.mutation.2 <- read.table(args$mut2_kelly, sep = "\t", header = T, comment.char = "")
##both Bick/Bolton
matches.2.c.p <- read.table(args$matches2, sep = "\t", header = T, comment.char = "")


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
    
    # "coding_sequence_variant" = "feature_truncation", # manually assigned based on length ref & alt
    # "protein_altering_variant" = "feature_truncation", # manually assigned based on length ref & alt
    "feature_truncation" = "feature_truncation",
    
    "3_prime_UTR_variant" = "3_prime_UTR_variant", # missing
    "5_prime_UTR_variant" = "5_prime_UTR_variant", # missing
    
    
    "stop_gained" = "stop_gained", #8
    "start_lost" = "start_lost", #9
    "stop_lost" = "missense_variant", # A sequence variant where at least one base of the terminator codon (stop) is changed, resulting in an elongated transcript
    "stop_retained_variant" = "synonymous_variant" # A sequence variant where at least one base in the terminator codon is changed, but the terminator remains
  )
  
  MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP","HGVSp_VEP","n.HGVSp","HGVSc_VEP","n.HGVSc",
                  "Consequence_VEP","SAMPLE","EXON_VEP","AAchange", 
                  "CosmicCount", "heme_cosmic_count", "myeloid_cosmic_count", "oncoKB", "isTSG")]
  # translate_type = c(
  #   'SNV' = 'SNV',
  #   'deletion' = 'Del',
  #   'insertion' = 'Ins')
  
  MUTS <- MUTS %>%
    mutate(VariantClass = case_when((Consequence_VEP=="coding_sequence_variant" | Consequence_VEP=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) > 0 ~ "inframe_deletion",
                                    (Consequence_VEP=="coding_sequence_variant" | Consequence_VEP=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) < 0 ~ "inframe_insertion",
                                    (Consequence_VEP=="coding_sequence_variant" | Consequence_VEP=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 != 0 ~ "frameshift_variant",
                                    (Consequence_VEP=="coding_sequence_variant" | Consequence_VEP=="protein_altering_variant") & nchar(REF) == nchar(ALT) ~ "missense_variant",
                                    TRUE ~ Consequence_VEP)) %>%
    mutate(VariantClass = str_split(Consequence_VEP,"&",simplify = TRUE)[,1]) %>%
    mutate(VariantClass = translate_consequence[VariantClass])

  
  MUTS$aa_ref <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1])
  MUTS$aa_alt <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][2])
  MUTS$aa_pos <- as.numeric(str_extract(MUTS$AAchange, "\\d+"))
  print("##########################################")
  # print(MUTS$aa_pos)
  print("##########################################")
  MUTS$var_key = paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  
  ###create list of tumor suppressor genes####
  ## file 1
  gene_census = read.table(args$TSG_file, comment.char = "", sep = "\t", quote = "", header = T)
  ## file 2
  oncoKB_curated = read.table(args$oncoKB_curated, comment.char = "", sep = "\t", quote = "", header = T)
  oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
  
  TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
  
  
  ###########annotate with COSMIC########
  ## need to add a parallel
  # message("annotating variants with COSMIC...")
  # MUTS <- cosmic.run(MUTS, "SAMPLE")
  # MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  # MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  # MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  # 
  # MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  # MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  # message("COSMIC done\n")
    
  
  
  #annotate oncokb
  # get_oncokb = function(mut) {
  # 
  #   # apiKey = "a83627cd-47e4-4be0-82dd-8f4cc9e4d6d0"
  #   apiKey = "da03f096-284a-490a-8b3f-95c5f1bf5666"
  #   request_url = paste0("https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
  #                        mut['SYMBOL_VEP'], "&alteration=", mut['AAchange'], "&consequence=", mut['VariantClass'])
  # 
  # 
  # 
  # 
  #   res=httr::content(httr::GET(request_url, httr::add_headers(Authorization = paste("Bearer", apiKey))))
  #   return(res$oncogenic)
  # }
  # 
  # message("annotating variants with oncoKB...")
  # h = curl::new_handle()
  # curl::handle_setopt(h, http_version = 2)
  # httr::set_config(httr::config(http_version = 0))
  # cl = parallel::makeCluster(4)
  # MUTS$oncoKB = parallel::parApply(cl, MUTS, 1, get_oncokb)
  # parallel::stopCluster(cl)
  # message("oncoKB done\n")
  # MUTS2 <- MUTS
  
  
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
  ## about ~71 genes
  # ch.my.pd
  ch.my.pd.genes <- unique(A$Gene[A$ch_my_pd==1])

  ch_my_variants = NULL
  pd_dfs = split(A, paste(A$keys))
  MUTS$aa_pos <- as.character(MUTS$aa_pos)
  MUTS$Exon <- sapply(MUTS$EXON_VEP, function(x) str_split(x,"/")[[1]][1])
  colnames(MUTS)[colnames(MUTS)=="SYMBOL_VEP"] <- c("Gene")
  colnames(MUTS)[colnames(MUTS)=="AAchange.y"] <- c("AAchange")
  MUTS.temp <- MUTS
  matched = list()

  message("annotating variants with PD Table...")
  for (i in 1:length(pd_dfs)) {
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
  message("PD Table done\n")
  
  ## easier to just not return anything from TC and call this below it
  MUTS = left_join(MUTS, ch_my_variants, by="var_key")
  ## for the ifelse statement below, ch_my_pd CANNOT be NA!!!!!!!!!!!
  MUTS$ch_my_pd <- fillna(MUTS$ch_my_pd, 0)
  ## update ch_my_pd for those ~71 genes listed in Slack if B/B hotspot (12/22/2021 9:41 PM)
  MUTS$n.HGVSp <- fillna(MUTS$n.HGVSp, 0) # mother F-word!! if MUTS$Gene %in% ch.my.pd.genes need this if MUTS$n.HGVSp is NA, example 1415947_23153_0_0.final.tsv
  MUTS$n.HGVSc <- fillna(MUTS$n.HGVSc, 0) # mother F-word!! if MUTS$Gene %in% ch.my.pd.genes need this if MUTS$n.HGVSp is NA, example 1415947_23153_0_0.final.tsv
  MUTS$ch_my_pd <- ifelse((MUTS$n.HGVSp>=2 | MUTS$n.HGVSc>=2) & MUTS$Gene %in% ch.my.pd.genes, 1, MUTS$ch_my_pd)  ## might need to rerun pd_table again!!!
  MUTS <- MUTS %>% 
    mutate(WHY_CH_ch_my_pd = ifelse(MUTS$ch_my_pd>0, "PD_table;", ""))
  
  
  # 6. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times
  
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, paste0(WHY_CH_ch_my_pd, " Cosmic_heme;"), WHY_CH_ch_my_pd))
  
  # 7. Any variant noted as potentially oncogenic in an in-house dataset of 7,000 individuals with myeloid neoplasm greater than or equal to 5 times
  panmyeloid_variant_counts = read.table(args$pan_myeloid, comment.char = "", sep = "\t", quote = "", header = T)
  panmyeloid_variant_counts$var_key = with(panmyeloid_variant_counts, paste(CHROM,POS,REF,ALT, sep = ":"))
  panmyeloid_variant_counts <- panmyeloid_variant_counts[panmyeloid_variant_counts$Annotation == "ONCOGENIC" | 
                                                           panmyeloid_variant_counts$Annotation == "SNP",]
  
  MUTS$Gene_HGVSp_VEP <- paste(MUTS$Gene, gsub(".*p.","", MUTS$HGVSp_VEP), sep = "_")
  panmyeloid_variant_counts$Gene_HGVSp_VEP <- paste(panmyeloid_variant_counts$SYMBOL_VEP, gsub(".*p.","", panmyeloid_variant_counts$HGVSp_VEP), sep = "_")
  
  message("annotating variants with panmyeloid...")
  # MUTS <- sqldf("SELECT l.*, r.MDS, r.AML, r.MPN, r.Annotation
  #           FROM `MUTS` as l
  #           LEFT JOIN `panmyeloid_variant_counts` as r
  #           on l.var_key = r.var_key OR l.Gene_HGVSp_VEP = r.Gene_HGVSp_VEP")
  message("Please don't crash here...")
  
  message("Using No SQL, Manual Binding")
  MUTS_pan_HGVSp <- left_join(MUTS, panmyeloid_variant_counts %>% dplyr::select(Gene_HGVSp_VEP, MDS, AML, MPN, Annotation), by = c("Gene_HGVSp_VEP"="Gene_HGVSp_VEP"))
  MUTS_pan_var_key <- MUTS_pan_HGVSp %>% filter(is.na(Annotation))
  MUTS_pan_HGVSp <- MUTS_pan_HGVSp %>% filter(!is.na(Annotation))
  MUTS_pan_var_key <- left_join(MUTS_pan_var_key %>% dplyr::select(-MDS, -AML, -MPN, -Annotation), panmyeloid_variant_counts %>% dplyr::select(var_key, MDS, AML, MPN, Annotation), by = c("var_key"="var_key"))
  MUTS <- rbind(MUTS_pan_var_key, MUTS_pan_HGVSp)
  MUTS <- MUTS %>% arrange(CHROM, POS, ALT, REF)
  
  MUTS <- MUTS %>%
    mutate(
      MDS = fillna(MDS, 0),
      AML = fillna(AML, 0),
      MPN = fillna(MPN, 0),
      Annotation = fillna(Annotation, ''),
      n_panmyeloid = ifelse(Annotation == 'ONCOGENIC', MDS + AML + MPN, 0))
  message("panmyeloid done\n")
  
  
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(n_panmyeloid>=5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(MUTS$n_panmyeloid>=5, paste0(WHY_CH_ch_my_pd, " Pan_Myeloid;"), WHY_CH_ch_my_pd))
  
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", 1, 0),
           isOncogenic = oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic",
           WHY_CH_ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", "OncoKB;", ""))
  MUTS$isTSG <- MUTS$Gene %in% TSG


  # 2. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  # Genes not listed in the cancer census or OncoKB were reviewed in the literature to determine if they were potentially tumor suppressor genes.# Annotate PD based on prevelance in cancer databases (COSMIC, Oncokb)
  # "The nonsense is known as stop_gained" https://www.biostars.org/p/202928/
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(isTSG & VariantClass %in% c("frameshift_variant", "splice_region_variant", "stop_gained"), 1, ch_pd),
           WHY_CH_ch_pd = ifelse(isTSG & VariantClass %in% c("frameshift_variant", "splice_region_variant", "stop_gained"), 
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

  MUTS$ch_pd = ifelse(MUTS$Gene %in% TSG &
                        MUTS$VariantClass %in% c("frameshift_variant", "splice_region_variant", "stop_gained"), 1, MUTS$ch_pd)

  
  annotate_PD.R = MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  sample=grep("SAMPLE",colnames(x))
  annotate.sample <- grep("SAMPLE",colnames(annotate_PD.R))
  
  if(all(annotate_PD.R[,c(1:4,annotate.sample)]==x[,c(1:4,sample)])) {
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
                      "CosmicCount","heme_cosmic_count","myeloid_cosmic_count",
                      "MDS", "AML","MPN","n_panmyeloid",
                      "oncoKB","isOncogenic","isTSG",
                      "ch_my_pd","ch_pd","ch_pd2",
                      "VariantClass","AAchange","Gene","WHY_CH"))
  } else {
    message("something wrong with dims")
  }

  # if(all(annotate_PD.R[,c(1:4,8)]==x[,c(1:4,sample)])) {
  #   message("good same dim as orig")
  #   return(
  #     annotate_PD.R %>% # annotate_PD.R is the MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  #       dplyr::mutate(loci = str_extract(AAchange, "[A-Z]\\d+"),
  #                     AAchange = str_remove(AAchange, "p."),
  #                     gene_aachange = paste(Gene,AAchange,sep = "_"),
  #                     gene_loci = paste(Gene,loci,sep = "_"),
  #                     ch_pd2 = case_when(ch_pd==1 ~ 1,
  #                                             gene_loci %in% topmed.loci.n$gene_loci & gene_aachange %in% topmed.mutation.1$gene_aachange & VariantClass=="Missense_Mutation" ~ 1,
  #                                             gene_loci %in% topmed.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
  #                                             gene_loci %in% kelly.loci.n$gene_loci & gene_aachange %in% kelly.mutation.1$gene_aachange & VariantClass=="Missense_Mutation" ~ 1,
  #                                             gene_loci %in% kelly.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
  #                                             TRUE ~ 0))
  # 
  #   )
  # } else {
  #   message("something wrong with dims")
  # }
}


annotation <- annotate.PD(input_clean)
if (all(input[,1:4] == annotation[,1:4])) {
  input[, c("ch_my_pd", "ch_pd", "ch_pd2")] = annotation[, c("ch_my_pd", "ch_pd", "ch_pd2")]
  # input$oncoKB <- annotation$oncoKB
  # input$isOncogenic <- annotation$isOncogenic
  # input$isTSG <- annotation$isTSG
  input$VariantClass <- annotation$VariantClass
  input$isOncogenic <- annotation$isOncogenic
  input$WHY_CH <- annotation$WHY_CH
  write.table(input, args$output, row.names = F, quote = F, sep = "\t")
  which(input$ch_pd2==1)
} else {
  print("Something went wrong")
}
