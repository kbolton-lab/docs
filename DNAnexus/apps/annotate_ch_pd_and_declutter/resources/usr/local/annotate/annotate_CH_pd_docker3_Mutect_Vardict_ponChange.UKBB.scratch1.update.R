#!/usr/local/bin/Rscript

library(argparser)

library(tidyverse)
library(stringr)
library(ggsci)
library(data.table)
library(httr)
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
parser <- add_argument(parser, "--input", type="character", help="input")
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

input <- read.table(args$input, comment.char = "", sep = "\t", quote="", header = T)
colnames(input)[colnames(input)=="AAchange.x"] <- "AAchange"

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
                  "oncoKB","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","COSMIC_ID",
                  "VariantClass","SAMPLE","EXON_VEP","AAchange")]
  
  # MUTS <- MUTS %>%
  #   mutate(VariantClass = str_split(Consequence_VEP,"&",simplify = TRUE)[,1]) %>%
  #   mutate(VariantClass = case_when((VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) > 0 ~ "inframe_deletion",
  #                                   (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 == 0 & nchar(REF)-nchar(ALT) < 0 ~ "inframe_insertion",
  #                                   (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF)-nchar(ALT) %% 3 != 0 ~ "frameshift_variant",
  #                                   (VariantClass=="coding_sequence_variant" | VariantClass=="protein_altering_variant") & nchar(REF) == nchar(ALT) ~ "missense_variant",
  #                                   TRUE ~ VariantClass)) %>%
  #   mutate(VariantClass = translate_consequence[VariantClass])
  
  MUTS$aa_ref <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1])
  MUTS$aa_alt <- sapply(MUTS$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][2])
  MUTS$aa_pos <- as.numeric(str_extract(MUTS$AAchange, "\\d+"))
  MUTS$var_key = paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  
  ###create list of tumor suppressor genes####
  ## file 1
  gene_census = read.table(args$TSG_file, comment.char = "", sep = "\t", quote = "", header = T)
  ## file 2
  oncoKB_curated = read.table(args$oncoKB_curated, comment.char = "", sep = "\t", quote = "", header = T)
  oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
  TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
  
  # message("annotating variants with oncoKB...")
  # json_query <- apply(MUTS, 1, function(x) {
  #   my_list <- list(
  #     alteration=x["AAchange"],
  #     consequence=x["VariantClass"],
  #     type="regular",
  #     gene=list(
  #       hugoSymbol=x["SYMBOL_VEP"]
  #     )
  #   )
  #   return(my_list)
  # })
  # 
  # POST_file <- paste0("/scratch1/fs1/bolton/brian/",args$sample_id,".oncoKB.json")
  # RETURN_file <- paste0("/scratch1/fs1/bolton/brian/",args$sample_id,".oncoKB.return.tsv")
  # write(jsonlite::toJSON(unname(json_query), auto_unbox = T), POST_file)
  # POST_command <- paste0('curl -X POST -H "Content-Type: application/json" -d@',
  #                        POST_file, 
  #                        ' https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange -H "Authorization: Bearer da03f096-284a-490a-8b3f-95c5f1bf5666" | jq ".[].oncogenic" > ',
  #                        RETURN_file)
  # system(command = POST_command, intern = T)
  # MUTS$oncoKB <- read.table(RETURN_file, quote='"', sep = "\t", header = F, blank.lines.skip=FALSE)$V1
  # system(command = paste0("rm ", POST_file), intern = T)
  # system(command = paste0("rm ", RETURN_file), intern = T)
  # 
  # message("oncoKB done\n")
  
  ###########annotate with COSMIC########
  ## MUTS MUST BE ORDERED by CHROM (and maybe POS? just to be safe?? Basically whatever column is paralleled needs to be sorted by...)
  # message("annotating variants with COSMIC...")
  # MUTS <- cosmic.run(MUTS, "SAMPLE", "var_key")
  # MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  # MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  # MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  # message("COSMIC done\n")
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
  truncating <- c("frameshift_variant", "splice_region_variant", "stop_gained")
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



annotation <- annotate.PD(input)
my_pd <- all(annotation$ch_my_pd >= input$ch_my_pd)
pd <- all(annotation$ch_pd >= input$ch_pd)
pd2 <- all(annotation$ch_pd2 >= input$ch_pd2)
print(paste0("####################################################### extra ch_my_pd:     ", sum(annotation$ch_my_pd > input$ch_my_pd)))
print(paste0("####################################################### extra ch_pd:     ", sum(annotation$ch_pd > input$ch_pd)))
print(paste0("####################################################### extra ch_pd2:     ", sum(annotation$ch_pd2 > input$ch_pd2)))

if (all(input[,1:4] == annotation[,1:4]) & my_pd & pd & pd2) {
  input[,c("ch_my_pd","ch_pd","ch_pd2")] <- annotation[,c("ch_my_pd","ch_pd","ch_pd2")]
  
  write.table(input, paste0(args$out_folder,"/",args$sample_id,".final.updated.tsv"), row.names = F, sep = "\t", quote = F)
}




