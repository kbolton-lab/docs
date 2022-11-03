#!/usr/bin/env Rscript
library(argparser)
library(sqldf)
library(jsonlite)


parser <- arg_parser("Process NGS variant calls")

parser <- add_argument(parser, "--TSG-file", type="character", help="tumor suppressor genes")
parser <- add_argument(parser, "--oncoKB-curated", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--pd-annotation-file", type="character", help="PD annoation files from papaemmanuil lab") # updated
parser <- add_argument(parser, "--pan-myeloid", type="character", help="panmyeloid_variant_counts.tsv") # updated
parser <- add_argument(parser, "--bolton-bick-vars", type="character", help="help") # updated
parser <- add_argument(parser, "--mut2-bick", type="character", help="") # updated
parser <- add_argument(parser, "--mut2-kelly", type="character", help="") # updated
parser <- add_argument(parser, "--matches2", type="character", help="") # updated
parser <- add_argument(parser, "--target-length", type="integer", help="total sum of target intervals", default=39000000)
parser <- add_argument(parser, "--input", type="character", help="input")
parser <- add_argument(parser, "--output", type="character", help="output")
parser <- add_argument(parser, "--col-order", type="character", help="--col-order")

args <- parse_args(parser)

# args <- parse_args(parser, list(c("-T","/Volumes/bolton/Active/projects/annotation_files/gene_census_TSG.txt"),
#                                 c("--oncoKB-curated","/Volumes/bolton/Active/projects/annotation_files/all_curated_genes_v2.0.tsv"),
#                                 c("-p","/Volumes/bolton/Active/projects/annotation_files/pd_table_kbreview_bick_trunc2.txt"), # new
#                                 c("--pan-myeloid","/Volumes/bolton/Active/projects/annotation_files/panmyeloid_variant_counts.vep.annotated.vcf.tsv"), # new
#                                 c("--bolton-bick-vars","/Volumes/bolton/Active/projects/annotation_files/bick.bolton.vars3.txt"), # new
#                                 c("--mut2-bick","/Volumes/bolton/Active/projects/annotation_files/topmed.n2.mutation.c.p.txt"), # new
#                                 c("--mut2-kelly","/Volumes/bolton/Active/projects/annotation_files/kelly.n2.mutation.c.p.txt"),
#                                 c("--matches2","/Volumes/bolton/Active/projects/annotation_files/matches.2.c.p.txt"), # MORE NEW!
#                                 c("--col-order","/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/1000051_23153_0_0.columns.txt"),
#                                 c("--target-length", 38997831),
#                                 c("--output","/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/already_ran/1200213_23153_0_0.updated.tsv"),
#                                 c("-i","/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/already_ran/1200213_23153_0_0.final.tsv")))
library(tidyverse)

'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

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

input <- read.table(args$input, sep = "\t", comment.char = "", quote = "", header = TRUE)
cols = read.table(args$col_order, header = F, sep = "\t", comment.char = "", quote="\"")
# cols = read.table("/Volumes/bolton/Active/projects/mocha/UKBB/exome/results/germline/1000051_23153_0_0.columns.txt",
#                   header = F, sep = "\t", comment.char = "", quote="\"")
output <- args$output
# PON
bf.corr <- .05/args$target_length
input$PON_FISHER <- fillna(input$PON_FISHER, 0)
input <- input[input$PON_FISHER <= bf.corr,]

for (col in c("ch_my_pd", "ch_pd2", "Gene", "AAchange")) {
  if (paste0(col,".x") %in% colnames(input)) {
    input <- input %>% dplyr::select(-paste0(col,".x"))
  }
  # if(paste0(col,".y") %in% colnames(test.ar2)){
  #   colnames(test.ar2)[colnames(test.ar2)==paste0(col,".y")] <- col
  # }
  colnames(input)[colnames(input)==paste0(col,".y")] <- col
}

cols.rem <- c("MDS", "AML","MPN","n_panmyeloid",
              "isOncogenic","isTSG",
              "ch_my_pd","ch_pd","ch_pd2",
              "VariantClass","AAchange","Gene","WHY_CH")
input <- input %>% dplyr::select(-cols.rem[cols.rem %in% colnames(input)])


#########################################################################################################################
#########################################################################################################################



Vardict.complex <- c("Vardict_calpos","Vardict_calref","Vardict_calalt","Vardict_calpos_best","Vardict_calref_best","Vardict_calalt_best")
input[,Vardict.complex[!Vardict.complex %in% colnames(input)]] <- NA


input$Mutect2_ROQ <- NA
input$Mutect2_samtools_DP <- NA



input$gnomAD_MAX.Stringent.005 <- (input$MAX_gnomAD_AF_VEP < 0.005 &
                                        input$MAX_gnomADe_AF_VEP < 0.005 &
                                        input$MAX_gnomADg_AF_VEP < 0.005)
input$gnomAD_MAX.lessStringent.005 <- (input$MAX_gnomAD_AF_VEP < 0.005 |
                                            input$MAX_gnomADe_AF_VEP < 0.005 |
                                            input$MAX_gnomADg_AF_VEP < 0.005)

# n_panmyeloid, myeloid_cosmic_count, isOncogenic, isTSG, WHY_CH in ANNOTATE


## ch_pd stuff
vars <- read.table(args$bolton_bick_vars, sep = "\t", header = T, comment.char = "")
vars$gene_aachange <- paste(vars$SYMBOL_VEP, vars$AAchange2, sep = "_")
vars$gene_cDNAchange <- paste(vars$SYMBOL_VEP, gsub(".*:","",vars$HGVSc_VEP), sep="_")
topmed.mutation.2 <- read.table(args$mut2_bick, sep = "\t", header = T, comment.char = "")
kelly.mutation.2 <- read.table(args$mut2_kelly, sep = "\t", header = T, comment.char = "")
##both Bick/Bolton
matches.2.c.p <- read.table(args$matches2, sep = "\t", header = T, comment.char = "")


input$AAchange <- gsub("(.*p\\.)(.*)", "\\2", input$HGVSp)
for (i in 1:length(AminoAcids)) {
  input$AAchange <- gsub(names(AminoAcids)[i], AminoAcids[i], input$AAchange)
}
input$gene_loci_p <- paste(input$SYMBOL_VEP,
                                  paste0(sapply(input$AAchange, function(x) str_split(x, "[0-9]+", n=2)[[1]][1]),
                                         as.numeric(str_extract(input$AAchange, "\\d+"))),
                                  sep = "_")
input$gene_loci_c <- paste(input$SYMBOL_VEP,
                                  gsub(".*:", "", input$HGVSc),
                                  sep = "_")
input$gene_loci_vep <- ifelse(is.na(input$gene_loci_p),input$gene_loci_c,input$gene_loci_p)
input$key <- with(input, paste(CHROM,POS,REF,ALT,sep = ":"))
input$gene_aachange <- with(input, paste(SYMBOL_VEP, AAchange, sep = "_"))
input$gene_cDNAchange <- paste(input$SYMBOL_VEP, gsub(".*:","",input$HGVSc_VEP), sep="_")

dims <- dim(input)[[1]]
input <- sqldf("SELECT l.*, r.`n.loci.vep`, r.`source.totals.loci`
            FROM `input` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
input <- input[!duplicated(input),]
## make sure aachange exists as in doesn't end with an '_'; example: DNMT3A_ for splice 
input <- sqldf("SELECT l.*, r.`n.HGVSp`, r.`source.totals.p`
            FROM `input` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.gene_aachange = r.gene_aachange) AND r.gene_aachange NOT LIKE '%_'")
input <- input[!duplicated(input),]
input <- sqldf("SELECT l.*, r.`n.HGVSc`, r.`source.totals.c`
            FROM `input` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.gene_cDNAchange = r.gene_cDNAchange) AND r.gene_cDNAchange NOT LIKE '%_'")
input <- input[!duplicated(input),]
dim(input)[[1]] == dims

input$n.loci.vep <- fillna(input$n.loci.vep, 0)
input$source.totals.loci <- fillna(input$source.totals.loci, 0)
input$n.HGVSp <- fillna(input$n.HGVSp, 0)
input$source.totals.p <- fillna(input$source.totals.p, 0)
input$n.HGVSc <- fillna(input$n.HGVSc, 0)
input$source.totals.c <- fillna(input$source.totals.c, 0)

vars.truncating <- vars[vars$truncating=="truncating",]
input <- sqldf("SELECT l.*, r.`truncating_loci_total`, r.`source.totals.loci.truncating`
            FROM `input` as l
            LEFT JOIN `vars.truncating` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep") # source.totals.loci
input <- input[!duplicated(input),]
input$truncating_loci_total <- fillna(input$truncating_loci_total, 0)
VAFS <- getVAFs(input)
input[,VAFS] <- apply(input[,VAFS], 2, function(x) fillna(x, 0))
input$max.under.0.35 <- apply(input[,VAFS],1,function(x) max(x, na.rm = T) < 0.35)




#########################################################################################################################
#########################################################################################################################

#'@returns dataframe with equal number of rows as input with binded counts columns
#'@param df dataframe: this is the MUTS dataframe; needs to have `HGVSp_VEP`,`var_key`,`sample` columns
#'@param sample.column string: column name for sample name/id
#'@example cosmic.run(MUTS, "SAMPLE", "var_key")
cosmic.run <- function(df, sample.column="SAMPLE", var_key.column="var_key") {
  # COSV52681947
  X1 <- split(df, df[,"CHROM"])
  ptm <- proc.time()
  X.1 <- parallel::mclapply(ls(X1), function(y) {
  # X.1 <- lapply(ls(X1), function(y) {
    print(y)
    df.chr <- X1[[y]]
    df.chr$HGVSp_VEP <- gsub(".*:","",df.chr$HGVSp_VEP)
    df.chr$Gene_HGVSp_VEP <- with(df.chr, paste(SYMBOL_VEP, HGVSp_VEP, sep = "_"))
    # df.chr$skey <- with(df.chr, paste(var_key, eid, sep = ":"))
    df.chr$skey <- paste(df.chr[[var_key.column]], df.chr[[sample.column]], sep = ":")
    
    cosmic <- read.table(paste0("/scratch1/fs1/bolton/UKBB/germline/already/cosmic/CosmicMutantExport.final.more.minimal.test.",y,".tsv"),
                         header = T, sep = "\t", comment.char = "", quote="")
                         # cosmic <- read.table(paste0("/scratch1/fs1/bolton/UKBB/germline/already/cosmic/CosmicMutantExport.final.more.minimal.test.",y,".tsv"),

    colnames(cosmic) <- c("COSMIC_ID","var_key","HGVSP","CosmicCount","heme_cosmic_count","myeloid_cosmic_count","Gene_HGVSp_VEP")
    # colnames(df.chr)[grep("key", colnames(df.chr))] <- "var_key"
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
    
    #cosmic.test$skey <- with(cosmic.test, paste(var_key, eid, sep = ":"))
    cosmic.test$skey <- paste(cosmic.test[[var_key.column]], cosmic.test[[sample.column]], sep = ":")
    row.names(cosmic.test) <- cosmic.test$skey
    # row.names(cosmic.test) <- cosmic.test$var_key
    ## put back into row order of df.chr
    cosmic.test <- cosmic.test[df.chr$skey,]
    # cosmic.test <- cosmic.test[df.chr$var_key,]
    
    if (all(cosmic.test[,c("CHROM","POS","REF","ALT",sample.column)]==df.chr[,c("CHROM","POS","REF","ALT",sample.column)])) {
      print(paste("good inner", y))
      return(cosmic.test)
    } else {
      print(paste("bad inner", y))
    }
  })
  proc.time() - ptm
  
  if (length(X1) == length(X.1)) { 
    names(X.1) <- ls(X1) 
  } else {
    stop("WRONG dimensions of df split by CHROM")
  }
  ## are all the names in the same order
  all.names <- all(names(X1)==names(X.1))
  all.match <- all(sapply(ls(X1), function(x) {
    dim(X1[[x]])[1] == dim(X.1[[x]])[1]
  }))
  if (all.match & all.names) {
    cosmic.final2 <- bind_rows(X.1)
  }
  if (all(cosmic.final2[,c("CHROM","POS","REF","ALT",sample.column)]==df[,c("CHROM","POS","REF","ALT",sample.column)])) {
    print("COSMIC good #################################################")
    MUTS <- cbind(df, cosmic.final2 %>% dplyr::select(COSMIC_ID,CosmicCount,heme_cosmic_count,myeloid_cosmic_count))
    return(MUTS)
  } else {
    print("COSMIC bad #################################################")
  }
}


# n_panmyeloid, myeloid_cosmic_count, isOncogenic, isTSG, WHY_CH in ANNOTATE
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
  
  
  # MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP","HGVSp_VEP","HGVSc_VEP","VariantClass","SAMPLE","EXON_VEP","AAchange.x",
  #                 "CosmicCount","heme_cosmic_count","myeloid_cosmic_count","eid")]
  MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP","HGVSp_VEP","HGVSc_VEP","VariantClass","SAMPLE","EXON_VEP","AAchange","SAMPLE")]
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
  
  
  ###########annotate with COSMIC########
  #source("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/cosmic.parallel.compute1.R")
  ## need to add a parallel
  ## MUTS MUST BE ORDERED by CHROM (and maybe POS? just to be safe?? Basically whatever column is paralleled needs to be sorted by...)
  # message("annotating variants with COSMIC...")
  # MUTS <- cosmic.run(MUTS, "SAMPLE", "var_key")
  # MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  # MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  # MUTS$myeloid_cosmic_count = fillna(MUTS$myeloid_cosmic_count, 0)
  # message("COSMIC done\n")
  
  
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
  message("oncoKB done\n")
  
  
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
  
  
  # MUTS2 = MUTS
  
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
  colnames(MUTS)[colnames(MUTS)=="SYMBOL_VEP"] <- c("Gene")
  MUTS.temp <- MUTS
  matched = list()
  
  message("annotating variants with PD Table...")
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
  message("PD Table done\n")
  
  ## easier to just not return anything from TC and call this below it
  MUTS = left_join(MUTS, ch_my_variants, by="var_key")
  MUTS$ch_my_pd = fillna(MUTS$ch_my_pd, 0)
  
  
  MUTS <- MUTS %>% 
    mutate(WHY_CH_ch_my_pd = ifelse(MUTS$ch_my_pd>0, "PD_table;", ""))
  
  
  # 6. Any variant occurring in the COSMIC “haematopoietic and lymphoid” category greater than or equal to 10 times
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(heme_cosmic_count >= 10 | myeloid_cosmic_count >= 5, paste0(WHY_CH_ch_my_pd, " Cosmic_heme;"), WHY_CH_ch_my_pd))
  # MUTS$ch_my_pd = ifelse(MUTS$heme_cosmic_count>=10,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  
  # 7. Any variant noted as potentially oncogenic in an in-house dataset of 7,000 individuals with myeloid neoplasm greater than or equal to 5 times
  panmyeloid_variant_counts = read.table(args$pan_myeloid, comment.char = "", sep = "\t", quote = "", header = T)
  panmyeloid_variant_counts$var_key = with(panmyeloid_variant_counts, paste(CHROM,POS,REF,ALT, sep = ":"))
  panmyeloid_variant_counts <- panmyeloid_variant_counts[panmyeloid_variant_counts$Annotation == "ONCOGENIC" | 
                                                           panmyeloid_variant_counts$Annotation == "SNP",]
  
  MUTS$Gene_HGVSp_VEP <- paste(MUTS$Gene, gsub(".*p.","", MUTS$HGVSp_VEP), sep = "_")
  panmyeloid_variant_counts$Gene_HGVSp_VEP <- paste(panmyeloid_variant_counts$SYMBOL_VEP, gsub(".*p.","", panmyeloid_variant_counts$HGVSp_VEP), sep = "_")
  
  message("annotating variants with panmyeloid...")
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
  message("panmyeloid done\n")
  
  MUTS <- MUTS %>%
    mutate(ch_my_pd = ifelse(n_panmyeloid>=5, ifelse(MUTS$ch_my_pd==2,2,1), MUTS$ch_my_pd),
           WHY_CH_ch_my_pd = ifelse(MUTS$n_panmyeloid>=5, paste0(WHY_CH_ch_my_pd, " Pan_Myeloid;"), WHY_CH_ch_my_pd))
  
  # MUTS$ch_my_pd = ifelse(MUTS$n_panmyeloid>=5,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  MUTS <- MUTS %>%
    mutate(ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", 1, 0),
           WHY_CH_ch_pd = ifelse(oncoKB=="Oncogenic" | oncoKB=="Likely Oncogenic" | oncoKB=="Predicted Oncogenic", "OncoKB;", ""))
  
  MUTS$isOncogenic <- MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic" | MUTS$oncoKB=="Predicted Oncogenic"
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
        # n_panmyeloid, myeloid_cosmic_count, isOncogenic, isTSG, WHY_CH in ANNOTATE
        # dplyr::select("CHROM","POS","REF","ALT","SAMPLE",
        #               "COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count",
        #               "MDS", "AML","MPN","n_panmyeloid",
        #               "oncoKB","isOncogenic","isTSG",
        #               "ch_my_pd","ch_pd","ch_pd2",
        #               "VariantClass","AAchange","Gene","WHY_CH") 
        dplyr::select("CHROM","POS","REF","ALT","SAMPLE",
                      "COSMIC_ID","CosmicCount","heme_cosmic_count","myeloid_cosmic_count",
                      "MDS", "AML","MPN","n_panmyeloid",
                      "isOncogenic","isTSG",
                      "ch_my_pd","ch_pd","ch_pd2",
                      "VariantClass","AAchange","Gene","WHY_CH") 
      
    )
  } else {
    message("Annotate PD failed: something wrong with dims")
  }
}


annotation <- annotate.PD(input)
setdiff(cols$V1, colnames(input))
input2 <-  left_join(input,
                    annotation,
                    by=c("CHROM","POS","REF","ALT","SAMPLE"))
setdiff(cols$V1, colnames(input2))
input2 <- input2[,cols$V1]
write.table(input2, output, sep = "\t", row.names = F, quote = F)

    