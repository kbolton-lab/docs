#!/usr/local/bin/Rscript

library(argparser)
library(vcfR)
library(tidyverse)
library(stringr)
library(data.table)
library(httr)
library(rtracklayer)

################################################
### Read in VCF files

################################################

'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

parser <- arg_parser("Process NGS variant calls")
parser <- add_argument(parser, "--impact-file", type="character", help="kelly's file")
parser <- add_argument(parser, "--bick-file", type="character", help="bick's file")
parser <- add_argument(parser, "--cosmic-file", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--variant-calls-files", type="character", help="variant caller files, comma separated")
parser <- add_argument(parser, "--TSG-file", type="character", help="tumor suppressor genes")
parser <- add_argument(parser, "--oncoKB-curated", type="character", help="oncoKB all_curated_genes_v2.0.tsv")
parser <- add_argument(parser, "--pd-annotation-file", type="character", help="PD annoation files from papaemmanuil lab")
parser <- add_argument(parser, "--pan-myeloid", type="character", help="panmyeloid_variant_counts.tsv")


## test
args <- parse_args(parser, list(c("-i","~/Bolton/data/hg38_mut_full_long_filtered_KB_deid_2.tsv"),
                                c("-b","~/Bolton/data/bick_topmed_variants.txt"),
                                c("-c","~/Bolton/data/hg38_cosmic78_parsed.sorted.txt"),
                                c("-v", paste("~/Bolton/UKBB/results/10/1000144_23153_0_0.mutect.final.annotated.vcf.gz",
                                              ## add more variant caller files here that are VEP Annotated,
                                              sep=",")),
                                c("-T","~/Bolton/data/gene_census_TSG.txt"),
                                c("--oncoKB-curated","~/Bolton/data/all_curated_genes_v2.0.tsv"),
                                c("-p","~/Bolton/data/pd_table_kbreview_bick_trunc.tsv"),
                                c("--pan-myeloid","~/Bolton/data/panmyeloid_variant_counts.tsv")))


variant.files <- str_split(args$variant_calls_files, ",")[[1]]

mutect2_vcf_file <- variant.files[which(grepl("mutect", variant.files))]
# vardict_vcf_file <- variant.files[which(grepl("vardict", variant.files))]
       
### read in     
mutect2 <- vcfR2tidy(read.vcfR(mutect2_vcf_file, verbose = FALSE),
                     single_frame = TRUE,
                     info_types = TRUE,
                     format_types = TRUE)
mutect.df <- mutect2[["dat"]]
## remove columns where all NA
na.mutect.cols <- names(mutect.df[,apply(mutect.df, 2, function(x) all(is.na(x)))])
mutect.df <- mutect.df[,-which(names(mutect.df) %in% c(na.mutect.cols))]

# vardict <- vcfR2tidy(read.vcfR(vardict_vcf_file, verbose = FALSE),
#                      single_frame = TRUE,
#                      info_types = TRUE,
#                      format_types = TRUE)
# vardict.df <- vardict[["dat"]]

## remove blank CSQ rows (non-Coding) and separate
mutect.df <- mutect.df[!is.na(mutect.df$CSQ),]
## repeat for other callers
## We only ran VEP for Mutect2 on this dataset so cannot run script.

# VEP CSQ
string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames <- str_split(string, "\\|")[[1]]
# CSQ has N columns, maybe we should suffix cols with VEP?
mutect.df <- mutect.df %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
## repeat for other callers
## We only ran VEP for Mutect2 on this dataset so cannot run script.


################################################
### previous MSK and Topmed analysis
################################################
mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- read.table(args$impact_file,
                                                         quote = "", header = T, sep = "\t")
bick.topmed <- read.table(args$bick_file,
                          quote = "", header = T, sep = "\t")
bick.topmed$CHROM <- paste0("chr", bick.topmed$CHROM)
## number of samples with variant at loci
num = 5

## combine to count only
# https://community.rstudio.com/t/summarise-max-but-keep-all-columns/52449/3
topmed.loci.n <- rbind(bick.topmed %>%
                         dplyr::filter(!grepl("c\\.", AAchange)) %>%
                         dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                         dplyr::group_by(Gene, loci) %>%
                         dplyr::mutate(n=dplyr::n()) %>%
                         dplyr::ungroup() %>%
                         dplyr::filter(n >= num) %>%
                         #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                         dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n),
                       bick.topmed %>%
                         dplyr::filter(grepl("c\\.", AAchange)) %>%
                         dplyr::mutate(loci = str_extract(AAchange, "c\\.\\d+[\\-\\+]\\d[TCGA-]")) %>% # c.769+2T>G and -> is for insertion
                         dplyr::group_by(Gene, loci) %>%
                         dplyr::mutate(n=dplyr::n()) %>%
                         dplyr::ungroup() %>%
                         dplyr::filter(n >= num) %>%
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
                        dplyr::filter(!grepl("c\\.", AAchange)) %>%
                        dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                        dplyr::group_by(Gene, loci) %>%
                        dplyr::mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
                        dplyr::ungroup() %>%
                        dplyr::filter(n >= num) %>%
                        #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                        dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,VariantClass),
                      mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>%
                        dplyr::filter(grepl("c\\.", AAchange)) %>%
                        dplyr::mutate(loci = str_extract(AAchange, "c\\.\\d+[\\-\\+]\\d[TCGA-]")) %>% # c.769+2T>G and -> is for insertion
                        dplyr::group_by(Gene, loci) %>%
                        dplyr::mutate(n=ifelse(is.na(loci),1,dplyr::n())) %>%
                        dplyr::ungroup() %>%
                        dplyr::filter(n >= num) %>%
                        #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                        dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,VariantClass)
)
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
                dplyr::filter(!grepl("c\\.", AAchange)) %>%
                dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                dplyr::group_by(Gene, loci) %>%
                dplyr::mutate(n=dplyr::n()) %>%
                dplyr::ungroup() %>%
                dplyr::filter(n >= num) %>%
                #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,ch_my_pd,ch_pd2,source),
              vars.loci.n %>%
                dplyr::filter(grepl("c\\.", AAchange)) %>%
                dplyr::mutate(loci = str_extract(AAchange, "c\\.\\d+[\\-\\+]\\d[TCGA-]")) %>% # c.769+2T>G and -> is for insertion
                dplyr::group_by(Gene, loci) %>%
                dplyr::mutate(n=dplyr::n()) %>%
                dplyr::ungroup() %>%
                dplyr::filter(n >= num) %>%
                #distinct(CHROM,POS,REF,ALT, .keep_all = TRUE) %>%
                dplyr::select(CHROM,POS,REF,ALT,Gene,loci,AAchange,n,ch_my_pd,ch_pd2,source)) 

vars <- vars[!duplicated(vars[,c("CHROM","POS","REF","ALT")]),]

mutect.df<- dplyr::left_join(mutect.df, vars,
                             by = c("CHROM"="CHROM", "POS"="POS", 
                                    "REF"="REF", "ALT"="ALT"))



################################################
### annotate.PD
#' @param x this is a dataframe that is annotated with VEP.  
################################################

AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               'Ter'= 'X', '='='=', 'del'='DEL')

annotate.PD <- function(x) {

  MUTS <- x
  translate_consequence = c(
    "frameshift_variant" = "Frame_Shift_", 
    "stop_gained" = "Nonsense_Mutation", 
    "splice_region_variant" = "Splice_Site",
    "inframe_deletion" = "In_Frame_Del", 
    "inframe_insertion" = "In_Frame_Ins",
    "synonymous_variant" = "synonymous_variant",
    "missense_variant" = "Missense",
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
  
  MUTS$HGVSp_VEP <- gsub("(.*:p\\.)(.*)(fs.*)", "\\2", MUTS$HGVSp_VEP)
  MUTS$HGVSp_VEP <- gsub("(.*:p\\.)(.*)", "\\2", MUTS$HGVSp_VEP)
  MUTS$HGVSp_VEP <- gsub("del", "", MUTS$HGVSp_VEP)
  MUTS$HGVSp_VEP <- gsub("%3D", "=", MUTS$HGVSp_VEP) 
  MUTS <- MUTS[,c("CHROM","POS","REF","ALT","SYMBOL_VEP","HGVSp_VEP","VariantClass","SAMPLE","EXON_VEP")]
  
  MUTS$aa_ref <- sapply(MUTS$HGVSp_VEP, function(x) str_split(x, "[0-9]+", n=2)[[1]][1])
  MUTS$aa_alt <- sapply(MUTS$HGVSp_VEP, function(x) str_split(x, "[0-9]+", n=2)[[1]][2])
  MUTS$aa_pos <- as.numeric(str_extract(MUTS$HGVSp_VEP, "\\d+"))
  
  MUTS <- MUTS %>% 
    mutate(AAchange=paste0(AminoAcids[aa_ref],aa_pos,AminoAcids[aa_alt]))
  MUTS$var_key = paste(MUTS$CHROM, MUTS$POS, MUTS$REF, MUTS$ALT, sep = ":")
  
  ###create list of tumor suppressor genes####
  ## file 1
  gene_census = read_delim(args$TSG_file, delim = "\t", guess_max = 5e6)
  ## file 2
  oncoKB_curated = read_delim(args$oncoKB_curated, delim = "\t", guess_max = 5e6)
  oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
  #oncoKB_variants = read_delim("data/all_annotated_variants_v2.0.tsv", delim = "\t", guess_max = 5e6)
  
  TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
  
  
  ###########annotate with COSMIC########
  COSMIChg38 <- read.table(args$cosmic_file, sep = "\t", quote = "", header = T)
  MUTS = left_join(MUTS, COSMIChg38 %>% dplyr::select(COSMIC78_ID, CosmicCount, heme_cosmic_count, var_key),by="var_key")

  MUTS$CosmicCount = fillna(MUTS$CosmicCount, 0)
  MUTS$heme_cosmic_count = fillna(MUTS$heme_cosmic_count, 0)
  
  ##read in PD annoation files from papaemmanuil lab##
  ## file 4
  #A <- readxl::read_xlsx(args$pd_annotation_file)[,1:9] # now tsv
  #A <- read.table(args$pd_annotation_file, sep = "\t", header = T, quote = "")[,1:9]
  A <- read.table("~/Bolton/data/pd_table_kbreview_bick_trunc.tsv", sep = "\t", header = T, quote = "")[,1:9]
  A$keys = apply(A, 1, function(x) {
    names(x)[x != '*'] %>% .[!. %in% c('source', 'ch_my_pd', 'ch_pancan_pd')]
  })
  pd_dfs = split(A, paste(A$keys))
  
  panmyeloid_variant_counts = read_delim(args$pan_myeloid, delim = "\t", guess_max = 5e6)
  panmyeloid_variant_counts$var_key = paste0("chr", gsub("_", ":", panmyeloid_variant_counts$ID_VARIANT))
  
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
  
  MUTS = left_join(MUTS, panmyeloid_variant_counts %>% dplyr::select(Annotation,MDS,AML,MPN,var_key), by = 'var_key') %>%
    mutate(
      MDS = fillna(MDS, 0),
      AML = fillna(AML, 0),
      MPN = fillna(MPN, 0),
      Annotation = fillna(Annotation, '')
    ) %>%
    mutate(n_panmyeloid = ifelse(Annotation == 'ONCOGENIC', MDS + AML + MPN, 0))
  
  MUTS$ch_my_pd = ifelse(MUTS$n_panmyeloid>=5,ifelse(MUTS$ch_my_pd==2,2,1),MUTS$ch_my_pd)
  # We annotated variants as oncogenic (CH-PD) if they fulfilled any of the following criteria: #reffered to in code as ch_pancan_pd
  # 1. Any variant noted as oncogenic or likely oncogenic in OncoKB
  
  MUTS$ch_pd = ifelse(MUTS$oncoKB=="Oncogenic" | MUTS$oncoKB=="Likely Oncogenic", 1,0)
  
  # 2. Any truncating mutations (nonsense, essential splice site or frameshift indel) in known tumor suppressor genes as per the Cancer Gene Census or OncoKB.
  # Genes not listed in the cancer census or OncoKB were reviewed in the literature to determine if they were potentially tumor suppressor genes.# Annotate PD based on prevelance in cancer databases (COSMIC, Oncokb)
  
  
  MUTS$ch_pd = ifelse(MUTS$Gene %in% TSG &
                        MUTS$VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site"), 1, MUTS$ch_pd)
  
  #3. Any variant reported as somatic at least 20 times in COSMIC
  MUTS$CosmicCount <- as.numeric(MUTS$CosmicCount)
  MUTS$ch_pd = ifelse(MUTS$CosmicCount>=20, 1, MUTS$ch_pd)
  
  #4. Any variant meeting criteria for CH-Myeloid-PD as above.
  
  MUTS = MUTS %>% mutate(ch_pd = ifelse(ch_my_pd == 1 | ch_my_pd == 2, 1, ch_pd))
  #MUTS$DMP_PATIENT_ID = gsub("-T.*", "", MUTS$DMP_ASSAY_ID)
  
  annotate_PD.R = MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
  
  sample=grep("SAMPLE",colnames(x))
  
  if(all(annotate_PD.R[,c(1:4,8)]==x[,c(1:4,sample)])) {
    message("good same dimensions as original dataframe")
    return(
      annotate_PD.R %>% # annotate_PD.R is the MUTS %>% distinct(CHROM,POS,REF,ALT,SAMPLE, .keep_all = TRUE)
        dplyr::mutate(loci = str_extract(AAchange, "[A-Z]\\d+"),
                      AAchange = str_remove(AAchange, "p."),
                      gene_aachange = paste(Gene,AAchange,sep = "_"),
                      gene_loci = paste(Gene,loci,sep = "_"),
                      ch_pd2 = case_when(ch_pd==1 ~ 1,
                                              gene_loci %in% topmed.loci.n$gene_loci & gene_aachange %in% topmed.mutation.1$gene_aachange & VariantClass=="Missense_Mutation" ~ 1,
                                              gene_loci %in% topmed.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
                                              gene_loci %in% kelly.loci.n$gene_loci & gene_aachange %in% kelly.mutation.1$gene_aachange & VariantClass=="Missense_Mutation" ~ 1,
                                              gene_loci %in% kelly.loci.n$gene_loci & VariantClass %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site") ~ 1,
                                              TRUE ~ 0))
    )
  } else {
    message("something wrong with dimensions of new dataframe")
  }
}


## run function for each caller
annotation <- annotate.PD(mutect.df)
## join back to original dataframe
final.mutect.df <- left_join(mutect.df,
                       annotation %>%
                         dplyr::select("CHROM","POS","REF","ALT","SAMPLE","CosmicCount","heme_cosmic_count","MDS",
                                       "AML","MPN","ch_my_pd","ch_pd","ch_pd2","VariantClass","AAchange","Gene"),
                       by=c("CHROM","POS","REF","ALT","SAMPLE"))
## repeat for other callers
## We only ran VEP for Mutect2 on this dataset so cannot run script.

