####################################################################################################################
# Produce wide MSK-IMPACT CH mutation dataframe
####################################################################################################################
library(tidyverse)
suppressMessages(library(dplyr))
library(ggplot2)
library(ggsignif)

library(magrittr)
#library(dplyr)
library(ggplot2)
library(reshape2)
#library(stringr)
#library(readr)
library(RColorBrewer)
library(ggpubr)
library(gridGraphics)
library(ggsci)
library(ggsignif)
library(grid)
library(kableExtra)
library(ggrepel)
library(sjPlot)
library(sjlabelled)
library(sjmisc)
library(forcats)
#library(survminer)
library(ggplotify)
#library(compareGroups)
library(patchwork)
#library(geepack)
library(scales)
library(table1)
library(data.table)
#library(metafor)
#library(imputeTS)
library(lubridate)


#' max function that returns 0 when empty set is given
#' @param x a numerical vector, could be of length 0
#' @return max of x
max0 = function(x) {
  if (length(x) == 0) {
    return(0)
  } else {
    return(max(x))
  }
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


#' Produce wide MSK-IMPACT CH mutation dataframe
#' @param M MSK-IMPACT CH mutation dataframe
#' @param P clinical dataframe
#' @param out_file outputfile, NULL by default
#' @param boolean_mutations whether to 
#' @param negative_sample whether include CH-negative samples to the dataframe
#' @return wide CH mutation dataframe with indicator and max VAF columns
get_mutcount_wide = function(M, P, out_file = NULL, boolean_mutations = TRUE, negative_sample = TRUE,
                             bin_vaf = TRUE, verbose = TRUE) {

  # Create indicator matrix for patient ~ gene
  # ================================================================================================================
  
  # patients who got IMPACT3
  # v3_patients = P %>%
  #   filter(Version == 'v3') %>%
  #   pull(MRN) %>%
  #   unique
  
  GENE = M %>%
    # filter(ch_pancan_pd == 1) %>%
    #filter(CH_nonsilent == 1) %>%
    reshape2::dcast(
      formula = eid ~ SYMBOL_VEP,
      value.var = 'eid',
      fun.aggregate = function(x) {min(1, length(x))}
    )
  
  if (verbose) {cat('produced gene counts for', nrow(GENE), 'patients', fill = T)}

  # CH binary indicators
  # CH = M %>% group_by(MRN) %>%
  #   summarise(
  #     ch_nonmy = any(CH_my == 0 & CH_nonsilent == 1), # important this precedes modification of CH_my
  #     ch_nonmy_nonpd = any(CH_my == 0 & ch_pancan_pd == 0 & CH_nonsilent == 1),
  #     ch_nonpd = any(ch_pancan_pd == 0 & CH_nonsilent == 1),
  #     ch_my_nonpd = any(CH_my == 1 & ch_pancan_pd == 0 & CH_nonsilent == 1),
  #     ch_nonmy_pd = any(CH_my == 0 & ch_pancan_pd == 1 & CH_nonsilent == 1),
  #     ch_tch_my_pd = any(ch_my_pd == 1 & (Gene=="PPM1D" | Gene=="TP53" | Gene=="CHEK2") & CH_nonsilent == 1),
  #     ch_tch = any((Gene=="PPM1D" | Gene=="TP53" | Gene=="CHEK2") & CH_nonsilent == 1),
  #     CH_my = any(CH_my == 1 & CH_nonsilent == 1),
  #     ch_my_pd = any(ch_my_pd == 1 & CH_nonsilent == 1),
  #     ch_pancan_pd = any(ch_pancan_pd == 1 & CH_nonsilent == 1),
  #     CH_all = any(CH_all == 1),
  #     CH_nonsilent = any(CH_nonsilent == 1),
  #     CH_silent = any(CH_silent == 1)
  #   )

  # Calculate max vafs per patient for each CH subclass
  # ================================================================================================================
  # MAXVAF = M %>% group_by(MRN) %>%
  #   summarise(
  #     VAF_all = max0(VAF_N[CH_all == 1]),
  #     VAF_nonsilent = max0(VAF_N[CH_nonsilent == 1]),
  #     VAF_silent = max0(VAF_N[CH_silent == 1]),
  #     VAF_my = max0(VAF_N[CH_my == 1 & (Gene=="PPM1D" | Gene=="TP53" | Gene=="CHEK2")]),
  #     VAF_tch_my = max0(VAF_N[CH_my == 1]),
  #     VAF_nonmy = max0(VAF_N[CH_my == 0 & CH_nonsilent == 1]),
  #     VAF_my_pd = max0(VAF_N[ch_my_pd == 1]),
  #     VAF_tch_my_pd = max0(VAF_N[ch_my_pd == 1 & (Gene=="PPM1D" | Gene=="TP53" | Gene=="CHEK2")]),
  #     VAF_pancan_pd = max0(VAF_N[ch_pancan_pd == 1]),
  #     VAF_nonmy_pd = max0(VAF_N[ch_nonmy_pd == 1]))
  
  MAXVAF = M %>% dplyr::group_by(eid) %>%
    dplyr::summarise(VAF_all = max0(MAX_VAF))
  
  # M %>% dplyr::group_by(SAMPLE) %>%
  #   dplyr::summarise(VAF_all = max0(MAX_VAF))
  
  # Calculate number of mutations per patient for each CH subclass
  # ================================================================================================================
  # MUTNUM = M %>% group_by(MRN) %>% # add option to not count multiple mutations on same gene count as one
  #   summarise(
  #     mutnum_all = length(unique(Gene[CH_all == 1])),
  #     mutnum_nonsilent = length(unique(Gene[CH_nonsilent == 1])),
  #     mutnum_my = length(unique(Gene[CH_my == 1 & CH_nonsilent == 1])),
  #     mutnum_nonmy = length(unique(Gene[CH_my == 0 & CH_nonsilent == 1])),
  #     mutnum_silent = length(unique(Gene[CH_silent == 1])),
  #     mutnum_pancan_pd = length(unique(Gene[CH_nonsilent == 1 & ch_pancan_pd == 1])),
  #     mutnum_my_pd = length(unique(Gene[CH_nonsilent == 1 & ch_my_pd == 1])),
  #     mutnum_nonmy_pd = length(unique(Gene[CH_nonsilent == 1 & ch_my_pd == 0 & ch_pancan_pd == 1])))
  
  MUTNUM = M %>% dplyr::group_by(eid) %>% # add option to not count multiple mutations on same gene count as one
    dplyr::summarise(
      mutnum_all = length(SYMBOL_VEP)
    )
  if (verbose) {cat('produced mutation counts for', nrow(MUTNUM), 'patients', fill = T)}
  
  # Create indicator matrix for gene ontology functional classes
  # ================================================================================================================
  # ONTOLOGY = M %>%
  #   reshape2::dcast(
  #     formula = MRN ~ ontology,
  #     value.var = 'Gene',
  #     fun.aggregate = function(x) {min(1, length(x))})
    
  # Create Ontology Functional Class Max VAF
  # ================================================================================================================
  # ONTOVAF = M %>%
  #   reshape2::dcast(
  #     formula = MRN ~ ontology,
  #     value.var = 'VAF_N',
  #     fun.aggregate = max, 
  #     fill = 0) %>%
  #   setNames(c('MRN', paste('VAF', colnames(.)[-1], sep = '_')))
  
  # Calculate max VAF per gene (nonsyn mutations)
  # ================================================================================================================
  # GENEVAF = M %>% filter(CH_nonsilent == 1) %>% 
  #   reshape2::dcast(
  #     formula = MRN ~ Gene,
  #     value.var = 'VAF_N',
  #     fun.aggregate = max,
  #     fill = 0) %>%
  #   setNames(c('MRN', paste('VAF', colnames(.)[-1], sep = '_')))

  GENEVAF = M %>% 
    reshape2::dcast(
      formula = eid ~ SYMBOL_VEP,
      value.var = 'MAX_VAF',
      fun.aggregate = max,
      fill = 0) %>%
    setNames(c('eid', paste('VAF', colnames(.)[-1], sep = '_')))
  
  if (verbose) {cat('produced max VAFs for', nrow(GENEVAF), 'patients', fill = T)}

  # Mutation substitution types counts for all CH
  # ================================================================================================================
  substitution_types = c(
    "G_T" = "C_A", "G_C" = "C_G", "G_A" = "C_T",
    "A_T" = "T_A", "A_G" = "T_C", "A_C" = "T_G",
    "C_A" = "C_A", "C_G" = "C_G", "C_T" = "C_T",
    "T_A" = "T_A", "T_C" = "T_C", "T_G" = "T_G")
  
  SUB = M %>% filter(VARIANT_CLASS_VEP == 'SNV') %>%
    mutate(sub = substitution_types[paste(REF, ALT, sep = '_')]) %>%
    reshape2::dcast(
      formula = eid ~ sub,
      value.var = 'sub',
      #fun.aggregate = function(x) {min(1, length(x))})
      fun.aggregate = length)

  # Mutation substitution Max VAFs
  # ================================================================================================================
  SUBVAF = M %>% filter(VARIANT_CLASS_VEP == 'SNV') %>%
    mutate(sub = substitution_types[paste(REF, ALT, sep = '_')]) %>%
    reshape2::dcast(
      formula = eid ~ sub,
      value.var = "MAX_VAF",
      fun.aggregate = max,
      fill = 0) %>%
    setNames(c('eid', paste('VAF', colnames(.)[-1], sep = '_')))

  # Add MaxVAF for Hotspots
  # ================================================================================================================
  # hotspots = c("dnmt3a_r882h", "ppm1d_hotspot_1", "ppm1d_hotspot_2", "asxl1_e635fs15", 
  #               "asxl1_r693", "sf3b1_k700e", "gnas_r201h", "cbl_hotspot", "srsf2_95")
  # 
  # HOTSPOT = M[, c("MRN","VAF_N",hotspots)] %>% 
  #   mutate_at(.vars = hotspots, .funs = funs(VAF_N*.)) %>%
  #   dplyr::select(-VAF_N) %>%
  #   aggregate(by = list(M$MRN), FUN = max) %>% 
  #   dplyr::select(-Group.1) %>%
  #   setNames(c("MRN", "VAF_dnmt3a_r882h", "VAF_ppm1d_1", "VAF_ppm1d_2", "VAF_asxl1_e635fs15", "VAF_asxl1_r693",
  #              "VAF_sf3b1_k700e", "VAF_gnas_r201h", "VAF_cbl_hotspot", "VAF_srsf2_95"))

  # Merge mutation counts with MAXVAF and mutation coutns for CH subclasses
  # ================================================================================================================
  M_wide = Reduce(
    #x = list(GENE, GENEVAF, MAXVAF, MUTNUM, ONTOLOGY, ONTOVAF, SUB, SUBVAF, HOTSPOT, CH),
    x = list(GENE, GENEVAF, MAXVAF, MUTNUM, SUB, SUBVAF),
    f = function(x, y) {dplyr::full_join(x, y, by = 'eid')})
  
  if (verbose) {cat('made merged wide dataframe with dimension', nrow(M_wide), 'x', ncol(M_wide), fill = T)}

  # Convert numerical VAFs into bins
  # ================================================================================================================
#  if (bin_vaf) {
    
    VAF_columns = M_wide %>% colnames() %>% grep(pattern = "VAF", value = TRUE)
    
    M_wide <- M_wide %>% mutate_at(vars(all_of(VAF_columns)), 
                                   funs(bin=cut(
                                     .,
                                     breaks = c(-Inf, 0, 0.05, 0.10, 0.15, 0.20, 0.25,Inf),
                                     labels = c(0, 1, 2, 3, 4, 5, 6)
                                   )))
  
  #   M_wide = M_wide %>% mutate_at(
  #     # .predicate = function(x) {grepl(x = x, pattern = "VAF")},
  #     .funs = function(x) {
  #       cut(
  #         x,
  #         breaks = c(-Inf, 0, 0.05, 0.10, 0.15, 0.20, 0.35),
  #         labels = c(0, 1, 2, 3, 4, 5)
  #       )
  #     },
  #     .vars = VAF_columns
  #   )
  #   if (verbose) {cat('binned', length(VAF_columns), 'VAF columns', fill = T)}
  # }

  # Create mutation matrix with CH-negative patients
  # ================================================================================================================
  if (negative_sample) {
    M_wide = dplyr::full_join(M_wide, P, by = 'eid') %>%
      mutate_at(colnames(M_wide), function(x){tidyr::replace_na(x, 0)}) #%>%
      # mutate(
      #   # If sample run in version 3 make those genes missing 
      #   PPM1D = ifelse(MRN %in% v3_patients, NA, PPM1D),
      #   SRSF2 = ifelse(MRN %in% v3_patients, NA, SRSF2),
      #   VAF_PPM1D = ifelse(MRN %in% v3_patients, NA, VAF_PPM1D),
      #   VAF_SRSF2 = ifelse(MRN %in% v3_patients, NA, VAF_SRSF2)
      # )
    M_long <- dplyr::left_join(M, P, by = c('eid')) %>%
      mutate_at(colnames(M), function(x){tidyr::replace_na(x, 0)})
    if (verbose) {cat('merged with CH negative patients with', nrow(M_wide), 'total patients', fill = T)}
  
  } else {
    
    M_wide = dplyr::left_join(M_wide, P, by = 'MRN') %>%
      mutate_at(colnames(M_wide), function(x){tidyr::replace_na(x, 0)}) %>%
      mutate(
        # If sample run in version 3 make those genes missing 
        PPM1D = ifelse(MRN %in% v3_patients, NA, PPM1D),
        SRSF2 = ifelse(MRN %in% v3_patients, NA, SRSF2),
        VAF_PPM1D = ifelse(MRN %in% v3_patients, NA, VAF_PPM1D),
        VAF_SRSF2 = ifelse(MRN %in% v3_patients, NA, VAF_SRSF2)
      )
    
    if (verbose) {cat('merged with clinical data for', nrow(M_wide), 'CH positive patients', fill = T)}
  }
  
  if (!is.null(out_file)) {
    if (verbose) {cat('writing final', nrow(M_wide), 'x', ncol(M_wide), 'dataframe to', out_file, fill = T)}
    write.table(M_wide, out_file, sep = "\t", row.names = F, quote = T)
  }
  
  return(M_wide)
}


######################
source("~/Bolton/R_stuff/toolbox.R")
# source("~/Bolton/HIV/metadata.R")

M <- read.table("/Users/brian/Bolton/UKBB/ch_pd2_transplant_KB_BW_FINAL.txt", sep = "\t", quote = "", header = T, comment.char = "")
M$ch_pd_final <- fillna(M$ch_pd_final,0)
M <- M[as.logical(M$ch_pd_final),]
M$eid <- as.integer(gsub("_.*","",M$eid))
getVAFs(M)
M$MAX_VAF <- apply(M[,getVAFs(M)],1,max)
M$MIN_VAF <- apply(M[,getVAFs(M)],1,min)
P <- read.table("~/Bolton/UKBB/transplant_match.txt")
P$transplant <- as.numeric(P$transplant)
P <- P[!duplicated(P$eid),]


M_wide <- get_mutcount_wide(M, P, out_file = NULL, boolean_mutations = TRUE, negative_sample = TRUE,
                         bin_vaf = TRUE, verbose = TRUE)
M <- M %>% left_join(P,by="eid")




panel_theme = theme_bw() + theme(
  panel.border = element_blank(),
  legend.position = "none",
  panel.grid.minor = element_blank(),
  plot.subtitle = element_text(hjust = 0.5, size = 8),
  plot.title = element_text(face = 'bold', size = 12, hjust = 0, vjust = -11),
  panel.grid.major = element_blank(),
  strip.background = element_blank(),
  strip.text = element_text(size = 6),
  axis.text.y = element_text(size = 6),
  axis.text.x = element_text(size = 6),
  axis.title = element_text(size = 8),
  axis.line = element_line(),
  plot.margin = unit(c(0,0,0,0), 'pt')
) 

cancer <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbCancerRegistry_09232021.rds")
lung <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbLungCancerCaCo_10132021.rds")
ukbb <- readRDS("/Volumes/yin.cao/Active/UKBBGenomicData/ukbbExome_09232021.rds")
## transplant
sum(M$eid[M$transplant] %in% cancer$eid)/length(M$eid[M$transplant]) # with CH 100%
sum(P$eid[as.logical(P$transplant)] %in% cancer$eid)/length(P$eid[as.logical(P$transplant)]) # total 41%
library(reshape2) # for melt
reshape2::melt(P,"transplant", value.name = "ageBaseline")
## no transplant
sum(M$eid[!M$transplant] %in% cancer$eid)/length(M$eid[!M$transplant]) # with CH 
sum(P$eid[!as.logical(P$transplant)] %in% cancer$eid)/length(P$eid[!as.logical(P$transplant)])


# therapy_colors = c('#0072B2','#D55E00')
therapy_colors = c('#0072B2','#D55E00','#B32D00')
# ggplot(P %>% group_by(ageBaseline, transplant) %>% summarise(n=n()), aes(fill=transplant, y=n, x=ageBaseline)) +
#   geom_bar(position="dodge", stat="identity")
# ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25)





# plot(density(M_wide$ageBaseline))
# abline(v=median(M_wide$ageBaseline))
transplant.meta <- read.csv("~/Bolton/UKBB/clinicaldata_ukbb_transplant.csv")
table(transplant.meta$transplant_type)
transplant.meta <- transplant.meta %>%
  mutate(transplant_type_class = case_when(transplant_type == "kidney" ~ "kidney",
                                           transplant_type == "liver" ~ "liver",
                                           transplant_type == "muscle" ~ "muscle",
                                           transplant_type == "lung" ~ "lung",
                                           TRUE ~ "Other"))
table(transplant.meta$transplant_type_class)
# pre before post
transplant.meta <- transplant.meta[order(transplant.meta$timing_blood, decreasing = T),] %>%
  group_by(eid) %>%
  dplyr::slice_max(eid, n=1, with_ties = F)
M_wide <- M_wide %>% left_join(transplant.meta %>% dplyr::select(eid, timing_blood), by="eid")
M_wide$timing_blood <- fillna(M_wide$timing_blood, "No Transplant")
M <- M %>% left_join(transplant.meta %>% dplyr::select(eid, timing_blood), by="eid")
M$timing_blood <- fillna(M$timing_blood, "No Transplant")
M$transplant = ifelse(M$timing_blood=="No Transplant","No Transplant",ifelse(M$timing_blood=="pre","Blood Draw Pre-Transplant","Blood Draw Post-Transplant"))

M_wide2 <- M_wide %>% 
  mutate(mutnum_all_r = case_when(mutnum_all == 0 ~ 0,
                                  mutnum_all == 1 ~ 1,
                                  mutnum_all >= 2 ~ 2),
         # transplant = ifelse(transplant==1,"Transplant","No Transplant"),
         transplant = ifelse(timing_blood=="No Transplant","No Transplant",ifelse(timing_blood=="pre","Blood Draw Pre-Transplant","Blood Draw Post-Transplant")),
         case_control_binary = relevel(factor(transplant), ref = "No Transplant"),
         age_bins = cut(x=ageBaseline, breaks=c(40,45,50,55,60,65,Inf), include.lowest = T),
         age_bins2 = case_when(
           age_bins == "[40,45]" ~ 1,
           age_bins == "(45,50]" ~ 2,
           age_bins == "(50,55]" ~ 3,
           age_bins == "(55,60]" ~ 4,
           age_bins == "(60,65]" ~ 5,
           age_bins == "(65,Inf]" ~ 6),
         ch_pd_binary = factor(ifelse(mutnum_all>0,1,0)),
         CH_all = ifelse(ch_pd_binary==1, 'CH+', 'CH-'),
         `CH+`=ifelse(CH_all=="CH+",1,0),
         `CH-`=ifelse(CH_all=="CH-",1,0)
  )

# df1 <- data.frame(M_wide$ageBaseline[M_wide$transplant=="Transplant"])
# colnames(df1) <- "Transplant"
# data <- reshape2::melt(df1)
# df2 <- data.frame(M_wide$ageBaseline[M_wide$transplant=="No Transplant"])
# colnames(df2) <- "No Transplant"
# data2 <- reshape2::melt(df2)
# 
# ggplot(data,
#        aes(x=value, fill=variable)) + geom_density(alpha=0.25) +
#   geom_density(aes(x=value, fill=variable), data=data2, alpha=0.25) +
#   ggtitle("Age Distributions") +
#   scale_fill_manual(name="Transplant Status",values = therapy_colors) +
#   theme(legend.title = element_text(size = 8), 
#         legend.text = element_text(size = 8)) +
#   theme(legend.key.size = unit(0.4, 'cm')) +
#   geom_vline(xintercept=mean(data$value), size=0.5, color="black") +
#   geom_vline(xintercept=mean(data2$value), size=0.5, color="black")

age_groups <- unique(M_wide2$age_bins)
font_size = 10
# age_curve_theme = 
#   theme(
#     legend.position = 'top',
#     legend.key.size = unit(5, 'mm'),
#     legend.title = element_blank(),
#     legend.direction = 'horizontal',
#     plot.title = element_text(hjust = -0.08),
#     axis.text.x = element_text(angle = 45, vjust = 0.5, size = font_size),
#     axis.text.y = element_text(size = font_size),
#     axis.title = element_text(size = font_size),
#     legend.text = element_text(size = font_size)
#   )

gene_list <- unique(M$SYMBOL_VEP)


## for M_long

M2 <- M
gene_list <- unique(M2$SYMBOL_VEP)
M_long2 <- M2 %>% 
  mutate(
         transplant = ifelse(timing_blood=="No Transplant","No Transplant",ifelse(timing_blood=="pre","Blood Draw Pre-Transplant","Blood Draw Post-Transplant")),
         case_control_binary = relevel(factor(transplant), ref = "No Transplant"),
         age_bins = cut(x=ageBaseline, breaks=c(40,45,50,55,60,65,Inf), include.lowest = T),
         age_bins2 = case_when(
           age_bins == "[40,45]" ~ 1,
           age_bins == "(45,50]" ~ 2,
           age_bins == "(50,55]" ~ 3,
           age_bins == "(55,60]" ~ 4,
           age_bins == "(60,65]" ~ 5,
           age_bins == "(65,Inf]" ~ 6),
         case_control_binary = relevel(factor(transplant), ref = "No Transplant")
  )

# tally
n_Cases.trans.pre.blood = M_wide2 %>% dplyr::count(transplant) %>% filter(transplant == 'Blood Draw Pre-Transplant') %>% pull(n)
n_Cases.trans.post.blood = M_wide2 %>% dplyr::count(transplant) %>% filter(transplant == 'Blood Draw Post-Transplant') %>% pull(n)
n_Control = M_wide2 %>% dplyr::count(transplant) %>% filter(transplant == 'No Transplant') %>% pull(n)
D = M %>%
  reshape2::dcast(
    formula = SYMBOL_VEP + transplant ~ .,
    value.var = 'eid',
    fun.aggregate = function(eid) {length(unique(eid))}) %>%
  dplyr::rename("n_patient" = ".") %>%
  mutate(
    prop_patient = case_when(
      transplant == "No Transplant" ~ n_patient/n_Control,
      transplant == "Blood Draw Pre-Transplant" ~ n_patient/n_Cases.trans.pre.blood,
      transplant == "Blood Draw Post-Transplant" ~ n_patient/n_Cases.trans.post.blood)) %>%
  filter(SYMBOL_VEP %in% gene_list) %>%
  arrange(SYMBOL_VEP)

D2 = M_wide %>%
  reshape2::dcast(
    formula = mutnum_all_r + transplant ~ .,
    value.var = 'eid',
    fun.aggregate = function(eid) {length(unique(eid))}) %>%
  dplyr::rename("n_patient" = ".") %>%
  mutate(prop_patient = case_when(
    transplant == "Transplant" ~ n_patient/n_Cases,
    transplant == "No Transplant" ~ n_patient/n_Control)) %>%
  arrange(desc(mutnum_all_r))

# ggplot(D2, aes(x = reorder(SYMBOL_VEP,-prop_patient), 
D2$mutnum_all_r[1]="2+"
D2 <- 
ggplot(D2, aes(x = reorder(mutnum_all_r,prop_patient), 
               y = prop_patient, fill = as.factor(transplant))) +
  geom_bar(stat = 'identity', position = "dodge", color = 'black', size = 0.25) +
  # panel_theme +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
        axis.line = element_line(colour = "black"), legend.title = element_blank(),
        legend.key.size = unit(5, 'mm'), legend.position = 'top', legend.direction = 'horizontal',
        axis.title = element_text(size = 12), 
        # axis.text.x = element_text(hjust = 1, size = 12),
        legend.text = element_text(size = 12))  +
  ylab("Proportion with CH Mutations") +
  xlab('Mutation Number') +
  scale_fill_manual(values = therapy_colors) +
  scale_color_manual(values = therapy_colors) +
  geom_text(aes(label=round(prop_patient,3)), position=position_dodge(width=0.9), vjust=-0.25) +
  annotate(
    geom = "text", x = 0.8, y = 0.16, 
    label = "(No Transplant)", hjust = 0, vjust = 1, size = 3
  )

# D2 = M_long2 %>% 
#   reshape2::dcast(
#     formula = SYMBOL_VEP  ~ .,
#     value.var = 'eid',
#     fun.aggregate = function(SAMPLE) {length(unique(SAMPLE))}
#   ) %>%
#   dplyr::rename("n_patient" = ".") %>%
#   mutate(
#     prop_patient = case_when(
#       case.ASCVD.score == 'CVD.High.ASCVD' ~ n_patient/n_Cases_High,
#       case.ASCVD.score == 'CVD.Low.ASCVD' ~ n_patient/n_Cases_Low,
#       case.ASCVD.score == 'No CVD' ~ n_patient/n_Control
#     )
#   ) %>%
#   filter(SYMBOL_VEP %in% gene_list) %>%
#   mutate(
#     SYMBOL_VEP = factor(SYMBOL_VEP, gene_list),
#     case.ASCVD.score = factor(case.ASCVD.score, c('No CVD', 'CVD.Low.ASCVD', 'CVD.High.ASCVD'))
#   ) %>%
#   arrange(SYMBOL_VEP)

library(ggsignif)
library(HSAUR2)
library(survival)
## may want to use kelly's asterick's formula cause I used conditional logit
## see Review_Code
# asterisks = sapply(factor(gene_list),
#  function(gene) {
#    model = clogit(
#      formula = as.formula(paste0(gene," ~ case_control_binary + strata(match_id)")),
#      data = M_wide
#    )
#    gene_pval = model %>% summary %>% coefficients %>% .['case_control_binaryTransplant', 'Pr(>|z|)']
#    gene_qval = p.adjust(gene_pval, method = 'fdr', n = length(gene_list))
#    return(signif.num(gene_qval, ns = F))
#  }
# )
asterisks = sapply(factor(gene_list),
 function(gene) {
   model = clogit(
     formula = as.formula(paste0(gene," ~ case_control_binary + strata(match_id)")),
     data = M_wide
   )
   gene_pval = model %>% summary %>% coefficients %>% .['case_control_binaryTransplant', 'Pr(>|z|)']
   gene_qval = p.adjust(gene_pval, method = 'fdr', n = length(gene_list))
   return(signif.num(gene_qval, ns = F))
 }
)

high.match_id <- M_wide2[M_wide2$case.ASCVD.score=="CVD.High.ASCVD","match_id"]
low.match_id <- M_wide2[M_wide2$case.ASCVD.score=="CVD.Low.ASCVD","match_id"]
asterisks.low = sapply(factor(gene_list), 
                       function(gene) {
                         model = clogit(
                           formula = as.formula(paste0(gene," ~ case_control_binary + ch_pd_binary + race_b + strata(match_id)")),
                           data = M_wide2[M_wide2$match_id %in% low.match_id & M_wide2$case.ASCVD.score != "CVD.High.ASCVD",]
                         )
                         gene_pval = model %>% summary %>% coefficients %>% .['case_control_binaryCVD', 'Pr(>|z|)']
                         gene_qval = p.adjust(gene_pval, method = 'fdr', n = length(gene_list))
                         return(signif.num(gene_qval, ns = F))
                       }
)

asterisks.high = sapply(factor(gene_list), 
                        function(gene) {
                          model = clogit(
                            formula = as.formula(paste0(gene," ~ case_control_binary + ch_pd_binary + race_b + strata(match_id)")),
                            data =  M_wide2[M_wide2$match_id %in% high.match_id & M_wide2$case.ASCVD.score != "CVD.Low.ASCVD",]
                          )
                          gene_pval = model %>% summary %>% coefficients %>% .['case_control_binaryCVD', 'Pr(>|z|)']
                          gene_qval = p.adjust(gene_pval, method = 'fdr', n = length(gene_list))
                          return(signif.num(gene_qval, ns = F))
                        }
) 




do_plot = function(p, f, w, h, r = 300, save_pdf = T) {
  ggsave(f, plot = p, width = w, height = h, dpi = r)
  if (save_pdf) {
    ggsave(paste0(str_remove(f, '\\..+'), '.pdf'), plot = p, width = w, height = h, dpi = r)
  }
  knitr::include_graphics(f)
}

sort(table(M$SYMBOL_VEP), decreasing = T)[1:10]
ggplot(data.frame(sort(table(M$SYMBOL_VEP), decreasing = T)[1:10]), aes(Var1,Freq)) +
  geom_bar(stat = "identity", fill="#CC0000") + 
  ggtitle("CH for top Genes") +
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  xlab('') 


  
  
###plots
p_hist = ggplot(D, aes(x = reorder(SYMBOL_VEP,-prop_patient), 
                       y = prop_patient, fill = as.factor(transplant))) +
  geom_bar(stat = 'identity', position = "dodge", color = 'black', size = 0.25) +
  panel_theme +
  theme(panel.grid.major = element_blank(), panel.border = element_blank(),
    axis.line = element_line(colour = "black"), legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'), legend.position = 'top', legend.direction = 'horizontal',
    axis.title = element_text(size = font_size), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = font_size),
    legend.text = element_text(size = font_size)) +
  #annotate('text', x = gene_list, y = 0.03, label = asterisks, size = 4) +
  ylab("Proportion with mutated Gene") +
  xlab('') +
  scale_fill_manual(values = therapy_colors) +
  scale_color_manual(values = therapy_colors)
p_hist
p_hist + labs(title = 'A')



df <- M_wide %>%
  dplyr::select(mutnum_all, age_bins, case_control_binary) %>%
  dplyr::group_by(age_bins, case_control_binary) %>%
  dplyr::summarise(prop = sum(mutnum_all)/length(mutnum_all),
            n=sum(mutnum_all),
            total=length(mutnum_all),
            lower = prop.test(n, total, conf.level = .95)$conf.int[1],
            upper = prop.test(n, total, conf.level = .95)$conf.int[2]
  ) %>%
  dplyr::mutate(lab=paste0(n,"(",round(prop,2),"%",")"))

df.single <- M_wide %>%
  dplyr::select(mutnum_all, age_bins, case_control_binary) %>%
  dplyr::group_by(age_bins) %>%
  dplyr::summarise(prop = sum(mutnum_all)/length(mutnum_all),
                   n=sum(mutnum_all),
                   total=length(mutnum_all),
                   lower = prop.test(n, total, conf.level = .95)$conf.int[1],
                   upper = prop.test(n, total, conf.level = .95)$conf.int[2]
  ) %>%
  dplyr::mutate(lab=paste0(n,"(",round(prop,2),"%",")"))


p_ribbon <- ggplot(df, aes(x=age_bins, y=prop, ymin = lower, ymax = upper,fill=case_control_binary, group=case_control_binary)) +
  geom_ribbon(alpha = 0.2, colour=NA) +
  geom_line(aes(color=case_control_binary)) +
  geom_point(aes(color=case_control_binary)) +
  xlab("Age") +
  ylab("Proportion with CH") +
  scale_fill_manual(values = therapy_colors) +
  scale_color_manual(values = therapy_colors) +
  panel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'),
    legend.position = 'top',
    legend.direction = 'horizontal',
    axis.title = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    legend.text = element_text(size = font_size)
  ) 
p_ribbon + labs(title = 'B')

M_wide3 <- M_wide %>% 
  mutate(mutnum_all_r = case_when(mutnum_all == 0 ~ 0,
                                  mutnum_all == 1 ~ 1,
                                  mutnum_all >= 2 ~ 2),
         transplant = ifelse(transplant==1,"Transplant","No Transplant"),
         case_control_binary = relevel(factor(transplant), ref = "No Transplant"),
         age_bins = cut(x=ageBaseline, breaks=c(40,42.5,45,47.5,50,52.5,55,57.5,60,62.5,65,67.5,Inf), include.lowest = T),
         age_bins2 = case_when(
           age_bins == "[40,42.5]" ~ 1,
           age_bins == "[42.5,45]" ~ 1,
           age_bins == "(45,47.5]" ~ 2,
           age_bins == "(47.5,50]" ~ 2,
           age_bins == "(50,52.5]" ~ 3,
           age_bins == "(52.5,55]" ~ 3,
           age_bins == "(55,57.5]" ~ 4,
           age_bins == "(57.5,60]" ~ 4,
           age_bins == "(60,62.5]" ~ 5,
           age_bins == "(62.5,65]" ~ 5,
           age_bins == "(65,67.5]" ~ 6,
           age_bins == "(67.5,Inf]" ~ 6),
         ch_pd_binary = factor(ifelse(mutnum_all>0,1,0)),
         CH_all = ifelse(ch_pd_binary==1, 'CH+', 'CH-'),
         `CH+`=ifelse(CH_all=="CH+",1,0),
         `CH-`=ifelse(CH_all=="CH-",1,0)
  )

df.single2 <- M_wide3 %>%
  dplyr::select(mutnum_all, age_bins, case_control_binary) %>%
  dplyr::group_by(age_bins) %>%
  dplyr::summarise(prop = sum(mutnum_all)/length(mutnum_all),
                   n=sum(mutnum_all),
                   total=length(mutnum_all),
                   lower = prop.test(n, total, conf.level = .95)$conf.int[1],
                   upper = prop.test(n, total, conf.level = .95)$conf.int[2]
  ) %>%
  dplyr::mutate(lab=paste0(n,"(",round(prop,2),"%",")"))


ggplot(df.single, aes(x=age_bins, y=prop, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2, colour=NA, group = 1) +
  geom_line(group = 1) +
  geom_point() +
  xlab("Age") +
  ylab("Proportion with CH") +
  # scale_fill_manual(values = therapy_colors) +
  # scale_color_manual(values = therapy_colors) +
  # panel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'),
    legend.position = 'top',
    legend.direction = 'horizontal',
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 10)
  )

ggplot(df.single2, aes(x=age_bins, y=prop, ymin = lower, ymax = upper)) +
  geom_ribbon(alpha = 0.2, colour=NA, group = 1) +
  geom_line(group = 1) +
  geom_point() +
  xlab("Age") +
  ylab("Proportion with CH") +
  # scale_fill_manual(values = therapy_colors) +
  # scale_color_manual(values = therapy_colors) +
  # panel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'),
    legend.position = 'top',
    legend.direction = 'horizontal',
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    legend.text = element_text(size = 10)
  )



  
p_vaf_mutation <- ggplot(M %>%
                           mutate(case_control = ifelse(transplant==0,"No Transplant","Transplant"))
                         , aes(x=case_control, y=MAX_VAF, fill=as.factor(case_control))) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank()) +
  ylab("VAF by mutation") +
  scale_fill_manual(name="CVD Status",values = therapy_colors) +
  scale_color_manual(name="CVD Status",values = therapy_colors) +
  panel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'),
    #legend.position = 'top',
    #legend.justification = "top",
    legend.direction = 'horizontal',
    axis.title = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    axis.title.x = element_blank(),
    legend.text = element_text(size = font_size)
  )
p_vaf_mutation




p_vaf_sample <- ggplot(M_wide %>% filter(mutnum_all>0), aes(x=case_control_binary, y=VAF_all, fill=case_control_binary)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x=element_blank()) +
  ylab("MAX VAF / sample") +
  scale_fill_manual(name="CVD Status",values = therapy_colors) +
  scale_color_manual(name="CVD Status",values = therapy_colors) +
  panel_theme +
  theme(
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'),
    #legend.position = 'top',
    #legend.direction = 'horizontal',
    axis.title.x = element_blank(),
    axis.title = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    legend.text = element_text(size = font_size)
  )
p_vaf_sample




tally_variable = function(D, baseline_cols, variable) {
  if (variable == 'Total') {
    res = cbind(
      D %>% reshape2::dcast(. ~ CH_all, fun.aggregate = length, value.var = 'eid') %>% dplyr::select(-.),
      D %>% reshape2::dcast(. ~ case_control_binary, fun.aggregate = length, value.var = 'eid') %>% dplyr::select(-.)
    ) %>%
      dplyr::mutate(Total = rowSums(.[1:2])) %>%
      dplyr::mutate(variable = 'Total', category = 'Total') %>%
      dplyr::select(baseline_cols) %>%
      dplyr::mutate(ref = "")
  } else {
    ref = levels(D[[variable]])[1]
    if (is.null(ref)) {
      ref = ""
    }
    res = cbind(
      D %>% reshape2::dcast(paste0(variable, ' ~ CH_all'), fun.aggregate = length, value.var = 'eid') %>%
        mutate(variable = variable) %>% rename(category = !!variable),
      D %>% reshape2::dcast(paste0(variable, ' ~ ch_pd_binary'), fun.aggregate = length, value.var = 'eid'),
      D %>% dplyr::count(get(variable)) %>% dplyr::select(n) %>% rename(Total = n)
    ) %>%
      dplyr::select(baseline_cols) %>%
      dplyr::arrange(category == 'Missing') %>% # sorting
      dplyr::mutate(ref = ref)
  }
  res = res %>%
    dplyr::mutate(`CH-` = paste0(`CH-`, ' (', signif(`CH-`* 100/Total, 2), '%)')) %>%
    dplyr::mutate(`CH+` = paste0(`CH+`, ' (', signif(`CH+`* 100/Total, 2), '%)'))
  return(res)
}
baseline_cols = c('variable', 'category', 'CH-', 'CH+', 'Total')

## table
display_kable = function(kable, file = 'test.pdf')  {
  
  kable %>% save_kable(file = file)
  kable %>% as_image()
}
# D <- M_wide %>% mutate(AGE=age_bins)
# rbind(D %>% tally_variable(baseline_cols = baseline_cols, variable = 'Total'),
# D %>% tally_variable(baseline_cols = baseline_cols, variable = 'gender'),
# D %>% tally_variable(baseline_cols = baseline_cols, variable = 'AGE'),
# D %>% tally_variable(baseline_cols = baseline_cols, variable = 'race_b')
# ) %>%
#   mutate(variable = format_variable(variable)) %>%
#   select(-ref) %>%
#   kable(format = "latex", booktabs = T) %>%
#   kable_styling(
#     latex_options = c("hold_position"),
#     full_width = T,
#     font_size = 10
#   ) %>%
#   column_spec(1, width = 10) %>%
#   column_spec(2:5, width = 70) %>%
#   collapse_rows(columns = 1:2, row_group_label_position = 'stack') %>%
#   knitr::display_kable('table1.1.pdf')




# M_wide_pairs_only <- M_wide %>%
#   group_by(`Case-control\r\nmatching ID`) %>%
#   dplyr::mutate(n=dplyr::n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n >= 2)
# 
# df2_pairs_only <- M_wide_pairs_only %>%
#   dplyr::select(mutnum_all, case_control_binary, ch_pd_binary, CH_all) %>%
#   dplyr::group_by(CH_all,case_control_binary) %>%
#   dplyr::summarise(n = dplyr::n()) %>%
#   dplyr::mutate(prop = n / sum(n))


################################################################################
library(ggsignif)
library(HSAUR2)
library(survival)



## with race
model = clogit(
  formula = I(case_control_binary=="Transplant") ~ ch_pd_binary + strata(match_id),
  data = M_wide
) 
pval <- round(model %>% summary %>% coefficients %>% .['ch_pd_binary1', 'Pr(>|z|)'],3)



# df3 <- M_wide %>%
#   dplyr::select(mutnum_all, case_control_binary, ch_pd_binary, CH_all) %>%
#   dplyr::group_by(CH_all,case_control_binary) %>%
#   dplyr::summarise(n = dplyr::n()) %>%
#   dplyr::mutate(prop = sum(n)) %>%
#   dplyr::mutate(prop = n / sum(n))
df3 <- M_wide %>%
  dplyr::select(mutnum_all, case_control_binary, ch_pd_binary, CH_all) %>%
  dplyr::group_by(case_control_binary, CH_all) %>%
  dplyr::summarise(n = dplyr::n()) %>%
  dplyr::mutate(prop =  n / sum(n))


dummy <- ggplot(df3 %>% filter(CH_all=="CH+"),
                aes(x=CH_all, y=prop, fill=case_control_binary)) + 
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(name="Transplant Status",values = therapy_colors) +
  scale_color_manual(name="Transplant Status",values = therapy_colors) +
  ylab("Proportion") + 
  ylim(c(0,0.05)) +
  panel_theme + 
  theme(
    panel.grid.major = element_blank(), 
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.title = element_blank(),
    legend.key.size = unit(5, 'mm'),
    legend.position = 'top',
    legend.direction = 'horizontal',
    axis.title = element_text(size = font_size),
    axis.text.x = element_text(size = font_size),
    legend.text = element_text(size = font_size),
    axis.title.x = element_text(size=8, hjust=0)
  ) + 
  xlab("*Conditional logistic regression stratified by matching set")
  # geom_signif(y_position=0.25, xmin=0.75, xmax=1.25, annotation=paste0("p = ",pval,"*"), tip_length=0.03,
  #             size=0.35, textsize = 3.0)

dummy + labs(title="D")

# df2 <- M_wide %>%
#   dplyr::select(mutnum_all, case_control_binary, ch_pd_binary, CH_all) %>%
#   dplyr::group_by(case_control_binary, CH_all) %>%
#   dplyr::summarise(n = dplyr::n()) %>%
#   dplyr::mutate(prop = n / sum(n))
# 
# M_wide_pairs_only <- M_wide %>%
#   group_by(`Case-control_matching_ID`) %>%
#   dplyr::mutate(n=dplyr::n()) %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(n >= 2)
# 
# df2_pairs_only <- M_wide_pairs_only %>%
#   dplyr::select(mutnum_all, case_control_binary, ch_pd_binary, CH_all) %>%
#   dplyr::group_by(CH_all,case_control_binary) %>%
#   dplyr::summarise(n = dplyr::n()) %>%
#   dplyr::mutate(prop = n / sum(n))
library(patchwork)
panel = ((p_hist + labs(title = 'A')) / (p_ribbon + labs(title = 'B')) |
           (((p_vaf_mutation + labs(title = 'C')) | p_vaf_sample) + plot_layout(guides="collect")) / (dummy + labs(title="D")))

do_plot(panel, "~/Bolton/UKBB/fig2_tranplant.png", 10, 6, save_pdf = T)
