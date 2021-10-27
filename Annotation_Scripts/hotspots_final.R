#!/usr/local/bin/Rscript

library(argparser)
library(tidyverse)
library(stringr)
library(VariantAnnotation)
library(vcfR)
library(sqldf)


'%!in%' = function(x, y){
  ! ('%in%'(x, y))
}

fillna <- function(x, r) {
  x[is.na(x)] = r
  return(x)
}

parser2 <- arg_parser("Process NGS variant calls")
parser2 <- add_argument(parser2, "--impact-file", type="character", help="kelly's file")
parser2 <- add_argument(parser2, "--bick-file", type="character", help="bick's file")
parser2 <- add_argument(parser2, "--number-samples-annot", type="integer", help="number samples in loci for variant", default=5)



# args2 <- parse_args(parser2, list(c("-i","~/Bolton/data/hg38_mut_full_long_filtered_KB_deid_2.tsv"),
#                                 c("-b","~/Bolton/data/bick_topmed_variants.txt")))
args <- parse_args(parser2, list(c("-i","~/Bolton/data/hg38_mut_full_long_filtered_KB_deid_2.tsv"),
                                  c("-b","~/Bolton/data/bick_topmed_variants.txt")))                               


########################################################################################
AminoAcids = c('Cys'= 'C', 'Asp'= 'D', 'Ser'= 'S', 'Gln'= 'Q', 'Lys'= 'K',
               'Ile'= 'I', 'Pro'= 'P', 'Thr'= 'T', 'Phe'= 'F', 'Asn'= 'N',
               'Gly'= 'G', 'His'= 'H', 'Leu'= 'L', 'Arg'= 'R', 'Trp'= 'W',
               'Ala'= 'A', 'Val'='V', 'Glu'= 'E', 'Tyr'= 'Y', 'Met'= 'M',
               '%3D'='=', '='='=')
string <- "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|FrameshiftSequence|WildtypeProtein|gnomADe|gnomADe_AF|gnomADe_AF_AFR|gnomADe_AF_AMR|gnomADe_AF_ASJ|gnomADe_AF_EAS|gnomADe_AF_FIN|gnomADe_AF_NFE|gnomADe_AF_OTH|gnomADe_AF_SAS|gnomADg|gnomADg_AF|gnomADg_AF_ami|gnomADg_AF_oth|gnomADg_AF_afr|gnomADg_AF_sas|gnomADg_AF_asj|gnomADg_AF_fin|gnomADg_AF_amr|gnomADg_AF_nfe|gnomADg_AF_eas|clinvar|clinvar_CLINSIGN|clinvar_PHENOTYPE|clinvar_SCORE|clinvar_RCVACC|clinvar_TESTEDINGTR|clinvar_PHENOTYPELIST|clinvar_NUMSUBMIT|clinvar_GUIDELINES"
CSQnames <- str_split(string, "\\|")[[1]]
NUMBER = args2$number_samples_annot
## ch_pd stuff
mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- read.table(args2$impact_file,
                                                         quote = "", header = T, sep = "\t")
mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- mut_full_long_filtered_KB_deid.ch_pd2.hg38[mut_full_long_filtered_KB_deid.ch_pd2.hg38$VariantClass!="Silent",]
mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- mut_full_long_filtered_KB_deid.ch_pd2.hg38[with(mut_full_long_filtered_KB_deid.ch_pd2.hg38, order(mut_full_long_filtered_KB_deid.ch_pd2.hg38$CHROM,
                                                                                                     mut_full_long_filtered_KB_deid.ch_pd2.hg38$POS)),]
mut_full_long_filtered_KB_deid.ch_pd2.hg38$AAchange <- gsub("\\*", "Ter", mut_full_long_filtered_KB_deid.ch_pd2.hg38$AAchange)

bick.topmed <- read.table(args2$bick_file,
                          quote = "", header = T, sep = "\t")
bick.topmed$CHROM <- paste0("chr", bick.topmed$CHROM)
bick.topmed <- bick.topmed[with(bick.topmed, order(CHROM,POS)),]
bick.topmed.vep.p <- read.table("/Users/brian/Bolton/data/bick_topmed_variants.brian.pick_order.vep.annotated.tsv")
colnames(bick.topmed.vep.p) <- c("CHROM","POS","REF","ALT","CSQ")
bick.topmed.vep.p <- bick.topmed.vep.p %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
if (all(bick.topmed[,1:4]==bick.topmed.vep.p[,1:4])) {
  bick.topmed <- cbind(bick.topmed, bick.topmed.vep.p %>% dplyr::select(-c(CHROM,POS,ALT,REF)))
}
bick.topmed <- bick.topmed %>% 
  dplyr::rename(VariantClass=ExonicFunc)
bick.topmed$AAchange <- gsub("X","Ter",bick.topmed$AAchange)
bick.topmed$AAchange2 <- gsub("(.*p\\.)(.*)", "\\2", bick.topmed$HGVSp)
for (i in 1:length(AminoAcids)) {
  bick.topmed$AAchange2 <- gsub(names(AminoAcids)[i], AminoAcids[i], bick.topmed$AAchange2)
}
tb <- bick.topmed %>% 
  dplyr::select(CHROM,POS,REF,ALT,VariantClass,Consequence_VEP,Gene,SYMBOL_VEP,AAchange,AAchange2,HGVSp_VEP,HGVSc_VEP)

topmed.loci.n <- rbind(tb %>%
                         dplyr::filter(!grepl("c\\.", AAchange)) %>%
                         dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                         dplyr::group_by(Gene, loci) %>%
                         dplyr::mutate(n.loci=ifelse(is.na(loci),1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(Gene, AAchange) %>%
                         dplyr::mutate(n.aa=ifelse(is.na(AAchange),1,dplyr::n())) %>%
                         dplyr::ungroup() %>%
                         dplyr::mutate(loci.vep = ifelse(grepl('[A-Z]\\d+',AAchange2),str_extract(AAchange2, '[A-Z]\\d+'), AAchange2)) %>%
                         dplyr::group_by(SYMBOL_VEP, loci.vep) %>%
                         dplyr::mutate(n.loci.vep=ifelse(loci.vep=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, HGVSp_VEP) %>%
                         dplyr::mutate(n.HGVSp=ifelse(HGVSp_VEP=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>%
                         dplyr::group_by(SYMBOL_VEP, HGVSc_VEP) %>%
                         dplyr::mutate(n.HGVSc=ifelse(HGVSc_VEP=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::filter(n.loci.vep >= NUMBER) %>%
                         dplyr::select(CHROM,POS,REF,ALT,Gene,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,loci.vep,n.loci.vep,n.HGVSp,n.HGVSc),
                       tb %>%
                         dplyr::filter(grepl("c\\.", AAchange)) %>%
                         dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                         dplyr::group_by(Gene, loci) %>%
                         dplyr::mutate(n.loci=ifelse(is.na(loci),1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(Gene, AAchange) %>%
                         dplyr::mutate(n.aa=ifelse(is.na(AAchange),1,dplyr::n())) %>%
                         dplyr::ungroup() %>%
                         dplyr::mutate(loci.vep = ifelse(grepl('[A-Z]\\d+',AAchange2),str_extract(AAchange2, '[A-Z]\\d+'), AAchange2)) %>%
                         dplyr::group_by(SYMBOL_VEP, loci.vep) %>%
                         dplyr::mutate(n.loci.vep=ifelse(loci.vep=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, HGVSp_VEP) %>%
                         dplyr::mutate(n.HGVSp=ifelse(HGVSp_VEP=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>%
                         dplyr::group_by(SYMBOL_VEP, HGVSc_VEP) %>%
                         dplyr::mutate(n.HGVSc=ifelse(HGVSc_VEP=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::filter(n.loci.vep >= NUMBER) %>%
                         dplyr::select(CHROM,POS,REF,ALT,Gene,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,loci.vep,n.loci.vep,n.HGVSp,n.HGVSc)
) %>% arrange(desc(n.loci.vep),desc(n.HGVSp))



## kelly
mut.vep.mr <- read.table("/Users/brian/Bolton/data/mut_full_long_filtered_KB_deid.brian.pick_order.mane.rank.vep.annotated.tsv")
colnames(mut.vep.mr) <- c("CHROM","POS","REF","ALT","CSQ")
mut.vep.mr <- mut.vep.mr %>% separate(CSQ, paste0(CSQnames, "_VEP"), sep="\\|", extra = "merge", fill = "right")
if (all(mut_full_long_filtered_KB_deid.ch_pd2.hg38[,1:4]==mut.vep.mr[,1:4])) {
  mut_full_long_filtered_KB_deid.ch_pd2.hg38 <- cbind(mut_full_long_filtered_KB_deid.ch_pd2.hg38, mut.vep.mr %>% dplyr::select(-c(CHROM,POS,ALT,REF)))
}
all(mut_full_long_filtered_KB_deid.ch_pd2.hg38$Gene==mut_full_long_filtered_KB_deid.ch_pd2.hg38$SYMBOL_VEP)
tk <- mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>% 
  dplyr::select(CHROM,POS,REF,ALT,VariantClass,Consequence_VEP,Gene,SYMBOL_VEP,AAchange,HGVSp_VEP,cDNAchange,HGVSc_VEP)
tk$AAchange2 <- gsub("(.*p\\.)(.*)", "\\2", tk$HGVSp)
for (i in 1:length(AminoAcids)) {
  tk$AAchange2 <- gsub(names(AminoAcids)[i], AminoAcids[i], tk$AAchange2)
}

#####
# gr <- GRanges(seqnames = mut_full_long_filtered_KB_deid.ch_pd2.hg38$CHROM,
#               ranges = IRanges(start=mut_full_long_filtered_KB_deid.ch_pd2.hg38$POS,
#                                end = mut_full_long_filtered_KB_deid.ch_pd2.hg38$POS),
#               REF = mut_full_long_filtered_KB_deid.ch_pd2.hg38$REF,
#               ALT = mut_full_long_filtered_KB_deid.ch_pd2.hg38$ALT)
# final.vcf <- VariantAnnotation::VCF(rowRanges = gr, fixed = DataFrame(REF=DNAStringSet(gr$REF),
#                                                                       ALT=DNAStringSetList(as.list(gr$ALT)),
#                                                                       QUAL = 1))
# VariantAnnotation::writeVcf(final.vcf, "~/Bolton/data/mut_full_long_filtered_KB_deid.ch_pd2.hg38.vcf")
# system("gsed -i '1 i\\##fileformat=VCFv4.2' ~/Bolton/data/mut_full_long_filtered_KB_deid.ch_pd2.hg38.vcf")
# system("bgzip -f ~/Bolton/data/mut_full_long_filtered_KB_deid.ch_pd2.hg38.vcf")
# system("tabix -f ~/Bolton/data/mut_full_long_filtered_KB_deid.ch_pd2.hg38.vcf.gz")


####

kelly.loci.n <- rbind(tk %>%
                   dplyr::filter(!grepl("c\\.", AAchange)) %>%
                   dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                   dplyr::group_by(Gene, loci) %>%
                   dplyr::mutate(n.loci=ifelse(is.na(loci),1,dplyr::n())) %>%
                   dplyr::ungroup() %>% 
                   dplyr::group_by(Gene, AAchange) %>%
                   dplyr::mutate(n.aa=ifelse(is.na(AAchange),1,dplyr::n())) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(loci.vep = ifelse(grepl('[A-Z]\\d+',AAchange2),str_extract(AAchange2, '[A-Z]\\d+'), AAchange2)) %>%
                   dplyr::group_by(SYMBOL_VEP, loci.vep) %>%
                   dplyr::mutate(n.loci.vep=ifelse(loci.vep=="",1,dplyr::n())) %>%
                   dplyr::ungroup() %>% 
                   dplyr::group_by(SYMBOL_VEP, HGVSp_VEP) %>%
                   dplyr::mutate(n.HGVSp=ifelse(HGVSp_VEP=="",1,dplyr::n())) %>%
                   dplyr::ungroup() %>%
                   dplyr::group_by(SYMBOL_VEP, HGVSc_VEP) %>%
                   dplyr::mutate(n.HGVSc=ifelse(HGVSc_VEP=="",1,dplyr::n())) %>%
                   dplyr::ungroup() %>% 
                   dplyr::filter(n.loci.vep >= NUMBER) %>%
                   dplyr::select(CHROM,POS,REF,ALT,Gene,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,loci.vep,n.loci.vep,n.HGVSp,n.HGVSc),
                 tk %>%
                   dplyr::filter(grepl("c\\.", AAchange)) %>%
                   dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                   dplyr::group_by(Gene, loci) %>%
                   dplyr::mutate(n.loci=ifelse(is.na(loci),1,dplyr::n())) %>%
                   dplyr::ungroup() %>% 
                   dplyr::group_by(Gene, AAchange) %>%
                   dplyr::mutate(n.aa=ifelse(is.na(AAchange),1,dplyr::n())) %>%
                   dplyr::ungroup() %>%
                   dplyr::mutate(loci.vep = ifelse(grepl('[A-Z]\\d+',AAchange2),str_extract(AAchange2, '[A-Z]\\d+'), AAchange2)) %>%
                   dplyr::group_by(SYMBOL_VEP, loci.vep) %>%
                   dplyr::mutate(n.loci.vep=ifelse(loci.vep=="",1,dplyr::n())) %>%
                   dplyr::ungroup() %>% 
                   dplyr::group_by(SYMBOL_VEP, HGVSp_VEP) %>%
                   dplyr::mutate(n.HGVSp=ifelse(HGVSp_VEP=="",1,dplyr::n())) %>%
                   dplyr::ungroup() %>%
                   dplyr::group_by(SYMBOL_VEP, HGVSc_VEP) %>%
                   dplyr::mutate(n.HGVSc=ifelse(HGVSc_VEP=="",1,dplyr::n())) %>%
                   dplyr::ungroup() %>% 
                   dplyr::filter(n.loci.vep >= NUMBER) %>%
                   dplyr::select(CHROM,POS,REF,ALT,Gene,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,n.loci,n.loci.vep,n.aa,n.HGVSp,n.HGVSc)
) %>% arrange(desc(n.loci.vep),desc(n.HGVSp))

####
##bick
# gr <- GRanges(seqnames = bick.topmed$CHROM,
#               ranges = IRanges(start=bick.topmed$POS,
#                                end = bick.topmed$POS),
#               REF = bick.topmed$REF,
#               ALT = bick.topmed$ALT)
# final.vcf <- VariantAnnotation::VCF(rowRanges = gr, fixed = DataFrame(REF=DNAStringSet(gr$REF),
#                                                                       ALT=DNAStringSetList(as.list(gr$ALT)),
#                                                                       QUAL = 1))
# VariantAnnotation::writeVcf(final.vcf, "~/Bolton/data/bick_topmed_variants.vcf")
# system("gsed -i '1 i\\##fileformat=VCFv4.2' ~/Bolton/data/bick_topmed_variants.vcf")
# system("bgzip -f ~/Bolton/data/bick_topmed_variants.vcf")
# system("tabix -f ~/Bolton/data/bick_topmed_variants.vcf.gz")
####


topmed.mutation.1.p <- tb %>% 
  group_by(SYMBOL_VEP, AAchange2) %>% 
  mutate(n.aa=ifelse(AAchange2=="",0,dplyr::n()),
         gene_aachange = paste(SYMBOL_VEP,AAchange2,sep = "_")) %>%
  filter(n.aa>=1) %>%
  dplyr::select(CHROM,POS,REF,ALT,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,n.aa,AAchange2,HGVSp_VEP,gene_aachange)
length(unique(topmed.mutation.1.p$gene_aachange))
topmed.mutation.1.c <- tb %>% 
  mutate(cDNAchange=gsub(".*:","",HGVSc_VEP)) %>%
  group_by(SYMBOL_VEP, cDNAchange) %>%
  mutate(n.cDNA=ifelse(cDNAchange=="",0,dplyr::n()),
         gene_cDNAchange = paste(SYMBOL_VEP,cDNAchange,sep = "_")) %>%
  filter(n.cDNA>=1) %>%
  dplyr::select(CHROM,POS,REF,ALT,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,n.cDNA,AAchange2,HGVSp_VEP,gene_cDNAchange)
length(unique(topmed.mutation.1.c$gene_cDNAchange))
topmed.mutation.1.c.p <- unique(c(unique(topmed.mutation.1.p$gene_aachange), unique(topmed.mutation.1.c$gene_cDNAchange)))
length(topmed.mutation.1.c.p)
write.table(setNames(data.frame(topmed.mutation.1.c.p), "exact_match"),
            "~/Bolton/data/topmed.mutation.1.c.p.txt", quote = F, row.names = F)



kelly.mutation.1.p <- tk %>% 
  group_by(SYMBOL_VEP, AAchange2) %>% 
  mutate(n.aa=ifelse(AAchange2=="",0,dplyr::n()),
         gene_aachange = paste(SYMBOL_VEP,AAchange2,sep = "_")) %>%
  filter(n.aa>=1) %>%
  dplyr::select(CHROM,POS,REF,ALT,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,n.aa,AAchange2,HGVSp_VEP,gene_aachange)
length(unique(kelly.mutation.1.p$gene_aachange))
kelly.mutation.1.c <- tk %>% 
  mutate(cDNAchange=gsub(".*:","",HGVSc_VEP)) %>%
  group_by(SYMBOL_VEP, cDNAchange) %>%
  mutate(n.cDNA=ifelse(cDNAchange=="",0,dplyr::n()),
         gene_cDNAchange = paste(SYMBOL_VEP,cDNAchange,sep = "_")) %>%
  filter(n.cDNA>=1) %>%
  dplyr::select(CHROM,POS,REF,ALT,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,n.cDNA,AAchange2,HGVSp_VEP,gene_cDNAchange)
length(unique(kelly.mutation.1.c$gene_cDNAchange))
kelly.mutation.1.c.p <- unique(c(unique(kelly.mutation.1.p$gene_aachange), unique(kelly.mutation.1.c$gene_cDNAchange)))
write.table(setNames(data.frame(kelly.mutation.1.c.p), "exact_match"),
            "~/Bolton/data/kelly.mutation.1.c.p.txt", quote = F, row.names = F)






# vars.loci.n <- rbind(mut_full_long_filtered_KB_deid.ch_pd2.hg38 %>%
#                        dplyr::select(CHROM,POS,REF,ALT,Gene,AAchange,ch_my_pd,ch_pd2) %>%
#                        mutate(source="kelly"),
#                      bick.topmed %>%
#                        dplyr::select(CHROM,POS,REF,ALT,Gene,AAchange) %>%
#                        mutate(ch_my_pd = NA,
#                               ch_pd2 = NA,
#                               source = "bick"))


# need: AAchange, Gene, AAchange2, SYMBOL_VEP, HGVSp_VEP, HGVSc_VEP
# get: CHROM,POS,REF,ALT,Gene,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,n.loci,n.loci.vep,n.aa,n.HGVSp,n.HGVSc
vars.loci.n <- rbind(tk %>% dplyr::select(CHROM,POS,ALT,REF,VariantClass,Consequence_VEP,AAchange,Gene,AAchange2,SYMBOL_VEP,HGVSp_VEP,HGVSc_VEP) %>% mutate(source="kelly"),
                     tb %>% dplyr::select(CHROM,POS,ALT,REF,VariantClass,Consequence_VEP,AAchange,Gene,AAchange2,SYMBOL_VEP,HGVSp_VEP,HGVSc_VEP) %>%  mutate(source = "bick"))


vars.loci.n <- rbind(vars.loci.n %>%
                         dplyr::filter(!grepl("c\\.", AAchange)) %>%
                         dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                         dplyr::group_by(Gene, loci) %>%
                         dplyr::mutate(n.loci=ifelse(is.na(loci),1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, AAchange, source) %>%
                         dplyr::mutate(n.aachange.source=paste(source,ifelse(AAchange=="",1,dplyr::n()), sep = ": ")) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, AAchange) %>%
                         mutate(source.totals.p = paste(unique(str_split(paste0(n.aachange.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                         dplyr::mutate(n.aa=ifelse(is.na(AAchange),1,dplyr::n())) %>%
                         dplyr::ungroup() %>%
                         dplyr::mutate(loci.vep = ifelse(grepl('[A-Z]\\d+',AAchange2),str_extract(AAchange2, '[A-Z]\\d+'), AAchange2)) %>%
                         ## just to get protein counts from bick and kelly
                         dplyr::group_by(SYMBOL_VEP, loci.vep, source) %>%
                         dplyr::mutate(n.loci.vep.source=paste(source,ifelse(loci.vep=="",1,dplyr::n()), sep = ": ")) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, loci.vep) %>%
                         mutate(source.totals.loci = paste(unique(str_split(paste0(n.loci.vep.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                         dplyr::mutate(n.loci.vep=ifelse(loci.vep=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, HGVSp_VEP) %>%
                         dplyr::mutate(n.HGVSp=ifelse(HGVSp_VEP=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>%
                         ## just to get cdna counts from bick and kelly
                         dplyr::group_by(SYMBOL_VEP, HGVSc_VEP, source) %>%
                         dplyr::mutate(n.loci.vep.cdna.source=paste(source,ifelse(HGVSc_VEP=="",1,dplyr::n()), sep = ": ")) %>%
                         dplyr::ungroup() %>% 
                         dplyr::group_by(SYMBOL_VEP, HGVSc_VEP) %>%
                         mutate(source.totals.c = paste(unique(str_split(paste0(n.loci.vep.cdna.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                         dplyr::mutate(n.HGVSc=ifelse(HGVSc_VEP=="",1,dplyr::n())) %>%
                         dplyr::ungroup() %>% 
                         dplyr::filter(n.loci.vep >= NUMBER | n.HGVSc >=5) %>%
                         mutate(cDNAchange=gsub(".*:","",HGVSc_VEP),
                                gene_loci_vep = ifelse(loci.vep=="",paste(SYMBOL_VEP, cDNAchange, sep="_"),paste(SYMBOL_VEP, loci.vep, sep="_"))
                                # ,source.totals = ifelse(n.loci.vep < 5, source.totals.c, source.totals.p)
                                ) %>%
                         dplyr::select(CHROM,POS,REF,ALT,Gene,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,
                                       loci.vep,n.loci.vep,n.HGVSp,n.HGVSc,gene_loci_vep,source.totals.loci, source.totals.p, source.totals.c),
                       vars.loci.n %>%
                         dplyr::filter(grepl("c\\.", AAchange)) %>%
                       dplyr::mutate(loci = ifelse(grepl('[A-Z]\\d+',AAchange),str_extract(AAchange, '[A-Z]\\d+'), AAchange)) %>%
                       dplyr::group_by(Gene, loci) %>%
                       dplyr::mutate(n.loci=ifelse(is.na(loci),1,dplyr::n())) %>%
                       dplyr::ungroup() %>% 
                       dplyr::group_by(Gene, AAchange, source) %>%
                       dplyr::mutate(n.aachange.source=paste(source,ifelse(AAchange=="",1,dplyr::n()), sep = ": ")) %>%
                       dplyr::ungroup() %>% 
                       dplyr::group_by(Gene, AAchange) %>%
                       mutate(source.totals.p = paste(unique(str_split(paste0(n.aachange.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                       dplyr::mutate(n.aa=ifelse(is.na(AAchange),1,dplyr::n())) %>%
                       dplyr::ungroup() %>%
                       dplyr::mutate(loci.vep = ifelse(grepl('[A-Z]\\d+',AAchange2),str_extract(AAchange2, '[A-Z]\\d+'), AAchange2)) %>%
                       ## just to get protein counts from bick and kelly
                       dplyr::group_by(SYMBOL_VEP, loci.vep, source) %>%
                       dplyr::mutate(n.loci.vep.source=paste(source,ifelse(loci.vep=="",1,dplyr::n()), sep = ": ")) %>%
                       dplyr::ungroup() %>% 
                       dplyr::group_by(SYMBOL_VEP, loci.vep) %>%
                       # mutate(source.totals.p = paste(unique(str_split(paste0(n.loci.vep.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                       mutate(source.totals.loci = paste(unique(str_split(paste0(n.loci.vep.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                       dplyr::mutate(n.loci.vep=ifelse(loci.vep=="",1,dplyr::n())) %>%
                       dplyr::ungroup() %>% 
                       dplyr::group_by(SYMBOL_VEP, HGVSp_VEP) %>%
                       dplyr::mutate(n.HGVSp=ifelse(HGVSp_VEP=="",1,dplyr::n())) %>%
                       dplyr::ungroup() %>%
                       ## just to get cdna counts from bick and kelly
                       dplyr::group_by(SYMBOL_VEP, HGVSc_VEP, source) %>%
                       dplyr::mutate(n.loci.vep.cdna.source=paste(source,ifelse(HGVSc_VEP=="",1,dplyr::n()), sep = ": ")) %>%
                       dplyr::ungroup() %>% 
                       dplyr::group_by(SYMBOL_VEP, HGVSc_VEP) %>%
                       mutate(source.totals.c = paste(unique(str_split(paste0(n.loci.vep.cdna.source, collapse = ","),",")[[1]]),collapse=", ")) %>%
                       dplyr::mutate(n.HGVSc=ifelse(HGVSc_VEP=="",1,dplyr::n())) %>%
                       dplyr::ungroup() %>% 
                       dplyr::filter(n.loci.vep >= NUMBER | n.HGVSc >=5) %>%
                       mutate(cDNAchange=gsub(".*:","",HGVSc_VEP),
                              gene_loci_vep = ifelse(loci.vep=="",paste(SYMBOL_VEP, cDNAchange, sep="_"),paste(SYMBOL_VEP, loci.vep, sep="_"))
                              # ,source.totals = ifelse(n.loci.vep < 5, source.totals.c, source.totals.p)
                       ) %>%
                       dplyr::select(CHROM,POS,REF,ALT,Gene,SYMBOL_VEP,VariantClass,Consequence_VEP,HGVSc_VEP,AAchange,AAchange2,HGVSp_VEP,
                                     loci.vep,n.loci.vep,n.HGVSp,n.HGVSc,gene_loci_vep,source.totals.loci, source.totals.p, source.totals.c)
) %>% arrange(desc(n.loci.vep),desc(n.HGVSp))


# vars.loci.final <- vars.loci.n[!duplicated(vars.loci.n[,c("gene_loci_vep")]),]
vars <- vars.loci.n[!duplicated(vars.loci.n[,c("CHROM","POS","REF","ALT")]),]

vars$key <- with(vars, paste(CHROM,POS,REF,ALT, sep = ":"))
write.table(vars,"~/Bolton/data/bick.bolton.vars2.txt", row.names = F, quote = F, sep = "\t")
########################################################################################
########################################################################################

