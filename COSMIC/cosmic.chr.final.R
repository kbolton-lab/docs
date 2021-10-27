#!/usr/local/bin/Rscript


library(tidyverse)

i <- commandArgs(trailingOnly = TRUE)

myeloid <- c("acute_myeloid_leukaemia",
             "acute_myeloid_leukaemia_therapy_related",
             "myelodysplastic_syndrome",
             "chronic_myelomonocytic_leukaemia",
             "polycythaemia_vera",
             "acute_leukaemia_of_ambiguous_lineage",
             "essential_thrombocythaemia", 
             "myelofibrosis",
             "chronic_myeloid_leukaemia",
             "acute_leukaemic_transformation_of_myeloproliferative_neoplasm",
             "mast_cell_neoplasm",
             "acute_myeloid_leukaemia_associated_with_MDS",
             "myelodysplastic_syndrome_therapy_related",
             "blast_phase_chronic_myeloid_leukaemia",
             "myelodysplastic-myeloproliferative_neoplasm-unclassifiable",
             "acute_leukaemic_transformation_of_primary_myelofibrosis",
             "myeloproliferative_neoplasm",
             "juvenile_myelomonocytic_leukaemia",
             "blastic_plasmacytoid_dendritic_cell_neoplasm",
             "acute_leukaemia_of_ambiguous_lineage_associated_with_MDS",
             "myelodysplastic-myeloproliferative_neoplasm",
             "granulocytic_sarcoma",
             "chronic_eosinophilic_leukaemia-hypereosinophilic_syndrome",
             "acute_leukaemic_transformation_of_polycythaemia_vera",
             "myeloid_neoplasm_unspecified_therapy_related",
             "acute_basophilic_leukaemia",
             "chronic_neutrophilic_leukaemia",
             "myeloproliferative_neoplasm_unclassifiable",
             "acute_leukaemic_transformation_of_essential_thrombocythaemia",
             "acute_leukaemic_transformation_of_chronic_myelomonocytic_leukaemia",
             "chronic_myelomonocytic_leukaemia_therapy_related")


cosmic94 <- readr::read_tsv(paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.chr",i,".tsv"),
                            col_names = c("gene name",
                                          "tumor_id",
                                          "Accession Number",
                                          "Primary site",
                                          "Primary histology",
                                          "Histology subtype 1",
                                          "Histology subtype 2",
                                          "Histology subtype 3",
                                          "GENOMIC_MUTATION_ID",
                                          "Mutation CDS",
                                          "Mutation AA",
                                          "Mutation Description",
                                          "GRCh",
                                          "Mutation genome position",
                                          "Mutation somatic status",
                                          "HGVSP",
                                          "HGVSC",
                                          "HGVSG"))

## remove duplicate [tumor_id, COSMIC ID, and HGVSC id}
## this will remove the transcript based duplicates but still keep all transcript annotions and HVGSP in case we need to join on a different ENSP protein
cosmic94.t <- cosmic94[!duplicated(cosmic94[,c("tumor_id","GENOMIC_MUTATION_ID","HGVSC")]),]

cosmic94.t <- cosmic94.t %>%
  group_by(GENOMIC_MUTATION_ID,`Accession Number`, .drop = F) %>%
  mutate(cosmic_count_chr=dplyr::n(),
         haematopoietic_and_lymphoid_tissue_count_chr=sum(`Primary site` == "haematopoietic_and_lymphoid_tissue"),
         myeloid_count_chr=sum(`Histology subtype 1` %in% myeloid)) %>%
  dplyr::select(GENOMIC_MUTATION_ID, #2
         `Accession Number`, #3
         tumor_id,
         HGVSG,
         HGVSC,
         HGVSP,
         `Mutation genome position`,
         cosmic_count_chr,
         haematopoietic_and_lymphoid_tissue_count_chr,
         myeloid_count_chr) %>% distinct()
cosmic94.t$CHROM <- paste0("chr",i)
if (i==1) {
  write.table(cosmic94.t[,c(11,1:10)], paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.chr",i,".final.parsed.test.tsv"), col.names=T,
              row.names=F, sep = "\t", quote = F)
} else {
  write.table(cosmic94.t[,c(11,1:10)], paste0("/storage1/fs1/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.chr",i,".final.parsed.test.tsv"), col.names=F,
              row.names=F, sep = "\t", quote = F)
}


# cat files
# cat $(ls CosmicMutantExport.chr*.final.parsed.test.tsv | sort -V) | awk -F'\t' '{print $9,$1,$2,$3,$4,$5,$6,$7,$8}' OFS='\t' > CosmicMutantExport.final.parsed.test.all.chrs.tsv
# cat $(ls CosmicMutantExport.chr*.final.parsed.test.tsv | sort -V) > CosmicMutantExport.final.parsed.test.all.chrs.tsv

## doing this in all with awk now
# cosmic.final <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.final.parsed.test.all.chrs.tsv")
# cosmic.final$Start <- gsub(".*:", "", cosmic.final$Mutation.genome.position)
# cosmic.final$Start <- gsub("-.*", "", cosmic.final$Start)
# cosmic.final$End <- gsub(".*-", "", cosmic.final$Mutation.genome.position)
# cosmic.final$Ref <- gsub("[0-9]*:g.[0-9]*([ACTG])>([ACTG])", "\\1", cosmic.final$HGVSG)
# cosmic.final$Alt <- gsub("[0-9]*:g.[0-9]*([ACTG])>([ACTG])", "\\2", cosmic.final$HGVSG)
# cosmic.final$var_key <- with(cosmic.final, paste(CHROM, Start, Ref, Alt, sep = ":"))
#   
# test <- cosmic.final %>% filter(var_key=="chr17:7673802:C:T")



