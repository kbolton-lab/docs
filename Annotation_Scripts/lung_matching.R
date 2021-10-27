#!/usr/local/bin/Rscript

library(argparser)
library(dplyr)
library(data.table)

parser2 <- arg_parser("Filter NSamples")
parser2 <- add_argument(parser2, "--number", type="integer", help="filter out > number")
parser2 <- add_argument(parser2, "--column", type="character", help="which column to group by, usually column with <chrom:pos:ref:alt>")
parser2 <- add_argument(parser2, "--input", type="character", help="in")
parser2 <- add_argument(parser2, "--output", type="character", help="out")
args2 <- parse_args(parser2)
args2 <- parse_args(parser2, list(c("-n","27"),
                                c("-c","key")))


# args2$input <- "/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.tsv"
args2$input <- "/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.tsv"
# args2$input <- "~/Bolton/HIV/final_all_rows_all_filters.tsv"
args2$output <- "/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.nsamp.tsv"

transplant.match.output.passed.pon <- fread(file=args2$input, header = T, sep = "\t", quote="")
transplant.match.output.passed.pon$key <- with(transplant.match.output.passed.pon, paste(CHROM,POS,REF,ALT,sep=":"))

total.passed.pon.nsamp <- transplant.match.output.passed.pon %>% 
  dplyr::group_by(key) %>% 
  dplyr::mutate(num=dplyr::n()) %>%
  ungroup() %>%
  filter(num<=args2$number)
  

write.table(total.passed.pon.nsamp, args2$output, sep = "\t", row.names = F, quote = F)

total.passed.pon.nsamp1 <- fread(file="/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.nsamp.tsv",
                                header = T, sep = "\t", quote="")
total.passed.pon.nsamp2 <- fread(file="/Volumes/bolton/Active/projects/mocha/UKBB/exome/test/test_annotate/combined/combined.final.reannotate.pon_pass.nsamp.updated.sorted.tsv",
                                header = T, sep = "\t", quote="")
total.passed.pon.nsamp2$pan_my <- apply(total.passed.pon.nsamp2[,c("MDS","AML","MPN")], 1, sum)
hotspots <- total.passed.pon.nsamp2[total.passed.pon.nsamp2$called & # called by both callers
                                      !total.passed.pon.nsamp2$Mutect2_PASS & # didn't pass mutect
                                      (total.passed.pon.nsamp2$heme_cosmic_count >= 10 | # either heme count >= 10
                                         total.passed.pon.nsamp2$myeloid_cosmic_count >= 5 | # or myeloid count >=5
                                         total.passed.pon.nsamp2$n.loci.vep >= 5 | # or bick/bolton loci >= 5
                                         total.passed.pon.nsamp2$pan_my >= 5),] # or pan myeloid >=5

transplant.ch_pd2 <- total.passed.pon.nsamp2[(total.passed.pon.nsamp2$ch_pd2 & # annotated ch_pd2 from the annotate_PD script
                           total.passed.pon.nsamp2$Mutect2_PASS & # Passed by Mutect
                           total.passed.pon.nsamp2$alt_strand_counts_min_1_caller_only & # Passed Mutect Strand Bias, at least >10% or <90% on both strands
                           total.passed.pon.nsamp2$max.over.0.02 & # At least 2% VAF
                           (total.passed.pon.nsamp2$max.under.0.35 | # Under 35% VAF 
                              (!total.passed.pon.nsamp2$max.under.0.35 & (!is.na(total.passed.pon.nsamp2$n.loci.vep)))) & # or if over 35% is bb hotspot
                           !total.passed.pon.nsamp2$Vardict_PON_2AT2_percent & # Passed PON2
                           !total.passed.pon.nsamp2$Mutect2_PON_2AT2_percent),] # Passed PON2
transplant.ch_pd2 <- transplant.ch_pd2[, c(1:252,272,270,271,257,256,259:269,253:255,258)]
hotspots <- hotspots[, c(1:252,272,270,271,257,256,259:269,253:255,258)]

# for (i in 1:length(colnames(transplant.ch_pd2))) { print(paste0(i,": ",colnames(transplant.ch_pd2)[i]))}

transplant.ch_pd2 <- sqldf("SELECT l.*, r.`source.totals`
            FROM `transplant.ch_pd2` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
transplant.ch_pd2 <- transplant.ch_pd2[!duplicated(transplant.ch_pd2),]
transplant.ch_pd2$n.loci.vep<- fillna(transplant.ch_pd2$n.loci.vep, 0)
transplant.ch_pd2 <- sqldf("SELECT l.*, r.`n.HGVSp`
            FROM `transplant.ch_pd2` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSp_VEP = r.HGVSp_VEP)")
transplant.ch_pd2 <- transplant.ch_pd2[!duplicated(transplant.ch_pd2),]
transplant.ch_pd2 <- sqldf("SELECT l.*, r.`n.HGVSc`
            FROM `transplant.ch_pd2` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSc_VEP = r.HGVSc_VEP)")
transplant.ch_pd2 <- transplant.ch_pd2[!duplicated(transplant.ch_pd2),]

hotspots <- sqldf("SELECT l.*, r.`source.totals`
            FROM `hotspots` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR l.gene_loci_vep = r.gene_loci_vep")
hotspots <- hotspots[!duplicated(hotspots),]
hotspots$n.loci.vep<- fillna(hotspots$n.loci.vep, 0)
hotspots <- sqldf("SELECT l.*, r.`n.HGVSp`
            FROM `hotspots` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSp_VEP = r.HGVSp_VEP)")
hotspots <- hotspots[!duplicated(hotspots),]
hotspots <- sqldf("SELECT l.*, r.`n.HGVSc`
            FROM `hotspots` as l
            LEFT JOIN `vars` as r
            on l.key = r.key OR (l.HGVSc_VEP = r.HGVSc_VEP)")
hotspots <- hotspots[!duplicated(hotspots),]


transplant.ch_pd2$gnomAD_MAX.Stringent.005 <- (transplant.ch_pd2$MAX_gnomAD_AF_VEP < 0.005 &
                                                 transplant.ch_pd2$MAX_gnomADe_AF_VEP < 0.005 &
                                                 transplant.ch_pd2$MAX_gnomADg_AF_VEP < 0.005)
hotspots$gnomAD_MAX.Stringent.005 <- (hotspots$MAX_gnomAD_AF_VEP < 0.005 &
                                        hotspots$MAX_gnomADe_AF_VEP < 0.005 &
                                        hotspots$MAX_gnomADg_AF_VEP < 0.005)

## if False means all 3 groups has a Max VAF > .0007 so this is less stringent for filtering out because
## all 3 needs to have a Max above .005
transplant.ch_pd2$gnomAD_MAX.lessStringent.005 <- (transplant.ch_pd2$MAX_gnomAD_AF_VEP < 0.005 |
                                                     transplant.ch_pd2$MAX_gnomADe_AF_VEP < 0.005 |
                                                     transplant.ch_pd2$MAX_gnomADg_AF_VEP < 0.005)
hotspots$gnomAD_MAX.lessStringent.005 <- (hotspots$MAX_gnomAD_AF_VEP < 0.005 |
                                            hotspots$MAX_gnomADe_AF_VEP < 0.005 |
                                            hotspots$MAX_gnomADg_AF_VEP < 0.005)




get_oncokb = function(mut) {
  
  apiKey = "a83627cd-47e4-4be0-82dd-8f4cc9e4d6d0"
  request_url = paste0("https://oncokb.org//api/v1/annotate/mutations/byProteinChange?hugoSymbol=",
                       mut['SYMBOL_VEP'], "&alteration=", mut['AAchange'], "&consequence=", mut['VariantClass'])
  
  res=httr::content(httr::GET(request_url, httr::add_headers(Authorization = paste("Bearer", apiKey))))
  return(res$oncogenic)
}
h = curl::new_handle()
curl::handle_setopt(h, http_version = 2)
httr::set_config(httr::config(http_version = 0))
cl = parallel::makeCluster(8)
transplant.ch_pd2$oncoKB <- parallel::parApply(cl, transplant.ch_pd2, 1, get_oncokb)
hotspots$oncoKB <- parallel::parApply(cl, hotspots, 1, get_oncokb)
parallel::stopCluster(cl)
transplant.ch_pd2$isOncogenic <- transplant.ch_pd2$oncoKB=="Oncogenic" | transplant.ch_pd2$oncoKB=="Likely Oncogenic"
hotspots$isOncogenic <- hotspots$oncoKB=="Oncogenic" | hotspots$oncoKB=="Likely Oncogenic"

###create list of tumor suppressor genes####
## file 1
gene_census = read.table(args$TSG_file, sep="\t", header=T)
## file 2
oncoKB_curated = read.table(args$oncoKB_curated, sep="\t", header=T)
oncoKB_curated.tsg = oncoKB_curated %>% filter(Is_Tumor_Suppressor=="Yes")
TSG = c(gene_census$Gene,oncoKB_curated.tsg$Hugo_Symbol) %>% unique()
transplant.ch_pd2$isTSG <- transplant.ch_pd2$SYMBOL_VEP %in% TSG
hotspots$isTSG <- hotspots$SYMBOL_VEP %in% TSG

ct.heme <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.high.heme.tsv", 
                      sep = "\t", quote = "", header = T)
ct.myeloid <- read.table("/Volumes/bolton/Active/projects/annotation_files/cosmic/CosmicMutantExport.high.myeloid.tsv", 
                         sep = "\t", quote = "", header = T)
ct <- rbind(ct.heme,ct.myeloid) %>% group_by(GENOMIC_MUTATION_ID) %>% mutate(n=n()) %>% distinct()
ct$gene <- gsub("_.*", "", ct$Gene_HGVSp_VEP)
ct$aa.pos <- gsub(".*_", "", ct$Gene_HGVSp_VEP)
ct$aa.pos <- as.numeric(str_extract(ct$aa.pos, "\\d+"))
transplant.ch_pd2$aa.pos <- as.numeric(str_extract(gsub(".*_", "",transplant.ch_pd2$gene_loci_vep), "\\d+"))
hotspots$aa.pos <- as.numeric(str_extract(gsub(".*_", "",hotspots$gene_loci_vep), "\\d+"))
colnames(ct)[2] <- "key"
ct$CHROM.POS <- unlist(lapply(ct$key, function(x) paste(str_split(x,":")[[1]][1],str_split(x,":")[[1]][2],sep = ":")))
ct$GENE.AA.POS <- with(ct, paste(gene, aa.pos, sep=":"))

library(jsonlite)
transplant.ch_pd2$near.cosmic.HS <- apply(transplant.ch_pd2[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange")], 1, function(x) {
  p = c(-3:-1,1:9)
  n = c(-9:-1,1:9)
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% ct$GENE.AA.POS
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% ct$CHROM.POS
  if (any(any.in.p)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p], sep="|"))
    )
  } else if (any(any.in.n)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$CHROM.POS[ct$CHROM.POS %in% vector.n], ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n], sep="|"))
    )
  } else {
    return()
  }
})
transplant.ch_pd2$near.cosmic.HS <- sapply(transplant.ch_pd2$near.cosmic.HS, toJSON)
transplant.ch_pd2$near.cosmic.HS <- ifelse(transplant.ch_pd2$near.cosmic.HS == "{}","",transplant.ch_pd2$near.cosmic.HS)

hotspots$near.cosmic.HS <- apply(hotspots[,c("CHROM","POS","SYMBOL_VEP","aa.pos","gene_aachange")], 1, function(x) {
  p = c(-3:-1,1:9)
  n = c(-9:-1,1:9)
  prot = p + as.integer(x["aa.pos"])
  vector.p = paste(x["SYMBOL_VEP"], prot ,sep = ":")
  any.in.p <- vector.p %in% ct$GENE.AA.POS
  nuc = n + as.integer(x["POS"])
  vector.n = paste(x["CHROM"], nuc ,sep = ":")
  any.in.n <- vector.n %in% ct$CHROM.POS
  if (any(any.in.p)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$Gene_HGVSp_VEP[ct$GENE.AA.POS %in% vector.p],ct$GENOMIC_MUTATION_ID[ct$GENE.AA.POS %in% vector.p], sep="|"))
    )
  } else if (any(any.in.n)) {
    # print(paste("Yes",x))
    return(
      c(x[["gene_aachange"]],
        paste(ct$CHROM.POS[ct$CHROM.POS %in% vector.n], ct$GENOMIC_MUTATION_ID[ct$CHROM.POS %in% vector.n], sep="|"))
    )
  } else {
    return()
  }
})
hotspots$near.cosmic.HS <- sapply(hotspots$near.cosmic.HS, toJSON)
hotspots$near.cosmic.HS <- ifelse(hotspots$near.cosmic.HS == "{}","",hotspots$near.cosmic.HS)



#####################################
bf.correction <- 0.05/38793002


transplant.ch_pd2$n.HGVSp<- fillna(transplant.ch_pd2$n.HGVSp, 0)
transplant.ch_pd2$n.HGVSc<- fillna(transplant.ch_pd2$n.HGVSc, 0)
transplant.ch_pd2$isBoltonBickHotspot <- transplant.ch_pd2$n.loci.vep >=5 | transplant.ch_pd2$n.HGVSp | transplant.ch_pd2$n.HGVSc
# colnames(transplant.ch_pd2)[colnames(transplant.ch_pd2)=="isHotSpot"] <- "isBoltonBickHotspot"
transplant.ch_pd2$Vardict_PASS <- ifelse(transplant.ch_pd2$Vardict_FILTER == "PASS" | transplant.ch_pd2$Vardict_FILTER == "BCBIO",1,0)
transplant.ch_pd2$Vardict_PASS <- fillna(transplant.ch_pd2$Vardict_PASS, 0)
transplant.ch_pd2$PON_FISHER <- fillna(transplant.ch_pd2$PON_FISHER, 0)
transplant.ch_pd2$failures2 <- apply(transplant.ch_pd2[,c("Vardict_PASS","Vardict_PON_2AT2_percent","Mutect2_PON_2AT2_percent",
                                                                    "alt_strand_counts_min_2_callers","max.under.0.35","max.over.0.02",
                                                                    "complexity_filters","passed","PON_FISHER","isHotSpot",
                                                                    "Mutect2_CALLER")], 
      1, function (x) {
        reasons <- c("Vardict_PON2:", "Mutect2_PON2:","Passed Both Callers:","SB both Callers:","max.under 35% & not HS:",
                     "max.over 2%:","Complexity filters:","PoN:")
        values <- c(x[["Vardict_PON_2AT2_percent"]], x[["Mutect2_PON_2AT2_percent"]], as.numeric(!x[["passed"]]),
                    as.numeric(!x[["alt_strand_counts_min_2_callers"]]), as.numeric(!x[["max.under.0.35"]]), 
                    as.numeric(!x[["max.over.0.02"]]), as.numeric(!x[["complexity_filters"]]), 
                    as.numeric(!(x[["PON_FISHER"]] <= bf.correction)))
        
        failures = list()
        for (i in 1:length(values)) {
          if (values[i]==1) {
            if (reasons[i]=="max.under 35% & not HS:") {
              if (!x["isHotSpot"]) {
                failures <- append(failures,paste(reasons[i],values[i]))
              }
            } else if (reasons[i]=="Passed Both Callers:") {
              if (!x["Vardict_PASS"]) {
                failures <- append(failures,paste(reasons[i],values[i]))
              }
            }
            else {
              failures <- append(failures,paste(reasons[i],values[i]))
            }
          }
        }
        failures.out <- paste(failures, collapse = ",  ")
        if (!x[["Mutect2_CALLER"]]) {
          failures.out <- paste0(failures.out, ", Vardict Complex (Not called by Mutect)")
        }
        return(failures.out)
      }
)

hotspots$n.HGVSp<- fillna(hotspots$n.HGVSp, 0)
hotspots$n.HGVSc<- fillna(hotspots$n.HGVSc, 0)
hotspots$isBoltonBickHotspot <- hotspots$n.loci.vep >=5 | hotspots$n.HGVSp | hotspots$n.HGVSc
hotspots$Vardict_PASS <- ifelse(hotspots$Vardict_FILTER == "PASS" | hotspots$Vardict_FILTER == "BCBIO",1,0)
hotspots$Vardict_PASS <- fillna(hotspots$Vardict_PASS, 0)
hotspots$PON_FISHER <- fillna(hotspots$PON_FISHER, 0)
ho
hotspots$failures2 <- apply(hotspots[,c("Vardict_PASS","Vardict_PON_2AT2_percent","Mutect2_PON_2AT2_percent",
                                                          "alt_strand_counts_min_2_callers","max.under.0.35","max.over.0.02",
                                                          "complexity_filters","passed","PON_FISHER","isBoltonBickHotspot",
                                                          "Mutect2_CALLER")], 
                                     1, function (x) {
                                       reasons <- c("Vardict_PON2:", "Mutect2_PON2:","Passed Both Callers:","SB both Callers:","max.under 35% & not BoltonBickHS:",
                                                    "max.over 2%:","Complexity filters:","PoN:")
                                       values <- c(x[["Vardict_PON_2AT2_percent"]], x[["Mutect2_PON_2AT2_percent"]], as.numeric(!x[["passed"]]),
                                                   as.numeric(!x[["alt_strand_counts_min_2_callers"]]), as.numeric(!x[["max.under.0.35"]]), 
                                                   as.numeric(!x[["max.over.0.02"]]), as.numeric(!x[["complexity_filters"]]), 
                                                   as.numeric(!(x[["PON_FISHER"]] <= bf.correction)))
                                       
                                       failures = list()
                                       for (i in 1:length(values)) {
                                         if (values[i]==1) {
                                           if (reasons[i]=="max.under 35% & not BoltonBickHS:") {
                                             if (!x["isBoltonBickHotspot"]) {
                                               failures <- append(failures,paste(reasons[i],values[i]))
                                             }
                                           } else if (reasons[i]=="Passed Both Callers:") {
                                             failures <- append(failures,paste(reasons[i],values[i]))
                                           }
                                           else {
                                             failures <- append(failures,paste(reasons[i],values[i]))
                                           }
                                         }
                                       }
                                       failures.out <- paste(failures, collapse = ",  ")
                                       if (!x[["Mutect2_CALLER"]]) {
                                         failures.out <- paste0(failures.out, ", Vardict Complex (Not called by Mutect)")
                                       }
                                       return(failures.out)
                                     }
)



transplant.ch_pd2$failuresGnomAD2 <- apply(transplant.ch_pd2[,c("gnomAD_MAX.Stringent.005","gnomAD_MAX.lessStringent.005")], 
                                          1, function (x) {
                                            reasons <- c("GnomAD Filter More Stringent:", 
                                                         "GnomAD Filter Less Stringent:")
                                            values <- c(as.numeric(!x[["gnomAD_MAX.Stringent.005"]]), as.numeric(!x[["gnomAD_MAX.lessStringent.005"]]))
                                            
                                            failures = list()
                                            for (i in 1:length(values)) {
                                              if (values[i]==1) {
                                                failures <- append(failures,paste(reasons[i],values[i]))
                                              }
                                            }
                                            paste(failures, collapse = ",  ")
                                          })

hotspots$failuresGnomAD2 <- apply(hotspots[,c("gnomAD_MAX.Stringent.005","gnomAD_MAX.lessStringent.005")], 
                                           1, function (x) {
                                             reasons <- c("GnomAD Filter More Stringent:", 
                                                          "GnomAD Filter Less Stringent:")
                                             values <- c(as.numeric(!x[["gnomAD_MAX.Stringent.005"]]), as.numeric(!x[["gnomAD_MAX.lessStringent.005"]]))
                                             
                                             failures = list()
                                             for (i in 1:length(values)) {
                                               if (values[i]==1) {
                                                 failures <- append(failures,paste(reasons[i],values[i]))
                                               }
                                             }
                                             paste(failures, collapse = ",  ")
                                           })

t1 <- readxl::read_excel("~/Bolton/UKBB/ch_pd2_transplant_KB_BW_FINAL_reannotate.xlsx")
t1$key <- with(t1, paste(CHROM,POS,REF,ALT,SAMPLE,sep=":"))
transplant.ch_pd2$samplekey <- paste(transplant.ch_pd2$key, transplant.ch_pd2$SAMPLE, sep = ":")
transplant.ch_pd2.already.have <- transplant.ch_pd2$samplekey %in% t1$key
sum(as.numeric(t1$ch_pd_final[(t1$key %in% transplant.ch_pd2$samplekey)]), na.rm = T)
sum(!transplant.ch_pd2$samplekey %in% t1$key)
new.transplant <- transplant.ch_pd2[!transplant.ch_pd2$samplekey %in% t1$key,]
write.table(new.transplant, "/Users/brian/Bolton/UKBB/new.transplant.passed.mutect.tsv", sep = "\t", quote = F, row.names = F)
write.table(hotspots, "/Users/brian/Bolton/UKBB/hotspots.failed.mutect.tsv", sep = "\t", quote = F, row.names = F)


source("~/Bolton/UKBB/clinvar.R")