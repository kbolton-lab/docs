library(XML)
library(jsonlite)


clinvar <- read.table("/Users/brian/Bolton/UKBB/hotspots.failed.mutect.tsv", sep = "\t", header = T, quote="", comment.char="")
clinvar$clinvar_RCVACC_VEP[is.na(clinvar$clinvar_RCVACC_VEP)] <- ""
UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
clinvar$clinvar_germline <- sapply(clinvar$clinvar_RCVACC_VEP, function(clinvar_RCVACC_VEP) {
  
  ids <- str_split(clinvar_RCVACC_VEP,"&")[[1]]
  if (length(ids) >= 1 & ids!="") {
    total=0
    my_list=list()
    for (id in ids) {
      id_total <- 0
      my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=clinvarset&id=", id)
      doc <- GET(my_url, user_agent(UA))
      print("finished")
      data.xml <- XML::xmlParse(content(doc, "text"))
      values <- xpathSApply(data.xml, "//ClinVarResult-Set/ClinVarSet/ClinVarAssertion/ObservedIn/Sample/Origin", xmlValue)
      if (!is.na(table(values)["germline"])) {
        id_total <- table(values)["germline"]
        total <- total + id_total
        my_list[[id]] <- id_total
      } else {
        message(paste0("no germline for: ",id))
      }
      Sys.sleep(2.3)
    }
    my_list$total <- total
    return(toJSON(my_list[order(-unlist(my_list))]))
  } else {
    return("")
  }
})

write.table(clinvar$clinvar_germline, "~/Bolton/UKBB/hotspots.failed.mutect.clinvar.tsv", sep = "\t", 
            row.names = F, quote = F)

clinvar$dbSNP.MAF[92] <- sapply(clinvar$Existing_variation_VEP[92], function(Existing_variation_VEP) {
  ids <- str_split(Existing_variation_VEP,"&")[[1]]
  ids <- ids[grep("^rs", ids)]
  if (length(ids) >= 1) {
    my_list=list()
    for (id in ids) {
      print(id)
      my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&report=XML&id=", id)
      doc <- GET(my_url, user_agent(UA))
      print("finished")
      data.xml <- XML::xmlParse(content(doc, "text"))
      value <- getNodeSet(doc=data.xml, "//n:DocumentSummary/n:GLOBAL_MAFS/n:MAF[n:STUDY='ALFA']/n:FREQ",xmlValue,
                           namespaces = c(xsi = "https://www.w3.org/2001/XMLSchema-instance",
                                          n = "https://www.ncbi.nlm.nih.gov/SNP/docsum"))[[1]]
      if (length(value) >= 1) {
        my_list[[id]] <- value
      } else {
        message(paste0("no MAF for: ",id))
        my_list[[id]] <- ""
      }
      Sys.sleep(2.3)
    }
    return(toJSON(my_list))
  } else {
    return("")
  }
})


write.table(clinvar$dbSNP.MAF, "~/Bolton/UKBB/hotspots.failed.mutect.dbsnp.tsv", sep = "\t", 
            row.names = F, quote = F)

