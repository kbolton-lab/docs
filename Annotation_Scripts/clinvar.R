clinvar <- read.table("/Users/brian/Bolton/UKBB/new.transplant.passed.mutect.tsv", sep = "\t", header = T, quote="", comment.char="")
rstudioapi::getSourceEditorContext()$path
PATH <- function() {
  rstudioapi::getSourceEditorContext()$path
}
library(XML)
library(jsonlite)
length(unique(clinvar$clinvar_ID))
ids <- unique(clinvar$clinvar_ID[!is.na(clinvar$clinvar_ID)])

h = curl::new_handle()
curl::handle_setopt(h, http_version = 2)
httr::set_config(httr::config(http_version = 0))


clinvar_germline1 <- sapply(ids[1:200], function(id) {
  id <- as.character(id)
  print(paste0(id,"\n"))
  UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
  my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id)
  doc <- GET(my_url, user_agent(UA))
  print("finished")
  data.xml <- XML::xmlParse(content(doc, "text"))
  

  values <- xpathSApply(data.xml, "//eSummaryResult/DocumentSummarySet/DocumentSummary/trait_set/trait/trait_name", xmlValue)
  sorted_val_count <- sort(table(values), decreasing = T)
  
  if (length(sorted_val_count) < 2) {
    df <- setNames(as.data.frame(sorted_val_count),c("count"))
    df$`Interpreted condition` <- row.names(df)
    row.names(df) <- NULL
    df <- df[,2:1]
    clinvar_ids[[id]] <- toJSON(df)
    Sys.sleep(3)
    return(clinvar_ids)
  } else {
    clinvar_ids[[id]] <- toJSON(setNames(as.data.frame(sorted_val_count),c("Interpreted condition","count")))
    Sys.sleep(3)
    return(clinvar_ids)
  }
})

clinvar_germline2 <- sapply(ids[201:400], function(id) {
  id <- as.character(id)
  print(paste0(id,"\n"))
  UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
  my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id)
  doc <- GET(my_url, user_agent(UA))
  print("finished")
  data.xml <- XML::xmlParse(content(doc, "text"))
  
  
  values <- xpathSApply(data.xml, "//eSummaryResult/DocumentSummarySet/DocumentSummary/trait_set/trait/trait_name", xmlValue)
  sorted_val_count <- sort(table(values), decreasing = T)
  
  if (length(sorted_val_count) < 2) {
    df <- setNames(as.data.frame(sorted_val_count),c("count"))
    df$`Interpreted condition` <- row.names(df)
    row.names(df) <- NULL
    df <- df[,2:1]
    clinvar_ids[[id]] <- toJSON(df)
    Sys.sleep(3)
    return(clinvar_ids)
  } else {
    clinvar_ids[[id]] <- toJSON(setNames(as.data.frame(sorted_val_count),c("Interpreted condition","count")))
    Sys.sleep(3)
    return(clinvar_ids)
  }
})

clinvar_germline3 <- sapply(ids[401:600], function(id) {
  id <- as.character(id)
  print(paste0(id,"\n"))
  UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
  my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id)
  doc <- GET(my_url, user_agent(UA))
  print("finished")
  data.xml <- XML::xmlParse(content(doc, "text"))
  
  
  values <- xpathSApply(data.xml, "//eSummaryResult/DocumentSummarySet/DocumentSummary/trait_set/trait/trait_name", xmlValue)
  sorted_val_count <- sort(table(values), decreasing = T)
  
  if (length(sorted_val_count) < 2) {
    df <- setNames(as.data.frame(sorted_val_count),c("count"))
    df$`Interpreted condition` <- row.names(df)
    row.names(df) <- NULL
    df <- df[,2:1]
    clinvar_ids[[id]] <- toJSON(df)
    Sys.sleep(3)
    return(clinvar_ids)
  } else {
    clinvar_ids[[id]] <- toJSON(setNames(as.data.frame(sorted_val_count),c("Interpreted condition","count")))
    Sys.sleep(3)
    return(clinvar_ids)
  }
})

clinvar_germline4 <- sapply(ids[601:800], function(id) {
  id <- as.character(id)
  print(paste0(id,"\n"))
  UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
  my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id)
  doc <- GET(my_url, user_agent(UA))
  print("finished")
  data.xml <- XML::xmlParse(content(doc, "text"))
  
  
  values <- xpathSApply(data.xml, "//eSummaryResult/DocumentSummarySet/DocumentSummary/trait_set/trait/trait_name", xmlValue)
  sorted_val_count <- sort(table(values), decreasing = T)
  
  if (length(sorted_val_count) < 2) {
    df <- setNames(as.data.frame(sorted_val_count),c("count"))
    df$`Interpreted condition` <- row.names(df)
    row.names(df) <- NULL
    df <- df[,2:1]
    clinvar_ids[[id]] <- toJSON(df)
    Sys.sleep(3)
    return(clinvar_ids)
  } else {
    clinvar_ids[[id]] <- toJSON(setNames(as.data.frame(sorted_val_count),c("Interpreted condition","count")))
    Sys.sleep(3)
    return(clinvar_ids)
  }
})

clinvar_germline5 <- sapply(ids[801:1000], function(id) {
  id <- as.character(id)
  print(paste0(id,"\n"))
  UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
  my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id)
  doc <- GET(my_url, user_agent(UA))
  print("finished")
  data.xml <- XML::xmlParse(content(doc, "text"))
  
  
  values <- xpathSApply(data.xml, "//eSummaryResult/DocumentSummarySet/DocumentSummary/trait_set/trait/trait_name", xmlValue)
  sorted_val_count <- sort(table(values), decreasing = T)
  
  if (length(sorted_val_count) < 2) {
    df <- setNames(as.data.frame(sorted_val_count),c("count"))
    df$`Interpreted condition` <- row.names(df)
    row.names(df) <- NULL
    df <- df[,2:1]
    clinvar_ids[[id]] <- toJSON(df)
    Sys.sleep(3)
    return(clinvar_ids)
  } else {
    clinvar_ids[[id]] <- toJSON(setNames(as.data.frame(sorted_val_count),c("Interpreted condition","count")))
    Sys.sleep(3)
    return(clinvar_ids)
  }
})

clinvar_germline6 <- sapply(ids[1001:length(ids)], function(id) {
  id <- as.character(id)
  print(paste0(id,"\n"))
  UA <- "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/41.0.2227.0 Safari/537.36"
  my_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=", id)
  doc <- GET(my_url, user_agent(UA))
  print("finished")
  data.xml <- XML::xmlParse(content(doc, "text"))
  
  
  values <- xpathSApply(data.xml, "//eSummaryResult/DocumentSummarySet/DocumentSummary/trait_set/trait/trait_name", xmlValue)
  sorted_val_count <- sort(table(values), decreasing = T)
  
  if (length(sorted_val_count) < 2) {
    df <- setNames(as.data.frame(sorted_val_count),c("count"))
    df$`Interpreted condition` <- row.names(df)
    row.names(df) <- NULL
    df <- df[,2:1]
    clinvar_ids[[id]] <- toJSON(df)
    Sys.sleep(3)
    return(clinvar_ids)
  } else {
    clinvar_ids[[id]] <- toJSON(setNames(as.data.frame(sorted_val_count),c("Interpreted condition","count")))
    Sys.sleep(3)
    return(clinvar_ids)
  }
})



length(clinvar_germline)
length(c(clinvar_germline,clinvar_germline))
class(c(clinvar_germline,clinvar_germline))
?rbindlist(c(clinvar_germline,clinvar_germline))

length(clinvar_germline1) == 200
length(clinvar_germline2) == 200
length(clinvar_germline3) == 200
length(clinvar_germline4) == 200
length(clinvar_germline5) == 200
length(clinvar_germline6) == 164
final.list <- c(clinvar_germline1, clinvar_germline2, clinvar_germline3,
                clinvar_germline4, clinvar_germline5, clinvar_germline6)
data.frame(final.list)
fromJSON(final.list) %>% as.data.frame
write.table(final.list,"clinvar.final.list")
lapply(final.list, write, "clinvar.final.list.txt", append=TRUE)
names(final.list)
values <- sapply(final.list, toString)
class(final.list[[1]])
toString(final.list[[1]])
t <- data.frame(cbind(names(final.list), sapply(final.list, toString)))
clinvar$clinvar_ID <- ifelse(!is.na(clinvar$clinvar_ID), as.character(clinvar$clinvar_ID), clinvar$clinvar_ID)
clinvar.final <- clinvar %>% left_join(t, by = c("clinvar_ID"="X1"))
colnames(clinvar.final)[30] <- "Interpreted_condition_count"
clinvar.final$Hereditary_cancer_predisposing_syndrome <- ifelse(grepl("Hereditary cancer-predisposing syndrome",clinvar.final$Interpreted_condition_count),1,0)
clinvar.final <- clinvar.final[order(clinvar.final$Hereditary_cancer_predisposing_syndrome, decreasing = T),]
t2 <- clinvar.final[grepl("Hereditary cancer-predisposing syndrome",clinvar.final$Interpreted_condition_count),]
sum(!is.na(clinvar.final$Interpreted_condition_count))
write.table(clinvar.final, "~/Bolton/UKBB/clinvar.final.tsv", sep = "\t", 
            row.names = F, quote = F)

fnlist <- function(x, fil){ z <- deparse(substitute(x))
# cat(z, "\n", file=fil)
nams=names(x) 
for (i in seq_along(x) ){ cat(nams[i], "\t",  x[[i]], "\n", 
                              file=fil, append=TRUE) }
}
fnlist(final.list, "clinvar.final.list2.txt")

