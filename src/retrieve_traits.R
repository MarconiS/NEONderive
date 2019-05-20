
stack_chemical_leaf_products <- function(prd=NULL){
  # now, harmonize the traits data.
  product_name = list.files(paste("./tmp/filesToStack", prd,"/stackedFiles", sep=""))
  product_name = product_name[-which(product_name %in% c("validation.csv", "variables.csv"))]
  final_traits_data = read_csv(paste("./tmp/filesToStack", prd,"/stackedFiles/", product_name[1], sep="")) %>%
    select(-one_of("ovenStartDate", "ovenEndDate", "uid", "analyticalRepNumber", "accuracyQF", "measurementQF", "ligninSampleID",
                   "remarks", "sampleType","laboratoryName", "testMethod", "analyzedBy", "reviewedBy", "analysisDate", 
                   "dataQF","plotID", "collectDate", "chlorophyllSampleID", "sampleCondition", "cnSampleID"))
  
  for(ff in product_name[-1]){
    product_data = read_csv(paste("./tmp/filesToStack", prd,"/stackedFiles/",ff, sep="")) %>%
      select(-one_of("ovenStartDate", "ovenEndDate", "uid", "analyticalRepNumber", "accuracyQF", "measurementQF","ligninSampleID",
                     "remarks", "laboratoryName", "testMethod", "analyzedBy",  "chlorophyllSampleID", "analysisDate", 
                     "sampleCondition", "cnSampleID", "reviewedBy", "dataQF","plotID", "collectDate"))
    colnames(product_data)[which(colnames(product_data) %in% "dryMass")]<-
      paste(gsub('.{4}$', '', ff), "_dryWeight", sep="")
    colnames(product_data)[which(colnames(product_data) %in% "freshMass")]<-
      paste(gsub('.{4}$', '', ff), "_freshWeight", sep="")
    final_traits_data = inner_join(product_data, final_traits_data, by = 
                    c("namedLocation", "domainID", "siteID", "plotType", "sampleID")) 
  }
  #final_traits_data <- final_traits_data[, -grep(".x.x", colnames(final_traits_data))]
  #final_traits_data <- final_traits_data[, -grep(".y", colnames(final_traits_data))]
  write_csv(final_traits_data, './out/TOS_outputs/chemical_data.csv')
}

stack_isotopes_leaf_products <- function(prd=NULL){
  product_name = list.files(paste("./tmp/filesToStack", prd,"/stackedFiles", sep=""))
  product_name = product_name[-which(product_name %in% c("validation.csv", "variables.csv"))]
  final_isotopes_data = read_csv(paste("./tmp/filesToStack", prd,"/stackedFiles/", product_name, sep="")) %>%
    select(c(sampleID, d15N, d13C)) %>%
    write_csv('./out/TOS_outputs/isotopes_data.csv')                          
}
