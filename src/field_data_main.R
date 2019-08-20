#main
field_data_main <- function(){
  #'
  #'
  #' @example 
  #' data products download
  #' chemical >>  DP1.10026.001
  #' isotopes >> DP1.10053.001
  #' vegetation structure >> DP1.10098.001
  
  library(tidyverse)
  library(devtools)
  library(neonUtilities)
  library(downloader)
  library(httr)
  library(jsonlite)
  

  #source files
  file.sources = paste("./src/",
    list.files("./src/", pattern="retrieve"), sep="/") 
  sapply(file.sources,source,.GlobalEnv)

  dir.create(file.path("out", "TOS_outputs"))#, showWarnings = FALSE)
  dir.create(file.path("tmp"))#, showWarnings = FALSE)
  
  retrieve_TOS_data(10026)
  retrieve_TOS_data(10053)
  #retrieve_TOS_data(10098)

  #harmonize the three data products to make a single database
  stack_chemical_leaf_products(10026)
  stack_isotopes_leaf_products(10053)

  # get coordinates and position of the vegetation structure trees
  retrieve_vegetation_structure()

  #now connect with field data and position
  retrieve_joint_dataset()
  results = readr::read_csv("./out/TOS_outputs/vegetation_structure_utm.csv") 
  results = get_lat_long(results)
  results <- sf::st_as_sf(results, coords = c("itc_lat", "itc_lon"), crs = 4326)
  a <- results %>% select(taxonID)
  sf::st_write(results, "./out/TOS_outputs/vegetation_structure.shp", delete_layer=TRUE)
  
  #system2("rm -fr ./tmp/*")
}

field_data_main()
