retrieve_field_measures <- function(){
  #' clean vegetation structure data and provide a data frame to be used later in the pipeline
  #'
  folder_f = list.files("./dat/TOS_inputs/NEON_struct-woody-plant")
  tree_data <- NULL
  for(ff in folder_f){
    tree_name = list.files(paste("./dat/TOS_inputs/NEON_struct-woody-plant",ff,sep="/"), pattern = "apparent")
    if(length(tree_name)!=0){
      file_tree_data = read_csv(paste("./dat/TOS_inputs/NEON_struct-woody-plant",ff, tree_name, sep="/"))%>%
        select(c("individualID","tagStatus","growthForm","plantStatus","stemDiameter",
                 "height","baseCrownHeight","breakHeight","breakDiameter","maxCrownDiameter",
                 "ninetyCrownDiameter","canopyPosition","shape", "basalStemDiameter",
                 "basalStemDiameterMsrmntHeight", "maxBaseCrownDiameter", "ninetyBaseCrownDiameter"))

      tree_data <- rbind(tree_data, file_tree_data)
    }
  }
  write_csv(tree_data, './out/TOS_outputs/tree_measurements.csv')
}
