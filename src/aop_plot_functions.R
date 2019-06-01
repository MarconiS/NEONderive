aop_chm_plot <- function(plots, data, tileID, epsg, paths, bff = 20, cores = 4){
  #' use lidR to clip the lidar data, create a very high resolution chm, and return predicted itcs
  #' given we know the tree tops from NEON TOS
  #' @param plots data.frame.
  #' 
  #' @param tileID data.frame.
  #' @param epsg 
  #' @param paths   
  
  
  dir.create(file.path("out", "AOP"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot", "CHM"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot", "ITCs"))#, showWarnings = FALSE)
  dir.create(file.path("out", "AOP", "plot", "spectra"))#, showWarnings = FALSE)
  
  library(foreach)
  library(doParallel)
  
  # registerDoSEQ()
  # cl <- makeCluster(cores)
  # registerDoParallel(cl)
  # clusterCall(cl, function(x) .libPaths(x), .libPaths())
  # 
  #results <- foreach(ii = 1:nrow(tileID)) %dopar% {
  for(ii in 1:nrow(tileID)){
    
    library(lidR)
    library(sf)
    library(tidyverse)
    library(exactextractr)
    source("./src/utilities.R")
    file.sources = paste("./src/",
                         list.files("./src/", pattern="aop"), sep="/") 
    sapply(file.sources,source,.GlobalEnv)
    
    #for(ii in 1:nrow(tileID)){
    lid_tile <- list.files(paths$pt, pattern = paste(tileID[ii,], collapse = "_"))
    las_tl <- readLAS(paste(paths$pt, lid_tile, sep="/"))
    
    for(jj in 1:nrow(plots)){
      crds <- plots[jj, c("easting", "northing")]
      las <- lasclipRectangle(las_tl, xleft = crds$easting - bff, 
                              crds$northing - bff, 
                              crds$easting + bff, 
                              crds$northing + bff)
      tryCatch({
        las <- lasnormalize(las, tin())
      }, error=function(cond) {
      })
      #laspl = lasnormalize(las, tin())
      tryCatch({
        
        thr <- c(0,2,5,10,15)
        edg <- c(0, 1.5)
        chm <- grid_canopy(las, 0.25, pitfree(thr, edg, subcircle = 0.17))
        chm <- stretch(chm, minq=0.05, maxq=0.95)
        
        raster::writeRaster(chm, filename=paste("./out/AOP/plot/CHM/",
                                                plots[jj,"plotID"], ".tif", sep=""), 
                            format="GTiff", overwrite=TRUE)
        writeLAS(las, paste("./out/AOP/plot/las/",
                            plots[jj,"plotID"], ".las", sep=""))
        #write_csv(sp_check, paste("./out/AOP/spectra/", plots[jj,"plotID"], "test.csv", sep = ""))
      }, error=function(cond) {
        print(plots[jj,"plotID"])
        print(cond)
      })
    }
  }
  #stopCluster(cl)
}

get_itcs <- function(data, site, plots, pt = "~/Documents/Data/plot/"){
  library(lidR)
  library(tidyverse)
  library(raster)
  library(rgdal)  # input/output, projections
  library(sf)
  library(lidR)
  library(tidyverse)
  library(exactextractr)
  source("./src/utilities.R")
  ii = "MLBS_072"
  list_plots <- plots[,"plotID"] %>% unique %>% unlist
  for(ii in list_plots){
    
    tryCatch({
      chm_rgba = raster(paste(pt,"false_col/grayscale_", ii,".tif", sep=""))
      chm = raster(paste(pt,"CHM/", ii,".tif", sep=""))
      las = readLAS(paste(pt,"las/", ii,".las", sep=""))
      color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
      
      rescale <- function(x, max_val) (x-min(as.matrix(x), na.rm=T))/
        (max(as.matrix(x), na.rm=T) - min(as.matrix(x), na.rm=T)) * max_val
      chm_rgba = rescale(x = chm_rgba, max_val = 50)
      
      treetops <- data %>% dplyr::filter(plotID == unlist(ii)) %>%
        dplyr::select("individualID", "UTM_E", "UTM_N", "height", "stemDiameter") %>% #,
        unique %>%
        dplyr::group_by(UTM_E, UTM_N) %>%
        summarise_all(funs(if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
        st_as_sf(coords = c("UTM_E", "UTM_N"), crs = epsg)
      
      #use plot extent plus 3m buffer
      source("./src/utilities.R")
      # settbuff=buffer_bbox(st_bbox(treetops),3)
      # attr(settbuff, "class") = "bbox"
      # attr(st_geometry(treetops), "bbox") = settbuff
      # #crop the chm to the extent of field data collection 
      # chm_itc <- raster::crop(chm, extent(treetops))
      # #use tree data allometry to apply dalponte 2016
      #itcs
      ttops <- tree_detection(las, lmf(2.6,2))
      #algo = dalponte2016(chm_rgba, treetops = ttops)
      silva = silva2016(chm_rgba, treetops = ttops, exclusion = 0.2, max_cr_factor = 0.7)
      
      las  = lastrees(las, silva)
      metric = tree_metrics(las, .stdtreemetrics)
      hulls  = tree_hulls(las, attribute = "treeID")
      hulls@data = dplyr::left_join(hulls@data, metric@data)
      #add itc names and save itcs
      proj4string(hulls) <- crs(chm_rgba)
      hulls <- st_join(st_as_sf(hulls), treetops) %>%
        filter(!is.na(individualID))
      hulls = st_as_sf(hulls)
      st_write(hulls, paste("./out/AOP/plot/ITCs/", ii, ".shp", sep=""), delete_layer=TRUE)
      
      #extract data from hiperspectral
      hps_f = list.files(paste(pt, "/itcTiff", sep=""), pattern = ii)
      hps <- raster::brick(paste(pt, "/itcTiff/", hps_f, sep=""))
      
      #rasterize polygons
      library(fasterize)
      r <- raster(hulls, res = 1)
      r <- fasterize(hulls, r, field = "treeID")
      spdf_2 <- as(r,'SpatialPolygonsDataFrame')
      crs(hps) <- crs(r)
      hps <- crop(hps, r, snap='near')
      extent(hps) <- alignExtent(hps, r)
      extent(r) <- alignExtent(hps, r)
      
      hps <- addLayer(r, hps)
      hps <- data.frame(as.matrix(hps))
      hulls <- hulls %>% dplyr::select(treeID, individualID)
      colnames(hulls)[1] <- "layer"
      hps<- right_join(hulls, hps)
      write_csv(hps,paste(pt,"spectra/reflectance_", ii,".csv", sep=""))
    }, error = function(err) {print(paste(ii, err))})
  }
}
