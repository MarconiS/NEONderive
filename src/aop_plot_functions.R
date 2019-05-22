aop_chm_plot <- function(plots, tileID, epsg, paths, bff = 25, cores = 4){
  #' use lidR to clip the lidar data, create a very high resolution chm, and return predicted itcs
  #' given we know the tree tops from NEON TOS
  #' @include lidR
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
  
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  results <- foreach(ii = 1:nrow(tileID)) %dopar% {
    
    library(lidR)
    library(sf)
    library(tidyverse)
    library(exactextractr)
    
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
      treetops <- data %>% dplyr::filter(plotID == unlist(plots[jj,"plotID"])) %>%
        dplyr::select("individualID", "UTM_E", "UTM_N", "height", "stemDiameter") %>% #, 
        #       "maxCrownDiameter", "ninetyCrownDiameter") %>%
        unique %>%
        dplyr::group_by(UTM_E, UTM_N) %>%
        summarise_all(list(~if(is.numeric(.)) mean(., na.rm = TRUE) else first(.))) %>%
        st_as_sf(coords = c("UTM_E", "UTM_N"), crs = epsg)
      
      #use plot extent plus 3m buffer
      settbuff=buffer_bbox(st_bbox(treetops),3)
      attr(settbuff, "class") = "bbox"
      attr(st_geometry(treetops), "bbox") = settbuff
      #crop the chm to the extent of field data collection 
      chm_itc <- raster::crop(chm, extent(treetops))
      #use tree data allometry to apply dalponte 2016
      #itcs
      algo = dalponte2016(chm_itc, treetops = as(treetops, "Spatial"))
      las  = lastrees(las, algo)
      metric = tree_metrics(las, .stdtreemetrics)
      hulls  = tree_hulls(las, attribute = "treeID")
      hulls@data = dplyr::left_join(hulls@data, metric@data)
      
      #add itc names and save itcs
      proj4string(hulls) <- crs(treetops)
      itcs <- st_join(st_as_sf(hulls), treetops) %>%
        filter(!is.na(individualID))
      st_write(itcs, paste("./out/AOP/plot/ITCs/", plots[jj,"plotID"], ".shp"), delete_layer=TRUE)
      
      #extract data from hiperspectral
      hps_f = list.files("./out/AOP/plot/itcTiff", pattern = unlist(plots[jj,"plotID"]))
      hps <- raster::brick(paste("./out/AOP/plot/itcTiff", hps_f, sep="/"))
      # vras <- velox(paste("./out/AOP/plot/itcTiff", hps_f, sep="/"), 
      #               extent = extent(hps), res=c(1,1), crs= crs(chm_itc))
      # #vras$crs= crs(chm_itc)
      #rasterize polygons
      library(fasterize)
      r <- raster(itcs, res = 1)
      r <- fasterize(itcs, r, field = "treeID")
      spdf_2 <- as(r,'SpatialPolygonsDataFrame')
      crs(hps) <- crs(r)
      hps <- crop(hps, r, snap='near')
      extent(hps) <- alignExtent(hps, r)
      extent(r) <- alignExtent(hps, r)
      
      hps <- addLayer(r, hps)
      hps <- data.frame(as.matrix(hps))
      itcs <- itcs %>% dplyr::select(treeID, individualID)
      colnames(itcs)[1] <- "layer"
      hps<- right_join(itcs, hps)
      # #itcs <- itcs %>% arrange(individualID)
      # itcs$ID
      # spectra <- vras$extract(sp= itcs)# , df = T, small = T)
      write_csv(hps, paste("./out/AOP/plot/spectra/", plots[jj,"plotID"], ".csv", sep = ""))
      #write_csv(sp_check, paste("./out/AOP/spectra/", plots[jj,"plotID"], "test.csv", sep = ""))
      }, error=function(cond) {
        warning(plots[jj,"plotID"])
      })
    }
  }
  stopCluster(cl)
}
