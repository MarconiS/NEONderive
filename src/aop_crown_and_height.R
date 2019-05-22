aop_canopy_height <- function(pt = NULL,wd = NULL, pttrn, cores = 4,segment = F, 
  pybin = "/home/s.marconi/.conda/envs/quetzal3/bin", epsg=NULL, chm_f = NULL, dtm_pt = NULL){
  library(foreach)
  library(doParallel)
  registerDoSEQ()
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function(x) .libPaths(x), .libPaths())
  
  results <- foreach(i = pttrn) %dopar% {
    #for(i in pattern){
    library(raster)
    library(lidR)
    #source(paste(wd, "src/polygonize.R", sep=""))
    tryCatch({
      if(!length(list.files(chm_f, pattern=i))>0){
        f = list.files(pt, pattern = i)
        dtm = raster(paste(dtm_pt, list.files(dtm_pt, pattern = i), sep=""))
        las = readLAS(paste(pt, f, sep=""))
        # normalization
        #dtm <- grid_terrain(las, 1, kriging(k = 10L))
        las <- lasnormalize(las, dtm)
	      thr <- c(0,2,5,10,15)
	      edg <- c(0, 1.5)
	      chm <- grid_canopy(las, 1, pitfree(thr, edg))
        writeRaster(chm, paste(chm_f, i, "_chm.tif",sep=""), format="GTiff", overwrite=TRUE)
      print(i)
      }else{
        print("tile exists")
    	}},error=function(e){message(paste(i,"resulted in", f, ": error!"))})
  }
  stopCluster(cl)
  return(results)
}

