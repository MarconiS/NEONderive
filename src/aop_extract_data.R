aop_extract_data <- function(x, f, itc.f, epsg, token, wd,pybin = "/home/s.marconi/.conda/envs/quetzal3/bin",  buffer = 20){
  library(raster)
  library(rgeos)
  library(rgdal)
  library(readr)
  #256196.6 4108745
  clip.xmin <- max((x$easting) - buffer, as.integer(x$easting/1000)*1000)
  clip.ymin <- max((x$northing) - buffer, as.integer(x$northing/1000)*1000)
  clip.xmax <- min(clip.xmin + 2*buffer, as.integer(x$easting/1000+1)*1000)
  clip.ymax <- min(clip.ymin + 2*buffer, as.integer(x$northing/1000+1)*1000)
  
  #256166.6 4108775
  entryID <- paste(token, x$siteID, x$taxonID,  x$individualID, sep="_")
  #launch python "/Users/sergiomarconi/anaconda3/bin/"
  Sys.setenv(PATH = paste(pybin, Sys.getenv("PATH"),sep=":"))
  system2("python3", args=(sprintf('"%1$s" "%2$s" "%3$s" "%4$f" "%5$f" "%6$f" "%7$f" "%8$d" "%9$s"',
                                   paste(wd, "src/extractCrown.py", sep=""), 
                                   f, entryID, clip.xmin, clip.xmax, clip.ymin, clip.ymax, as.numeric(epsg), wd)))
  #clip<- stack(paste(wd, "/itcTiff/", entryID, '.tif', sep=""))
}

