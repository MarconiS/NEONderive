get_epsg_from_utm <- function(utm){
  dictionary <- cbind(32616, 32615, 32617, 32616, 32616, 32616, 32612, 32613, 32617, 32617, 
                      32614, 32618, 32616, 32619, 32617, 32615) 
  colnames(dictionary) <- c("STEI", "CHEQ", "SCBI", "GRSM", "ORNL", "TALL", "MOAB", 
                            "JORN", "OSBS", "MLBS", "KONZ", "HARV", "LENO", "GUAN", "DSNY", "UKFS")
  return(dictionary[colnames(dictionary)==utm])
}

