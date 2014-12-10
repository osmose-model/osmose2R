runOsmose = function(osmose=NULL, java="java", input="input/config.csv", output="output/",
                     options=NULL, log="osmose.log", verbose=NULL, clean=TRUE) {
  
  if(is.null(verbose))  verbose = interactive()
  
  if(is.null(osmose) & interactive())  stop("No default OSMOSE java executable")
  if(is.null(osmose) & !interactive()) osmose = "osmose.jar"

  if(isTRUE(clean)) file.remove(file.path(output, dir(path=output, recursive=TRUE)))
  
  if(is.null(options)) options = ""
  
  run.osmose = paste(java, options, "-jar", osmose, input, output)
  if(!isTRUE(verbose)) run.osmose = paste(run.osmose, ">", log, "2>", log)
  
  system(run.osmose, wait=TRUE)
  
  return(invisible(run.osmose))
  
}


readOsmoseFiles = function(path, type, bySpecies=FALSE, ...) {
  
  xclass = paste("osmose", type, sep=".")
  
  
    allFiles = dir(path=path, recursive=TRUE, include.dirs=FALSE)
    csvFiles = allFiles[grepl(".csv", allFiles)]
    
    if(!isTRUE(bySpecies)) {
      
      type_  = paste0(type, "_")
      files  = csvFiles[grepl(type_, csvFiles)]
      output = .readFilesList(files=files, path=path, type=type, ...)
      
    } else {
      
      type_  = paste0(type, "-")
      files  = csvFiles[grepl(type_, csvFiles)]
      files  = .bySpecies(files=files)
      output = lapply(files, FUN=.readFilesList, path=path, type=type, ...)
      
    }
    
  if(!is.null(output)) class(output) = xclass
  
    return(output)
  
  }


getSizeSpectrum = function(file, sep=",") {
  # use readOsmoseCsv
  sizeSpectrum = read.table(file, sep=sep, dec=".", skip=1,
                            header=TRUE)
  nsp = ncol(sizeSpectrum) - 2
  times = unique(sizeSpectrum$Time)
  lengths = unique(sizeSpectrum$Size)
  
  out = array(dim = c(length(times), length(lengths), nsp))
  
  for(t in seq_along(times)) {
    out[t,,]  = as.matrix(sizeSpectrum[sizeSpectrum$Time==times[t],-(1:2)])
  }
  colnames(out) = lengths
  rownames(out) = round(times,3)
  dimnames(out)[[3]] = paste0("sp.", seq(nsp)-1)
  return(out)
}



getMortality = function(x, stage="adults", type="total") {
  .calcMort = function(x) {
    x = as.data.frame(x)
    x$natural = x$pred + x$starv + x$other + x$out
    x$total = x$natural + x$fishing
    return(x)
  }
  .getZ = function(x, stage, type) {
    x = x[[stage]]
    x = apply(x, 1:2, mean, na.rm=TRUE)
    x = .calcMort(x)
    x = x[, type]
    return(x)
  }
  
  out = sapply(x, .getZ, stage=stage, type=type)
  return(out)
}


getAverageMortality = function(x, stage="adults", freq=12) {
  
  .getZ = function(x, stage) {
    x = x[[stage]]
    x = apply(x, 1:2, mean, na.rm=TRUE)
    x = freq*colMeans(x, na.rm=TRUE)
    return(x)
  }
  
  out = sapply(x, .getZ, stage=stage)
  return(out)
}

getMortalityDeviation = function(x, stage, type, pars=NULL) {
  x     = getMortality(x=x, stage=stage, type=type)
  if(!is.null(pars)) {
    proxy = pars$dt.save*pars$M.proxy/pars$dt    
  } else {
    proxy = colMeans(x)
  }
  out   = t(apply(x, 1, "-", proxy))
  return(out)
}

