.bySpecies = function(files, sep=c("_", "-")) {
    out = NULL
    if(length(files)>0) {
      sp  = sapply(sapply(files, FUN=.strsplit2v, sep[1],
                          USE.NAMES=FALSE)[2,], FUN=.strsplit2v, sep[2],
                   USE.NAMES=FALSE)[2,]
      out = as.list(tapply(files, INDEX=sp, FUN=identity))
    }
    # change names for all species
    return(out)
  }

.strsplit2v = function(...) {
    out = matrix(unlist(strsplit(...)), ncol=1)
    names(out) = NULL
    return(out)
  }

.getModelName <-
  function(path) {
    strsplit(dir(path=path, pattern="_biomass_")[1],"_")[[1]][1]
  }

.readOsmoseCsv = function(file, sep=",", skip=1, row.names=1, 
                          na.strings=c("NA", "NaN"), rm=1, ...) {
  out = read.csv(file=file, sep=sep, skip=skip, 
                 row.names=row.names, na.strings=na.strings, ...)
#   mat = as.matrix(out)
#   out[is.na(mat) & !is.nan(mat)] = Inf
#   out[is.na(mat) & is.nan(mat)] = NA
}


.readMortalityCsv = function(file, sep=",", skip=1, row.names=1, 
                             na.strings=c("NA", "NaN"), rm=1, ...) {
  
  x = readLines(file)
  .subSep = function(x) gsub(";", ",", x)

  x = lapply(x, .subSep)
  legend = x[[1]]
  headers = x[2:3]
  x = x[-c(1:3)]
  x = paste(x, collapse="\n")
  x = read.csv(text=x, header=FALSE, na.strings=na.strings)
  times = x[,1]
  x = as.matrix(x[, -c(1,17)])
  x = array(as.numeric(x), dim=c(nrow(x), 3, 5))
  rownames(x) = round(times, 2)
  dimnames(x)[[3]] = c("pred", "starv", "other", "fishing", "out")
  
  return(x)
}

.errorReadingOutputs = function(e) {
  warning(e)
  return(invisible(NULL))
}

.warningReadingOutputs = function(type) {
  e = sprintf("File type '%s' is not recognized by Osmose2R", type)
  warning(e)
  return(invisible(NULL))
}


.readFilesList = function(files, path, type, ...) {
  output = tryCatch(
    switch(type,
           abundance       =  .read_1D(files=files, path=path, ...),
           biomass         =  .read_1D(files=files, path=path, ...),
           yield           =  .read_1D(files=files, path=path, ...),
           yieldN          =  .read_1D(files=files, path=path, ...),
           meanTL          =  .read_1D(files=files, path=path, ...),
           meanTLCatch     =  .read_1D(files=files, path=path, ...),
           meanSize        =  .read_1D(files=files, path=path, ...),
           meanSizeCatch   =  .read_1D(files=files, path=path, ...),
           biomassPredPreyIni =  .read_1D(files=files, path=path, ...),
           predatorPressure          = .read_2D(files=files, path=path, ...),
           dietMatrix                = .read_2D(files=files, path=path, ...),
           AgeSpectrumSpeciesB       = .read_2D(files=files, path=path, ...),
           AgeSpectrumSpeciesN       = .read_2D(files=files, path=path, ...),
           AgeSpectrumSpeciesYield   = .read_2D(files=files, path=path, ...),
           AgeSpectrumSpeciesYieldN  = .read_2D(files=files, path=path, ...),
           SizeSpectrumSpeciesB      = .read_2D(files=files, path=path, ...),
           SizeSpectrumSpeciesN      = .read_2D(files=files, path=path, ...),
           SizeSpectrumSpeciesYield  = .read_2D(files=files, path=path, ...),
           SizeSpectrumSpeciesYieldN = .read_2D(files=files, path=path, ...),
           mortalityRate             = .read_MortStage(files=files, path=path, ...),
           mortalityRateDistribByAge = .read_MortStagebyAgeorSize(files=files, path=path, ...),
           mortalityRateDistribBySize = .read_MortStagebyAgeorSize(files=files, path=path, ...),
           # osmose 3r1
           abundanceDistribBySize         = .read_2D(files=files, path=path, ...),
           biomasDistribBySize            = .read_2D(files=files, path=path, ...),
           naturalMortalityDistribBySize  = .read_2D(files=files, path=path, ...),
           naturalMortalityNDistribBySize = .read_2D(files=files, path=path, ...),
           yieldDistribBySize             = .read_2D(files=files, path=path, ...),
           yieldNDistribBySize            = .read_2D(files=files, path=path, ...),
           abundanceDistribByAge          = .read_2D(files=files, path=path, ...),
           biomasDistribByAge             = .read_2D(files=files, path=path, ...),
           meanSizeDistribByAge           = .read_2D(files=files, path=path, ...),
           naturalMortalityDistribByAge   = .read_2D(files=files, path=path, ...),
           naturalMortalityNDistribByAge  = .read_2D(files=files, path=path, ...),
           yieldDistribByAge              = .read_2D(files=files, path=path, ...),
           yieldNDistribByAge             = .read_2D(files=files, path=path, ...),
           biomasDistribByTL              = .read_2D(files=files, path=path, ...),
           dietMatrixbyAge                = .read_2D_ByAgeorSize(files=files, path=path, ...),
           dietMatrixbySize               = .read_2D_ByAgeorSize(files=files, path=path, ...),
           meanTLDistribByAge             = .read_2D(files=files, path=path, ...),
           meanTLDistribBySize            = .read_2D(files=files, path=path, ...),
           predatorPressureDistribByAge   = .read_2D_ByAgeorSize(files=files, path=path, ...),
           predatorPressureDistribBySize  = .read_2D_ByAgeorSize(files=files, path=path, ...),
           .warningReadingOutputs(type)), 
    error = .errorReadingOutputs)
  
  return(output)
}

.read_1D = function(files, path, ...) {
    # TO_DO: change for the unified approach! species as list
    if(length(files)!=0) {
      x = .readOsmoseCsv(file.path(path, files[1]), ...)
      species = names(x)
      times   = rownames(x)
      
      output = array(dim=c(dim(x),length(files)))
      output[,,1] = as.matrix(x)
      if(length(files)>1) {
        for(i in seq_along(files[-1])) {
          x = .readOsmoseCsv(file.path(path, files[i+1]), ...)
          output[,,i+1]= as.matrix(x)
        }
      }
      rownames(output) = times
      colnames(output) = species
    } else {
      output = NULL
    }
    
    return(output)
  }

.read_2D = function(files, path, ...) {
  
  if(length(files)!=0) {

    x = .readOsmoseCsv(file.path(path, files[1]), row.names=NULL, ...)
    
    rows    = unique(x[,1])
    cols    = unique(x[,2])
    slices  = names(x)[-(1:2)]
    
    x = .reshapeOsmoseTable(x)

    out = array(dim = c(dim(x), length(files)))
    
    out[, , , 1] = x
    
    if(length(files)>1) {
      for(i in seq_along(files[-1])) {
        x = .readOsmoseCsv(file.path(path, files[i+1]), row.names=NULL, ...)
        x = .reshapeOsmoseTable(x)
        out[, , , i+1]= x
      }
    }
    
    out = aperm(out, c(1,2,4,3))
    
    rownames(out) = rows
    colnames(out) = cols
    
    nsp = dim(out)[4]
    
    output=list()
    
    for(i in seq_len(nsp)) {
      y = out[, , , i, drop=FALSE]
      dnames = dimnames(y)[1:3]
      dim(y) = dim(y)[-length(dim(y))]
      dimnames(y) = dnames
      output[[i]] = y
    }
    
    names(output) = slices

  } else {
    output = NULL
  }
  
  return(output)
}

.read_MortStage = function(files, path, ...) {
  
  if(length(files)!=0) {
    
    x = .readMortalityCsv(file.path(path, files[1]), row.names=NULL, ...)
    
    rows = row.names(x)
    cols = c("pred", "starv", "other", "fishing", "out")
    
    out = array(dim = c(dim(x), length(files)))
    
    out[, , , 1] = x
    
    if(length(files)>1) {
      for(i in seq_along(files[-1])) {
        x = .readMortalityCsv(file.path(path, files[i+1]), row.names=NULL, ...)
        out[, , , i+1]= x
      }
    }
    
    rownames(out) = rows
    dimnames(out)[[3]] = cols
    
    output=list()
    
    output$eggs      = out[, 1, ,]
    output$juveniles = out[, 2, ,]
    output$adults    = out[, 3, ,]
    
  } else {
    output = NULL
  }
  
  return(output)
}


.reshapeOsmoseTable = function(x) {
  
  rows    = unique(x[,1])
  cols    = unique(x[,2])
  slices  = names(x)[-(1:2)]
  
  x       = as.matrix(x[,-(1:2)])
  dim(x)  = c(length(cols), length(rows), length(slices))
  x       = aperm(x, c(2,1,3))
  
  dimnames(x) = list(rows, cols, slices)
  
  return(x)
}


.rewriteOutputs = function(path) {
  dir.create(file.path(path, "osmose2R"))
  # not finished
}

.countOnes = function(files, ...) {
  
  out = numeric(length(files))
  
  for(i in seq_along(files)) {
    
    x = read.csv(files[i], header=FALSE, ...)
    out[i] = sum(x>0, na.rm=TRUE)
    
  }
  
  out = c(min=min(out), mean=mean(out), 
          median=median(out), max=max(out))
  return(out)
  
}

.plotCI = function(x, y, prob, col, replicates=FALSE, nrep=3, lwd=2.5, alpha=0.1) {
  
  if(dim(x)[3]==1) {
    lines(x=y, y=apply(x, 1, mean, na.rm=TRUE), col=col)
    return(invisible(NULL))
  }
  
  x.inf = apply(x, 1, quantile, prob=prob/2)
  x.sup = apply(x, 1, quantile, prob=1-prob/2)
  x.50  = apply(x, 1, median)
  x.pol = c(y, rev(y), y[1])
  y.pol = c(x.inf, rev(x.sup), x.inf[1])
  polygon(x.pol, y.pol, col=makeTransparent(col=col, alpha=alpha), border=NA)
  if(isTRUE(replicates)) {
    nrep = max(min(nrep, dim(x)[3]),2)
    matplot(y, x[,,seq_len(nrep)], add=TRUE, type="l", lty=1, 
            col=makeTransparent(col=col, alpha=(alpha + 2)/3))
  }
  lines(y, x.50, col=col, lwd=lwd)
  return(invisible(NULL))
  
}

plot.osmose.biomass = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE,
                               freq=12, alpha=0.5, col="black", xlim=NULL, ylim=NULL, nrep=3,
                               aggregate=FALSE, ...) {
  
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  
  par(oma=c(1,1,1,1), mar=c(3,4,1,1))
  par(mfrow=getmfrow(ncol(object)))
  
  if(isTRUE(aggregate)) {
    .plotAverageBiomass(object, col=col, ...)
    return(invisible())
  }
  
  species = colnames(object)
  start   = if(is.null(start)) as.numeric(rownames(object)[1]) else start
  
  for(sp in species) {
    .plotBiomass(x=object, sp=sp, start=start, conf=conf, factor=factor, 
                 replicates=replicates, nrep=nrep, freq=freq, col=col, alpha=alpha, 
                 xlim=xlim, ylim=xlim) 
      
  }
    
  return(invisible())
}

.removeZeros = function(object) {
  remove = apply(object, 2, function(x) all(x==0))
  object = object[ , !remove, ]
  return(object)
}

plot.osmose.yield = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE, nrep=3,
                               freq=12, alpha=0.5, col="black", xlim=NULL, ylim=NULL, 
                               aggregate=FALSE, zeros=TRUE, ...) {
  

  if(!isTRUE(zeros)) object = .removeZeros(object)
    
  opar = par(no.readonly = TRUE)
  on.exit(par(opar))
  
  if(isTRUE(aggregate)) {
    .plotAverageYield(object, col=col, ...)
    return(invisible())
  }
  
  par(oma=c(1,1,1,1), mar=c(3,4,1,1))
  par(mfrow=getmfrow(ncol(object)))
  
  species = colnames(object)
  start   = if(is.null(start)) as.numeric(rownames(object)[1]) else start
  
  for(sp in species) {
    .plotBiomass(x=object, sp=sp, start=start, conf=conf, factor=factor, freq=freq, nrep=nrep,
                 col=col, alpha=alpha, xlim=xlim, ylim=xlim, replicates=replicates) 
    
  }
  
  return(invisible())
}


.plotAverageBiomass = function(x, col="grey", factor=1e-6, border=NA, ...) {
  
  x = factor*apply(x, 2, mean, na.rm=TRUE)
  barplot(x, border=border, col=col, ...)
  return(invisible())
  i
}

.plotAverageYield = function(x, col="grey", factor=1e-6, ...) {
  
  x[x==0] = NA
  x = as.data.frame(factor*apply(x, 1:2, mean, na.rm=TRUE))
  boxplot(x, col=col, ...)
  return(invisible())
}

.plotBiomass = function(x, sp, start, conf=0.95, factor=1e-6, freq=12, replicates=FALSE, nrep=3,
                                      col="black", alpha=0.5, xlim=NULL, ylim=NULL) {
  
  prob = 1 - conf
  
  x = factor*x[, sp, , drop=FALSE]
  times = seq(from=start + 0.5/freq, by=1/freq, len=nrow(x))
  xlim = if(is.null(xlim)) range(times)
  ylim = if(is.null(ylim)) c(0.75, 1.25)*range(x)
  
  plot.new()
  plot.window(xlim=xlim, ylim=ylim)
  
  .plotCI(x=x, y=times, prob=prob, col=col, alpha=alpha, replicates=replicates, nrep=nrep)
  mtext(toupper(sp), 3, line=-1.5, adj=0.05, cex=0.75)
  axis(1)
  axis(2, las=2)
  box()
  
  return(invisible())
}


.getmfrow = function(n) {
  m1 = floor(sqrt(n))
  m2 = ceiling(n/m1)
  out = rev(sort(c(m1, m2)))
  return(out)
}

getmfrow = function(n) .getmfrow(n=n)

makeTransparent = function(..., alpha=0.5) {
  
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
}
