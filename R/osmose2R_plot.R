plot.osmose.abundance = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE,
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



plot.osmose.yieldN = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE, nrep=3,
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


plot.osmose.meanTL = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE,
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

plot.osmose.meanTLCatch = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE,
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

plot.osmose.meanSize = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE,
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

plot.osmose.meanSizeCatch = function(object, start=NULL, conf=0.95, factor=1e-6, replicates=FALSE,
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
