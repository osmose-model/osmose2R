calculateSeasonalDistribution = function(x, by=12, log=FALSE, fill=FALSE) {
  
  n0 = length(x)
  if(any(x<0, na.rm=TRUE)) stop("only positive values are allowed in x.")
  if(n0%%by != 0) {
    n = by - n0%%by
    x = c(x, rep(NA, n))
  } 
  x = matrix(as.numeric(x), nrow=by)
  total = colSums(x)
  total[total==0] = 1e-3
  output = as.numeric(t(t(x)/total))
  clim = calculateClimatology(output, by=by, fill=fill)
  if(isTRUE(fill)) output[is.na(output)] = clim[is.na(output)]
  length(output) = n0
  if(isTRUE(log)) output = log(output)
  return(output)
  
}

calculateClimatology = function(x, by=12, fill=FALSE) {
  
  n0 = length(x)
  n = 0
  if(n0%%by != 0) {
    n = by - n0%%by 
    x = c(x, rep(NA, n))
  } 
  x = matrix(as.numeric(x), nrow=by)
  output = apply(x, 1, median, na.rm=TRUE)
  if(isTRUE(fill)) {
    output = rep(output, n + n0)
    length(output) = n0
  }
  return(output)
  
}


aggregatebyFreq = function(x, freq=NULL, FUN=sum) {
  FUN = match.fun(FUN)
  if(!is.null(freq)) {
    ind = rep(seq_len(nrow(x)), each=freq, len = nrow(x)) 
    new2 = as.matrix(aggregate(x, by=list(ind), FUN=FUN)[,-1])
    return(new2)
  } else {
    return(x)
  }
}

aggregatebySize = function(x, by=NULL, FUN=sum) {
  if(!is.null(by)) {
    ind = rep(seq_len(ncol(x)), each=by, len = ncol(x)) 
    new2 = t(as.matrix(aggregate(t(x), by=list(ind), FUN=FUN)[,-1]))
    sizes = as.numeric(colnames(x))
    newsizes = aggregate(sizes, by=list(ind), FUN=min)[,-1]
    marks = newsizes + diff(newsizes)[1]/2
    colnames(new2) = marks
    return(new2)
  } else {
    return(x)
  }
}


formatCaL = function(x, crop=NULL, bin=NULL, freq=NULL, by=NULL, factor=1, start=0, ...) {
  
  sizes = as.numeric(colnames(x))
  times = as.numeric(rownames(x))
  
  if(!is.null(crop)) {
    
    if(!is.null(bin)) crop = crop + c(-bin, bin)/2 # to convert to class marks
    
    s1 = which(sizes < crop[1])
    s2 = which(sizes >= crop[1] & sizes < crop[2])
    s3 = which(sizes >= crop[2])
    
    s1 = apply(x[, s1, drop=FALSE], 1, sum)
    s2 = x[, s2, drop=FALSE]
    s3 = apply(x[, s3, drop=FALSE], 1, sum)
    
    s2[,1]         = s2[,1] + s1
    s2[, ncol(s2)] = s2[, ncol(s2)] + s3
    
    new = s2
    
    } else new = x

  new = aggregatebyFreq(x=new, freq=freq)
  new = aggregatebySize(x=new, by=by)
  
  rownames(new) = if(is.null(freq)) 
    start + times else seq(from=start, by=freq/12, len=nrow(new))      
  
  return(factor*new)

}

formatCaA = function(x, crop=NULL, bin=NULL, freq=NULL, by=NULL, factor=1, start=0, ...) {
  
  sizes = as.numeric(colnames(x))
  times = as.numeric(rownames(x))
  
  if(!is.null(crop)) {
    
    if(!is.null(bin)) crop = crop + c(-bin, bin)/2 # to convert to class marks
    
    s1 = which(sizes < crop[1])
    s2 = which(sizes >= crop[1] & sizes <= crop[2])
    s3 = which(sizes > crop[2])
    
    s1 = apply(x[, s1, drop=FALSE], 1, sum)
    s2 = x[, s2, drop=FALSE]
    s3 = apply(x[, s3, drop=FALSE], 1, sum)
    
    s2[,1]         = s2[,1] + s1
    s2[, ncol(s2)] = s2[, ncol(s2)] + s3
    
    new = s2
    
  } else new = x
  
  new = aggregatebyFreq(x=new, freq=freq)
  new = aggregatebySize(x=new, by=by)
  
  rownames(new) = if(is.null(freq)) 
    start + times else seq(from=start, by=freq/12, len=nrow(new))      
  
  return(factor*new)
  
}


getCatchability = function(sp, sim, obs, var="biomass", nlog=FALSE, byYear=FALSE) {
  bio.sim = sim$biomass[, sp]
  bio.obs = obs[[paste0(sp, ".", var)]]
  q.year = bio.obs/bio.sim
  q = mean(q.year, na.rm=TRUE)
  if(isTRUE(byYear)) q = q.year
  if(isTRUE(nlog)) q = -log10(q)
  return(q)
}

countMapsArea = function(path, pos=2) {
  maps = dir(path=path, pattern="csv")
  spMap = gsub("\\.csv", "", do.call(rbind, strsplit(maps, "-"))[,pos])
  spMap = as.factor(spMap)
  out = do.call(cbind, tapply(maps, INDEX=spMap, FUN=.countOnes, sep=";"))
  return(out)
}