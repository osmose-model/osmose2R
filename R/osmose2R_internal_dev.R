.readMortalitybyAgeorSizeCsv = function(file, sep=",", skip=1, row.names=1, 
                                  na.strings=c("NA", "NaN"), rm=1, ...) {
  
  x = readLines(file)
  .subSep = function(x) gsub(";", ",", x)
  
  x = lapply(x, .subSep) 
  legend = x[[1]]
  headers = x[2]
  x = x[-c(1:2)]
  x = paste(x, collapse="\n")
  x = read.csv(text=x, header=FALSE, na.strings=na.strings)
  x = x[order(x[,2]),]
  times = unique(x[,1])
  Age = unique(x[,2])
  x = as.matrix(x[, -c(1,2)])
  x = array(as.numeric(x), dim=c(length(times), length(Age), 5))
  rownames(x) = round(times, 2)
  colnames(x) = Age
  dimnames(x)[[3]] = c("pred", "starv", "other", "fishing", "out")
  
  return(x)
}

.read_MortStagebyAgeorSize = function(files, path, ...) {
  
  if(length(files)!=0) {
    
    x = .readMortalitybyAgeorSizeCsv(file.path(path, files[1]), row.names=NULL, ...)
    
    rows = row.names(x)
    cols = dimnames(x)[[3]]
    AgeorSize = colnames(x)
    
    out = array(dim = c(dim(x), length(files)))
    
    out[, , , 1] = x
    
    if(length(files)>1) {
      for(i in seq_along(files[-1])) {
        x = .readMortalitybyAgeorSizeCsv(file.path(path, files[i+1]), row.names=NULL, ...)
        out[, , , i+1]= x
      }
    }
    
    rownames(out) = rows
    dimnames(out)[[3]] = cols
    
    output= lapply(seq(length(AgeorSize)),function(x) out[,x,,])
    names(output)=AgeorSize
    
  } else {
    output = NULL
  }
  
  return(output)
}


.read_2D_ByAgeorSize = function(files, path, ...) {
  
  if(length(files)!=0) {
    
    x = .readOsmoseCsv(file.path(path, files[1]), row.names=NULL)
    
    rows    = unique(x[,1])
    cols    = unique(x[,2])
    slices  = names(x)[-(1:2)]
    
    x = .reshapeOsmoseTableByAgeorSize(x)
    
    out = array(dim = c(dim(x), length(files)))
    
    out[, , , 1] = x
    
    if(length(files)>1) {
      for(i in seq_along(files[-1])) {
        x = .readOsmoseCsv(file.path(path, files[i+1]), row.names=NULL)
        x = .reshapeOsmoseTable(x)
        out[, , , i+1]= x
      }
    }
    
    #out = aperm(out, c(1,2,4,3))
    
    rownames(out) = rows
    colnames(out) = slices
    dimnames(out)[[3]] = cols
    nsp = dim(out)[4]
    
    output=list()
    
    for(i in seq_len(nsp)) {
      y = out[, , , i, drop=FALSE]
      dnames = dimnames(y)[1:3]
      dim(y) = dim(y)[-length(dim(y))]
      dimnames(y) = dnames
      output[[i]] = y
    }
    
    names(output) = names(files)
    
  } else {
    output = NULL
  }
  
  return(output)
}

.reshapeOsmoseTableByAgeorSize = function(x) {
  
  rows    = unique(x[,1])
  cols    = unique(x[,2])
  slices  = names(x)[-(1:2)]
  
  x       = as.matrix(x[,-(1:2)])
  dim(x)  = c(length(cols), length(rows), length(slices))
  x       = aperm(x, c(2,3,1))
  
  dimnames(x) = list(rows, slices, cols)
  
  return(x)
}