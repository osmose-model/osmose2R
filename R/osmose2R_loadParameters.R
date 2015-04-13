

loadOsmoseParameters = function(File,path=NULL){
    
  require(stringr)
  require(R.utils)
  
  if(is.null(path)){
    
    path = normalizePath(dirname(File)) 
  
  } else {
    
    if(!isAbsolutePath(File)){ 
      
      File = file.path(path,File)
    }
  }
  
  Lines = readLines(File)
  
  Lines_trim = lapply(Lines,str_trim)
  Lines_trim[grep("^[[:punct:]]",Lines_trim)]=NULL
  Lines_trim = Lines_trim[nchar(Lines_trim)!=0]
  Separators = sapply(Lines_trim,.guessSeparator)
  KeySeparator = sapply(Separators,function(x) x=x[1])
    
  Key = mapply(.getKey,Lines_trim,KeySeparator)
  
  Values = mapply(.getValues,Lines_trim,KeySeparator)
  
  names(Values) = tolower(Key)
  
  ValuesDef = Values
  ValuesDef[grep("osmose.configuration",Key)] = NULL
  
  if(length(grep("osmose.configuration",Key))>0){
    
    for (i in grep("osmose.configuration",Key)){
    
      ValuesRec = lapply(Values[[i]],function(x) loadOsmoseParameters(x,path))
    
      ValuesDef = c(ValuesDef,ValuesRec[[1]])
      
    }
  }
  
  if(sum(table(names(ValuesDef))>1)>0) warning(paste0("'",names(which(table(names(ValuesDef))>1)),"' has been described twice, only the first element will be used \n",sep=""))
  return(ValuesDef)
  }
  
  
###########################################################################
###Internal function

.guessSeparator=function(Line){
  
  SEPARATORS = c(equal="=",semicolon=";",coma=",",colon=":",tab="\t")
    
  separator = SEPARATORS[lapply((str_split(Line,SEPARATORS)),length)>1]
  
  return(separator)
  
}

.getKey = function(Line,KeySeparator){
  
  Key = str_split(Line,KeySeparator)[[1]][1]
  return(str_trim(Key))
}

  
  .getValues = function(Line,KeySeparator){
        
    Values = str_sub(Line,gregexpr(KeySeparator,Line)[[1]][1]+1,nchar(Line))
  
    ValueSeparator = .guessSeparator(Values)
  
    if(length(ValueSeparator)==0){ValueSeparator="NA"}else{ValueSeparator=ValueSeparator}

    Value = str_trim(str_split(Values,ValueSeparator)[[1]])

    Value = Value[nchar(Value)!=0]
    
    return(list(Value))
  }

getParameters<-function(OsmoseParameters,par,what=NULL){
  
  if(!is.null(what)){
    
    if(what=="numeric") as.numeric(OsmoseParameters[[tolower(par)]]) else OsmoseParameters[[tolower(par)]]

  }else{
  
  OsmoseParameters[[tolower(par)]]
  
  }
}