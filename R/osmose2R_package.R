# osmose2R ----------------------------------------------------------------
# main function, creates 'osmose' class objects

osmose2R =  function(path=NULL, version="v3r2", species.names=NULL, ...) {
  if(is.null(path) & interactive()) {
    path = choose.dir(caption="Select OSMOSE outputs folder")
  }
  if(is.null(path)) stop("No path has been provided.")
  
  output = switch(version, 
                  v3r0 = osmose2R.v3r0(path=path, species.names=species.names, ...),
                  v3r1 = osmose2R.v3r1(path=path, species.names=species.names, ...),
                  v3r2 = osmose2R.v3r2(path=path, species.names=species.names, ...),
                  stop(sprintf("Incorrect osmose version %s", version))
  )
  class(output) = "osmose"
  return(output)
}


# methods for 'osmose' class ----------------------------------------------
 

plot.osmose = function(x, type="biomass", ...) {
  
  switch(type,
         biomass = plot(object=x$global$biomass, ...),
         yield   = plot(object=x$global$yield, ...),
         abundance   = plot(object=x$global$abundance, ...),
         yieldN   = plot(object=x$global$yieldN, ...),
         error("Plot type not defined."))
  
  return(invisible())
}

print.osmose =
  function(x, ...) {
    cat(paste0("OSMOSE v.", x$model$version,"\n"))
    cat("Model", sQuote(x$model$model),"\n")
    cat(x$model$sp, " species modeled (",x$model$simus,
        " simulations):\n", sep="")
    cat(paste(x$species, collapse=", "),".\n", sep="")
  }

print.summary.osmose =
  function(x, ...) {
    cat(paste0("OSMOSE v.", x$version,"\n"))
    cat("Model", sQuote(x$model),"\n")
    cat(x$sp, "species modeled:\n")
    cat(paste(x$species, collapse=", "),".\n", sep="")
    cat("Main indicators:\n")
    print(x$resumen)
  }

getVar =
  function(object, var, ...) {
    UseMethod("getVar")
  }

getVar.osmose =
  function(object, var, type="global", expected=TRUE, ...) {
    out = object[[type]][[var]]
     
    xclass = "list" %in% class(out)
    if(isTRUE(!xclass) & isTRUE(expected))
      out = apply(out, c(1,2), mean, na.rm=TRUE)
    
    return(out)
  }


summary.osmose =
  function(object, ...) {
    
    output = object$model
    output$species = object$species
    biomass = apply(object$global$biomass, 2, mean, na.rm=TRUE)
    yield = apply(object$global$yield, 2, mean, na.rm=TRUE)
    resumen = data.frame(biomass=biomass,
                         yield = yield)
    rownames(resumen) = object$species
    output$resumen = resumen
    
    class(output) = "summary.osmose"
    return(output)
    
  }
