
osmose2R <-
  function(path=NULL, species.names=NULL, ...) {
    if(is.null(path) & interactive()) {
      path = choose.dir(caption="Select OSMOSE outputs folder")
    }
    if(is.null(path)) stop("No path has been provided.")
    
    # General  
    pop = list(
      biomass    = readOsmoseFiles(path=path, type="biomass"),
      abundance  = readOsmoseFiles(path=path, type="abundance"),
      yield      = readOsmoseFiles(path=path, type="yield"),
      catch      = readOsmoseFiles(path=path, type="yieldN"),
      mortality  = readOsmoseFiles(path=path, type="mortalityRate", bySpecies=TRUE)
    )
    
    # Trophic
    Trophic = list(
      #   dietMatrix  = readOsmoseFiles(path=path, type="dietMatrix"),
         meanTL      = readOsmoseFiles(path=path, type="meanTL"),
         meanTLCatch = readOsmoseFiles(path=path, type="meanTLCatch"),
         predatorPressure = readOsmoseFiles(path=path, type="predatorPressure"),
         predPreyIni = readOsmoseFiles(path=path, type="biomassPredPreyIni")
      #   TLDistrib   = readOsmoseFiles(path=path, type="TLDistrib")
    )
    
    # Size indicators
    Size = list(
         meanSize      = readOsmoseFiles(path=path, type="meanSize"),
         meanSizeCatch = readOsmoseFiles(path=path, type="meanSizeCatch"),
         SizeSpectrumN = readOsmoseFiles(path=path, type="SizeSpectrumSpeciesN"),
         SizeSpectrumB = readOsmoseFiles(path=path, type="SizeSpectrumSpeciesB"),
         SizeSpectrumC = readOsmoseFiles(path=path, type="SizeSpectrumSpeciesYield"),
         SizeSpectrumY = readOsmoseFiles(path=path, type="SizeSpectrumSpeciesYieldN")
    )
    
    # Age indicators
    Age = list(
      AgeSpectrumN = readOsmoseFiles(path=path, type="AgeSpectrumSpeciesN"),
      AgeSpectrumB = readOsmoseFiles(path=path, type="AgeSpectrumSpeciesB"),
      AgeSpectrumC = readOsmoseFiles(path=path, type="AgeSpectrumSpeciesYield"),
      AgeSpectrumY = readOsmoseFiles(path=path, type="AgeSpectrumSpeciesYieldN")
    )
    
    model = list(
      version  = "3.0b",
      model    = .getModelName(path=path),
      simus    = dim(pop$biomass)[3],
      times    = as.numeric(row.names(pop$biomass)),
      T        = nrow(pop$biomass),
      start    = as.numeric(row.names(pop$biomass))[1],
      nsp      = ncol(pop$biomass),
      lspecies = if(!is.null(species.names)) species.names else colnames(pop$biomass)
    )
    
    
    output = list(model   = model,
                  species = colnames(pop$biomass),
                  global  = pop,
                  trophic = Trophic,
                  size    = Size,
                  age     = Age
    )
    
    class(output) = c("osmose")
    
    return(output)
  }

plot.osmose = function(x, type, ...) {
  
  switch(type,
         biomass = plot(object=x$global$biomass, ...),
         yield   = plot(object=x$global$yield, ...),
         error("Plot type not defined."))
  
  return(invisible())
}
    

print.osmose <-
  function(x, ...) {
    cat(paste0("OSMOSE v.", x$model$version,"\n"))
    cat("Model", sQuote(x$model$model),"\n")
    cat(x$model$sp, " species modeled (",x$model$simus,
        " simulations):\n", sep="")
    cat(paste(x$species, collapse=", "),".\n", sep="")
  }

print.summary.osmose <-
  function(x, ...) {
    cat(paste0("OSMOSE v.", x$version,"\n"))
    cat("Model", sQuote(x$model),"\n")
    cat(x$sp, "species modeled:\n")
    cat(paste(x$species, collapse=", "),".\n", sep="")
    cat("Main indicators:\n")
    print(x$resumen)
  }

getVar <-
  function(object, var, ...) {
    UseMethod("getVar")
  }

getVar.osmose <-
  function(object, var, type="global", expected=TRUE, ...) {
    out = object[[type]][[var]]
    if(expected) out = apply(out, c(1,2), mean, na.rm=TRUE)
    return(out)
  }


summary.osmose <-
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
