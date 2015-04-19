
# osmose2R for specific versions ------------------------------------------

osmose2R.v3r0 = function(path=NULL, species.names=NULL, ...) {
  
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
    meanTL      = readOsmoseFiles(path=path, type="meanTL"),
    meanTLCatch = readOsmoseFiles(path=path, type="meanTLCatch"),
    predatorPressure = readOsmoseFiles(path=path, type="predatorPressure"),
    predPreyIni = readOsmoseFiles(path=path, type="biomassPredPreyIni")
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
  
  return(output)
  
}

osmose2R.v3r1 = function(path=NULL, species.names=NULL, ...) {
  
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
    dietMatrix  = readOsmoseFiles(path=path, type="dietMatrix"),
    meanTL      = readOsmoseFiles(path=path, type="meanTL"),
    meanTLCatch = readOsmoseFiles(path=path, type="meanTLCatch"),
    predatorPressure = readOsmoseFiles(path=path, type="predatorPressure"),
    predPreyIni = readOsmoseFiles(path=path, type="biomassPredPreyIni"),
    TLDistrib   = readOsmoseFiles(path=path, type="TLDistrib")
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
    version  = "3u1",
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
  
  return(output)
  
}

osmose2R.v3r2 = function (path = NULL, species.names = NULL, ...) {
    
    
    pop = list(biomass = readOsmoseFiles(path = path, type = "biomass"), 
               abundance = readOsmoseFiles(path = path, type = "abundance"), 
               yield = readOsmoseFiles(path = path, type = "yield"), 
               yieldN = readOsmoseFiles(path = path, type = "yieldN"),
               mortality = readOsmoseFiles(path = path, type = "mortalityRate", bySpecies = TRUE))
    
    Trophic = list(meanTL = readOsmoseFiles(path = path, type = "meanTL"), 
                   meanTLCatch = readOsmoseFiles(path = path, type = "meanTLCatch"),
                   biomassByTL = readOsmoseFiles(path = path, type = "biomasDistribByTL"),
                   predatorPressure = readOsmoseFiles(path = path, type = "predatorPressure"), 
                   predPreyIni = readOsmoseFiles(path = path, type = "biomassPredPreyIni"),
                   dietMatrix = readOsmoseFiles(path = path, type = "dietMatrix"))
    
    Size = list(meanSize = readOsmoseFiles(path = path, type = "meanSize"),            
                meanSizeCatch = readOsmoseFiles(path = path, type = "meanSizeCatch"),
                abundanceBySize = readOsmoseFiles(path = path, type = "abundanceDistribBySize"),
                biomassBySize = readOsmoseFiles(path = path, type = "biomasDistribBySize"),
                yieldBySize = readOsmoseFiles(path = path, type = "yieldDistribBySize"),
                yieldNBySize = readOsmoseFiles(path = path, type = "yieldNDistribBySize"),
                meanTLBySize = readOsmoseFiles(path = path, type = "meanTLDistribBySize"),
                mortalityBySize = readOsmoseFiles(path = path, type = "mortalityRateDistribBySize", bySpecies = TRUE),
                dietMatrixBySize = readOsmoseFiles(path = path, type = "dietMatrixbySize", bySpecies = TRUE),
                predatorPressureBySize = readOsmoseFiles(path = path, type = "predatorPressureDistribBySize", bySpecies = TRUE))
    
    
    Age = list(abundanceByAge = readOsmoseFiles(path = path, type = "abundanceDistribByAge"),
               biomassByAge = readOsmoseFiles(path = path, type = "biomasDistribByAge"),
               yieldByAge = readOsmoseFiles(path = path, type = "yieldDistribByAge"),
               yieldNByAge = readOsmoseFiles(path = path, type = "yieldNDistribByAge"),
               meanSizeByAge = readOsmoseFiles(path = path, type = "meanSizeDistribByAge"),
               meanTLByAge = readOsmoseFiles(path = path, type = "meanTLDistribByAge"),
               mortalityByAge = readOsmoseFiles(path = path, type = "mortalityRateDistribByAge", bySpecies = TRUE),
               dietMatrixByAge = readOsmoseFiles(path = path, type = "dietMatrixbyAge", bySpecies = TRUE),
               predatorPressureByAge = readOsmoseFiles(path = path, type = "predatorPressureDistribByAge", bySpecies = TRUE))
    
    
    model = list(version = "3u2",
                 model = .getModelName(path = path), 
                 simus = dim(pop$biomass)[3], 
                 times = as.numeric(row.names(pop$biomass)), 
                 T = nrow(pop$biomass), 
                 start = as.numeric(row.names(pop$biomass))[1], 
                 nsp = ncol(pop$biomass), 
                 lspecies = if (!is.null(species.names)) species.names else colnames(pop$biomass))
    
    output = list(model = model, species = colnames(pop$biomass), 
                  global = pop, trophic = Trophic, size = Size, age = Age)
    
    
    return(output)
  }


