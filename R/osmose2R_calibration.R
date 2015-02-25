getOsmoseParameters = function(parameters, constants) {
  # make more general the use of LTL parameters
  
  output = constants
  
  if(is.null(output$scale)) output$scale = 1e6
  
  output = within(output, {
    
    end         = start + T - 1
    plankton    = LTL
    inter       = interannual
    
    M.base      = parameters["M.base", species]
    Linf        = parameters["Linf", species] # von bertalanffy iinf
    k           = parameters["k", species] # von bertalanffy k
    t0          = parameters["t0", species] # von bertalanffy t0
    vb.thr      = parameters["vb.thr", species] # von bertalanffy thr
    b           = parameters["b", species] # allometric power
    c           = parameters["c", species] # 'condition factor'
    B0          = scale*parameters["B0", species]
    starv       = parameters["starv", species] 
    maturity    = parameters["maturity", species] # length at maturity
    
    # accesibilities
    A.plankton    = 10^(-parameters["plankton.a", seq_along(plankton)])
    names(A.plankton) = plankton
    
    if(isTRUE(inter)) {
      
      A.plankton.ts = parameters[c("plankton.b", "plankton.c"), seq_along(plankton)] 
      colnames(A.plankton.ts) = plankton
    } else A.plankton.ts = NULL 
    
    # fishing
    F.base      = log(parameters["F.base", species])
    L50         = parameters["F.L50", species]
    L75         = parameters["F.L75", species]
    F.year      = getFishingByYear(parameters, species, start, end, inter)
    F.month     = getFishingByMonth(parameters, species, start, end, inter)
    
    # larval
    l.base     = grepl(pattern="L.base", x=rownames(parameters))
    l.shift    = grepl(pattern="L.shift", x=rownames(parameters))
    
    L.base      = log(parameters[l.base, species, drop=FALSE])
    L.shift     = if(any(l.shift)) parameters[l.shift, species, drop=FALSE] else NULL
    L.year      = getLarvalByYear(parameters, species, start, end, inter)
    L.seasonal  = parameters[c("L.seasonal1", "L.seasonal2", "L.seasonal3"), species]
    
    # flux inmigration
    flux.mass   = 1e6*parameters["flux.mass", species]
    flux.par    = parameters[c("flux.par1", "flux.par2"), species]
    
    # catchability
    q           = 10^(-parameters[c("q_index1", "q_index2"), species])
    
    access = if(all(c(species, plankton) %in% rownames(parameters))) {
      parameters[c(species, plankton), species]
    } else NULL
    
    
  })
  
  output$flux.length = getMigrationLength(par = output)
  names(output$longevity) = output$species
  names(output$egg.size) = output$species
  names(output$flux.type) = output$species
  names(output$selectivity.type) = output$species
  names(output$selectivity.by) = output$species
  names(output$fishery) = output$species
  # migration parameters
  names(output$migration.by) = output$species
  names(output$migration.type) = output$species
  names(output$M50) = output$species
  names(output$M75) = output$species
  names(output$Mmin) = output$species
  names(output$Mmax) = output$species

  
  return(output)
  
}

setCalibrationParameters = function(par, conf) {
  
  pars.sp = c("population.initialization.biomass", "mortality.natural.rate",
              "mortality.starvation.rate.max", "species.K", "species.lInf",
              "species.vonbertalanffy.threshold.age", "species.t0", 
              "species.length2weight.allometric.power", 
              "species.length2weight.condition.factor",
              "species.maturity.size", "flux.incoming.annual.biomass",
              "flux.incoming.age", "flux.incoming.size")
  
  pars.plk = "plankton.accessibility2fish"
  
  for(isp in par$species) {
    
    n = getSpNo(isp, par$species)
    # write species parameters
    conf[[paste0(pars.sp[1], ".sp", n)]] = as.numeric(par$B0[isp])
    conf[[paste0(pars.sp[2], ".sp", n)]] = as.numeric(par$M.base[isp])
    conf[[paste0(pars.sp[3], ".sp", n)]] = as.numeric(par$starv[isp])
    conf[[paste0(pars.sp[4], ".sp", n)]] = as.numeric(par$k[isp])
    conf[[paste0(pars.sp[5], ".sp", n)]] = as.numeric(par$Linf[isp])
    conf[[paste0(pars.sp[6], ".sp", n)]] = as.numeric(par$vb.thr[isp])
    conf[[paste0(pars.sp[7], ".sp", n)]] = as.numeric(par$t0[isp])
    conf[[paste0(pars.sp[8], ".sp", n)]] = as.numeric(par$b[isp])
    conf[[paste0(pars.sp[9], ".sp", n)]] = as.numeric(par$c[isp])
    conf[[paste0(pars.sp[10], ".sp", n)]] = as.numeric(par$maturity[isp])
    conf[[paste0(pars.sp[11], ".sp", n)]] = setFluxBiomass(par$flux.mass[isp])
    conf[[paste0(pars.sp[12], ".sp", n)]] = as.numeric(par$flux.age[isp])
    conf[[paste0(pars.sp[13], ".sp", n)]] = as.numeric(par$flux.length[isp])
    
  }
  
  
#   for(iplk in seq_along(par$plankton)) {
#     conf[[paste0(pars.plk, ".plk", iplk-1)]] = as.numeric(par$A.plankton[iplk])  
#   }
  
  conf = .setInitialFile(par=par, conf=conf)
  
  return(conf)  
  
}

.setInitialFile = function(par, conf) {
  if(par$initial != "null" & isTRUE(par$inter)) {
    initial.n = sample(seq_len(par$initial)-1, 1)
    initial.file = paste0(par$model, "_initial_", initial.n, ".nc")
    initial.file = file.path(par$initial.path, initial.file)
    conf[["simulation.restart.file"]] = initial.file 
    cat("Using initial file ", initial.file, ".\n", sep="")
  } else {
    conf[["simulation.restart.file"]] = "null"
  }
  return(conf)
}

updateOsmoseConfiguration = function(par, file) {
  
  conf = readOsmoseParameters(file=file)
  conf = setCalibrationParameters(par=par, conf=conf)
  writeOsmoseParameters(conf=conf, file=file)
  
  return(invisible(NULL))
  
}

osmoseGrowth = function(sp, par, n=100, plot=FALSE, add=FALSE, ...) {
  
  par = getGrowthParameters(par=par, sp=sp)
  lifespan = par$longevity
  maturity = par$maturity
  age = seq(from=0, to=1.1*lifespan, len=n)
  l = .osmoseGrowth(age=age, par=par)

  if(isTRUE(plot)) {
    if(!isTRUE(add)) 
      plot(age, l, type="l", lwd=1.5, xlab="Age (years)", ylab="Length (cm)", ...) else
        lines(age, l, lwd=1.5, ...)
    abline(h=maturity, col="red", lty=2)
  }
  return(invisible(data.frame(age=age,length=l)))

}

.osmoseGrowth = function(age, par) {
  
  Linf    = par$Linf
  k       = par$k
  t0      = par$t0
  vb.thr  = par$vb.thr
  egg  = if(!is.null(par$egg.size)) par$egg.size else 0
  
  l    = Linf*(1-exp(-k*(age-t0)))
  lthr = Linf*(1-exp(-k*(vb.thr-t0)))
  l2   = egg + age*(lthr-egg)/vb.thr
  
  l[age<=vb.thr] = l2[age<=vb.thr]
  
  return(l)
  
}

.osmoseGrowthInv = function(length, par) {
  
  Linf    = par$Linf
  k       = par$k
  t0      = par$t0
  vb.thr  = par$vb.thr
  egg     = if(!is.null(par$egg.size)) par$egg.size else 0
  lthr    = Linf*(1-exp(-k*(vb.thr-t0)))
  
  age     = t0 - (1/k)*log(1-length/Linf)
  age2    = vb.thr*(length-egg)/(lthr-egg)
  
  age[length<=lthr] = age2[length<=lthr]
  
  return(age)
  
}

getLength = function(age, par, sp) {
  par = getGrowthParameters(par=par, sp=sp)
  output = .osmoseGrowth(age=age, par=par)
 return(output) 
}

getMigrationLength = function(par) {
  
  length = NULL
  for(isp in par$species) {
    length = c(length, getLength(age=par$flux.age[isp], par=par, sp=isp)) 
  }
  return(length)
}

setFluxBiomass = function(x) {
  x = 0 + !!as.numeric(x)
  return(x)
}


refine.season = function(x, T=24, tiny=1) {
  if(all(is.na(x)) | all(x==0)) {
    nx = numeric(T)
  } else {
    freq = length(x)
    x = c(x[length(x)], x, x[1])
    time = (seq_along(x) - 1.5)/freq
    hfun = splinefun(x=time, y=log(x+tiny))
    ntime = seq(from=0, to=1, len=T+1)[-1] - 0.5/T
    nx    = pmax(exp(hfun(ntime)) - tiny, 0)
    nx    = nx/sum(nx, na.rm=TRUE)    
  }
  return(nx)
}

newTimeSteps = function(x, T=24, dt=2, refine=TRUE) {

if(refine) {
  cat("Going to ", dt*T, " time steps from ", T,".\n", sep="")
  N = length(x)
  newSteps = NULL
  for(n in seq_len(N)) {
    newSteps = c(newSteps, seq(from=dt*x[n], len=dt))
  }
} else {
  if(T%%dt!=0) stop("dt has to be a divisor of T")
  cat("Going to ", round(T/dt,0), " time steps from ", T,".\n", sep="")
  newSteps = unique(floor(x/dt))
} 
  return(newSteps)
}



num3 = function(x) {
  xmin  = min(x, na.rm=TRUE)
  xmax  = max(x, na.rm=TRUE)
  xmean = mean(x, na.rm=TRUE)
  return(c(mean=xmean, min=xmin, max=xmax))
}

.selectivity.edge = function(x, L50) {

  selec = numeric(length(x))
  selec[x >= L50] = 1
  names(selec) = x
  return(selec)
}

.selectivity.log = function(x, L50, L75, tiny=1e-6) {

  s1 = (L50*log(3))/(L75-L50)
  s2 = s1/L50
  selec = 1/(1+exp(s1-(s2*x)))
  selec[selec<tiny] = 0
  names(selec) = x
  return(selec)
  
}

.selectivity.norm = function(x, L50, L75, tiny=1e-6) {
  
  sd = (L75-L50)/qnorm(0.75)
  mean = L50
  selec = dnorm(x, mean=mean, sd=sd)
  selec = selec/max(selec, na.rm=TRUE)
  selec[selec<tiny] = 0
  names(selec) = x
  return(selec)
  
}

.selectivity.lnorm = function(x, L50, L75, tiny=1e-6) {
  
  sd = log(L75/L50)/qnorm(0.75)
  mean = log(L50)
  selec = dlnorm(x, mean=mean, sd=sd)
  selec = selec/max(selec, na.rm=TRUE)
  selec[selec<tiny] = 0
  names(selec) = x
  return(selec)
  
}


.selectivity.equilibrium = function(x, M, tiny=1e-6) {
  selec = exp(-M*x)
  names(selec) = x
  return(selec)
}

calculateSelectivity = function(par, n=1, tiny=1e-6) {

  x = switch(par$by, 
            size   = pretty(c(0, 1.1*par$Linf), n=60*n),
            age    = pretty(c(0, 1.1*par$longevity), n=50*n),
            stop("Invalid selectivity type: by must be 'age' or 'length' ")
            )
  
  par$L75 = max(1.01*par$L50, par$L75)
  
  out = switch(par$type,
               log  = .selectivity.log(x=x, L50=par$L50, L75=par$L75, tiny=tiny),
               norm = .selectivity.norm(x=x, L50=par$L50, L75=par$L75, tiny=tiny),
               lnorm = .selectivity.lnorm(x=x, L50=par$L50, L75=par$L75, tiny=tiny),
               edge = .selectivity.edge(x=x, L50=par$L50),
               stop("Invalid selectivity 'type': currently implemented 'log' and 'edge'. See help.")
               )
    
  return(out)
  
}

getLarvalMortality = function(par, sp, type="modelA", log=FALSE) {
  
  l = timeSeriesParameter(base=par$L.base[, sp], year=par$L.year[, sp], 
                          seasonal=par$L.seasonal[, sp], type=type, T=par$T, 
                          dt=par$dt, shift=par$L.shift[, sp], log=log)
  return(l)
}

getFishingMortality = function(par, sp, type="modelB", log=FALSE) {
  f = timeSeriesParameter(base=par$F.base[sp], year=par$F.year[, sp], 
                          seasonal=par$F.month[, sp], type=type, T=par$T, 
                          dt=par$dt, log=log)
  f = f/par$dt
  
  return(f)
}


getPlanktonAccessibility = function(par, sp, type="proxy", proxy=NA) {
  a = calculatePlanktonAccessibility(base=par$A.plankton[sp], par = par$A.plankton.ts[, sp], 
                                    type=type, T=par$T, dt=par$dt, proxy=proxy)
  return(a)
}

getMigrationFlux = function(par, sp) {
  flux = calculateMigrationFlux(base=par$flux.mass[sp], par=par$flux.par[, sp], 
                                          type=par$flux.type[sp], T=par$T, dt=par$dt)
  return(flux)
}

getSelectivityParameters = function(par, sp) {
  
  output = list()
  output = within(output, {
    
    by        = par$selectivity.by[sp]
    type      = par$selectivity.type[sp]
    longevity = par$longevity[sp]
    Linf      = par$Linf[sp]
    L50       = par$L50[sp]
    L75       = par$L75[sp]
    
  })
  
  return(output)
  
}

getMigrationDistributionParameters = function(par, sp) {
  
  output = list()
  output = within(output, {
    # migration model
    by         = par$migration.by[sp]
    type       = par$migration.type[sp]
    # life history
    longevity  = par$longevity[sp]
    Linf       = par$Linf[sp]
    k          = par$k[sp]
    t0         = par$t0[sp]
    vb.thr     = par$vb.thr[sp]
    egg.size   = par$egg.size[sp]
    mortality  = par$M.proxy # to change for a more standard name
    a          = par$c[sp]
    b          = par$b[sp]
    # migration distribution
    M50       = par$M50[sp]
    M75       = par$M75[sp]
    Mmin      = par$Mmin[sp]
    Mmax      = par$Mmax[sp]
    
  })
  
  return(output)
  
}

getGrowthParameters = function(par, sp) {
  
  output = list()
  output = within(output, {
    
    longevity = par$longevity[sp]
    Linf      = par$Linf[sp]
    k         = par$k[sp]
    t0        = par$t0[sp]
    vb.thr    = par$vb.thr[sp]
    egg.size  = par$egg.size[sp]
    maturity  = par$maturity[sp] # check, only work with sizes! TODO
    
  })
  
  return(output)
  
}

getSelectivity = function(par, sp, n=1, tiny=1e-6) {
  out = calculateSelectivity(par=getSelectivityParameters(par=par, sp=sp), n=n, tiny=tiny)
  return(out)
}


isHarvested = function(par, sp) {
  output = par$fishery[sp]
  return(output)
}

writeFishingFiles = function(par, output="input/fishing", n=1, tiny=1e-6, type="modelB") {
  
  fishing = NULL
  for(isp in par$species) {
    
    if(isHarvested(par=par, sp=isp)) {
      
      f           = getFishingMortality(par=par, sp=isp, type=type) 
      selectivity = getSelectivity(par=par, sp=isp, n=n, tiny=tiny)
      Fs          = f %o% selectivity
      
      fishing = cbind(fishing, f) 
      write.osmose(Fs, file=file.path(output, paste0("F-", isp, ".csv")))
      
    } else {
      
      fishing = cbind(fishing, 0)
      
    }
    
  }
  
  fishing = as.data.frame(fishing)
  colnames(fishing) = par$species
  
  return(invisible(fishing))
}

writeMigrationFluxFiles = function(par, output="input/flux", n=1, tiny=1e-6) {
  
  if(is.null(par$flux.mass)) return(invisible())

  migration = NULL
  
  for(isp in par$species) {
    
    if(par$flux.mass[isp] != 0) {
      
      B.flux     = getMigrationFlux(par=par, sp=isp) 
      flux.dist  = getFluxDistribution(par=par, sp=isp, n=n, tiny=tiny)
      Flux       = B.flux %o% flux.dist
    
      migration  = cbind(migration, B.flux)
      
      write.osmose(Flux, file=file.path(output, paste0("flux-", isp, ".csv")))
      
    } else {
      
      migration  = cbind(migration, rep(0, par$T*par$dt))
      
    }
    
  }
  
  migration = as.data.frame(migration)
  colnames(migration) = par$species
  
  return(invisible(migration))
}

getFluxDistribution = function(par, sp, n=1, tiny=1e-6) {
  out = calculateFluxDistribution(par=getMigrationDistributionParameters(par=par, sp=sp), 
                                  n=n, tiny=tiny)
  return(out)
}

calculateFluxDistribution = function(par, n=1, tiny=1e-6) {

  x = switch(par$by, 
             size   = head(pretty(c(0, par$Linf), n=20*n), -1),
             age    = head(pretty(c(0, par$longevity), n=10*n), -1),
             stop("Invalid selectivity type: by must be 'age' or 'size' ")
  )

  if(par$by=="age")  {
    age  = x
    size = .osmoseGrowth(age=age, par=par)
  }
  
  if(par$by=="size") {
    age  = .osmoseGrowthInv(length=x, par=par) 
    size = x
  }

  size2 = (size + c(size[-1], par$Linf))/2
  age2  = (age + c(age[-1], par$longevity))/2
  weight = par$a*size2^par$b
  
  par$M75 = max(1.01*par$M50, par$M75)
  
  out = switch(par$type,
               log         = .selectivity.log(x=size2, L50=par$M50, L75=par$M75, tiny=tiny),
               norm        = .selectivity.norm(x=size2, L50=par$M50, L75=par$M75, tiny=tiny),
               lnorm       = .selectivity.lnorm(x=size2, L50=par$M50, L75=par$M75, tiny=tiny),
               edge        = .selectivity.edge(x=size2, L50=par$M50),
               uniform     = .selectivity.edge(x=size2, L50=0),
               equilibrium = .selectivity.equilibrium(x=age2, M=par$M),
               null        = rep(1, len=length(x)),
               stop("Invalid migration 'type'. See help.")
  )
  
  names(out) = x
  
  out[x<par$Mmin | x>par$Mmax] = 0
  out = out*weight
  out = out/sum(out, na.rm=TRUE)
  
  return(out)
  
}


writeLarvalMortalityFiles = function(par, output="input/larval", type="modelA") {
  
  larval = NULL
  for(isp in par$species) {
    Ls = getLarvalMortality(par=par, sp=isp, type=type) 
    
    larval = cbind(larval, Ls)
    write.osmose(Ls, file=file.path(output, paste0("larval_mortality-", isp, ".csv")))
    
  }
  
  larval = as.data.frame(larval)
  colnames(larval) = par$species
  
  return(invisible(larval))
}

writePlanktonAccessibilityFiles = function(par, output="input/plankton", type="proxy", proxy=NA) {
  
  plkaccess = NULL
  for(iplk in par$plankton) {
    
    As = getPlanktonAccessibility(par=par, sp=iplk, type=type, proxy=proxy) 
    
    plkaccess = cbind(plkaccess, As)
    write.osmose(As, file=file.path(output, paste0("plankton_accessibility-", iplk, ".csv")))
    
  }
  
  plkaccess = as.data.frame(plkaccess)
  colnames(plkaccess) = par$plankton
  
  return(invisible(plkaccess))
}


writeAccesibilityFile = function(par, file="input/predation/predation-accessibility.csv") {
  if(!is.null(par$access)) write.osmose(x=par$access, file=file)
  return(invisible(!is.null(par$access)))
}


# writeFishingFiles = function(species, fishing, output="input/fishing", tiny=0.01) {
#   species.names = names(species)
#   for(isp in species.names) {
#     pars = species[[isp]]
#     out = fishing[,isp] %o% selectivity(pars, tiny=tiny)
#     write.osmose(out, file=file.path(output, paste0("F-", isp, ".csv")))
#   }  
# }
# 
# 
# createFishingFile = function(sp, species, fishing, output="input/fishing") {
#   pars = species[[sp]]
#   out = fishing[,sp] %o% selectivity(pars)
#   write.osmose(out, file=file.path(output, paste0("F-", sp, ".csv")))
# }



.interpolateParameters = function(input, T, positive=TRUE) {
  time = seq(from=0, to=1, len=length(input))
  newtime = seq(from=0, to=1, len=T+2)
  newtime = newtime[-c(1, length(newtime))]
  sfun = splinefun(x=time, y=input)
  newinput = sfun(x=newtime)
  if(positive) newinput = pmax(newinput, 0)
  return(newinput)
}

interpolateParameters = function(input, T, positive=TRUE) {
  out = apply(input, 2, .interpolateParameters, T=T, positive=positive)
  return(out)
}




write.osmose = function(x, file)   {
  write.table(x=x, file=file, sep=";", col.names=NA, quote=FALSE)
}



# function creating the actual maps and one-set conf file
setOsmoseMaps = function(object, type, interannual=FALSE, start=NULL, end=NULL,
                         sp, lifespan, ages=NULL, frequency=24, toPA=TRUE,
                         prob=TRUE, criteria="MinROCdist",
                         normalize=TRUE, lat=NULL, lon=NULL) {
  
  if(is.null(ages)) ages = seq_len(ceiling(lifespan) + 1) - 1
  
  object = window(object, start=start, end=end) 
  
  nstep = switch(type,
                 climatology = frequency/12,
                 seasonal    = frequency/4,
                 annual      = frequency,
                 frequency)
  
  # take care of different ages
  if(interannual) {
    
    map = switch(type,
                 climatology = object$prediction,
                 seasonal    = getSeasonalMap(object),
                 annual      = climatology(object$prediction, 
                                           object$info$time$year),
                 object$prediction)  
    
    index = switch(type,
                   climatology = "month",
                   seasonal    = "season",
                   annual      = "year",
                   "month")
    
    year.min = getYearMin(object, index)
    year.max = getYearMax(object, index)
    
    nmap = if(is.matrix(map)) 1 else dim(map)[3]
    
    
  } else {
    map = switch(type,
                 climatology = object$climatology,
                 seasonal    = object$season,
                 annual      = object$mean,
                 object$mean)
    
    year.min = 0
    year.max = max(object$info$time$year) - min(object$info$time$year) 
    
    nmap = if(is.matrix(map)) 1 else dim(map)[3]
    
    year.min = rep(year.min, nmap)
    year.max = rep(year.max, nmap)
    
  }
  
  map = removeSector(map, coords=object$info$coords, lat=lat, lon=lon)
  
  thr = kali::getThreshold(x=object, criteria=criteria)
  if(toPA) map = toPA(map, thr, prob=prob)
  if(normalize) map = normalize(map)
  
  
  mapNames = paste(sp, niceSeq(nmap), sep="-")
  # generalize 'steps'
  steps = matrix(seq_len(frequency)-1, ncol=nstep, nrow=nmap, byrow=TRUE)
  
  
  info = list(nmap = nmap, sp=sp, lifespan=lifespan, ages=ages,
              frequency=frequency, steps=steps, year.min=year.min, 
              year.max=year.max)
  
  output = list(map = map, files=mapNames, info=info)
  return(output)
}

createOsmoseMaps = function(..., outdir=NULL, confdir=NULL, write=TRUE, 
                            confile = "maps-parameters.csv") {
  
  .par = function(par, i=NULL, type=NULL) {
    if(!is.null(type)) par = paste(par, type, sep=".")
    sprintf(paste("movement", "map%1d", par, sep="."), i)
  }
  
  maps = list(...)
  if(length(maps)==1 & is.list(maps[[1]])) maps = maps[[1]]
  
  nmap = 0
  conf = list()
  for(i in seq_along(maps)) {
    obj = maps[[i]]
    if(write) writeMaps(obj, outdir=outdir)
    n = obj$info$nmap
    for(j in seq_len(n)) {
      conf[[.par("age.max", nmap)]] = max(obj$info$ages)
      conf[[.par("age.min", nmap)]] = min(obj$info$ages)
      conf[[.par("file", nmap)]]    = file.path("maps", paste0(obj$files[j], ".csv"))
      conf[[.par("season", nmap)]]  = obj$info$steps[j,]
      conf[[.par("species", nmap)]] = obj$info$sp
      conf[[.par("year.min", nmap)]] = obj$info$year.min[j]
      conf[[.par("year.max", nmap)]] = obj$info$year.max[j]
      nmap = nmap + 1
    }
  }
  writeOsmoseParameters(conf=conf, file=file.path(confdir, confile))
  DateStamp(nmap, "maps written.")
  
  return(invisible())
}

niceSeq = function(x, zero=FALSE) {
  ncode = floor(log10(x)) + 1 
  out = sprintf(paste0("%0", ncode, "d"), seq_len(x) - zero)
  return(out)
}

# function acting in object of maps, return invisible the number of maps
# and the conf lines, while printing the maps

writeMaps = function(object, files=NULL, outdir="maps") {
  if(is.null(files)) files = object$files
  if(!file.exists(outdir)) dir.create(outdir)
  maps = rotate(object$map, direction="anticlockwise")
  nmap = if(is.matrix(maps)) 1 else dim(maps)[3]
  for(i in seq_len(nmap)) {
    file = file.path(outdir, paste0(files[i], ".csv"))
    png(filename=file.path(outdir, paste0(files[i], ".png")))
    image.plot(.rotate(maps[,,i]))
    dev.off()
    write.table(maps[,,i], file=file, sep=";", quote=FALSE,
                row.names=FALSE, col.names=FALSE, na="-99")
  }
  return(invisible())
}

# function printing all the maps, concatenating the confs and printing the 
# conf file


writeSpeciesConf = function(data, file) {
  .par = function(par, i, type=NULL) {
    if(!is.null(type)) par = paste(par, type, sep=".")
    paste(par, paste0("sp",i), sep=".")
  }
  output = list()
  for(i in seq(from=0, len=nrow(data))) {
    output[[.par("name", i)]] = as.character(data$name[i+1])
    output[[.par("lifespan", i)]] = data$lifespan[i+1]
    output[[.par("lInf", i)]] = data$Linf[i+1]
    output[[.par("K", i)]] = data$K[i+1]
    output[[.par("t0", i)]] = data$t0[i+1]
    output[[.par("vonbertalanffy.threshold.age", i)]] = data$vb.thr[i+1]
    output[[.par("length2weight.condition.factor", i)]] = data$a[i+1]
    output[[.par("length2weight.allometric.power", i)]] = data$b[i+1]
    output[[.par("maturity", i, data$maturity.type[i+1])]] = data$maturity.thr[i+1]
    output[[.par("relativefecundity", i)]] = data$relative.fecundity[i+1]
    output[[.par("sexratio", i)]] = data$sex.ratio[i+1]
    output[[.par("egg.size", i)]] = data$egg.size[i+1]
    output[[.par("egg.weight", i)]] = data$egg.weight[i+1]
  }
  names(output) = paste("species", names(output), sep=".")
  writeOsmoseParameters(conf=output, file=file)
  return(invisible(output))
}

writeOsmoseParameters = function(conf, file, sep=";") {
  .writeParameter = function(x) {
    out = paste(names(x),paste(x, collapse=sep), sep=sep)
    return(out)
  }
  out = sapply(conf, .writeParameter)
  vars = names(out)
  ind = sort(vars, index.return=TRUE)$ix
  dim(out) = c(length(out), 1)
  out = out[ind,, drop=FALSE]
  rownames(out) = vars[ind]
  write.table(out, file=file, sep="", quote=FALSE, col.names=FALSE)
  return(invisible(out))
}

readOsmoseParameters = function(file, sep=";") {
  
  readConfigurationList(file=file, sep=sep, skip=0)
  
}

getSpNo = function(sp, species) which(species == sp) - 1

writeInputData = function(x, var, input="INPUT") {
  file = file.path(input, "DATA", paste0(var, ".csv"))
  write.csv(x[[var]], file=file)
  return(invisible())
}

getJavaExec = function(par) {
  java.path = par$java
  java = if(file.exists(java.path)) java.path else "java"
  return(java)
}

getInputConf = function(par) {
  input.path = par$input
  input = if(!is.null(input.path)) input.path else "input/config.csv"
  return(input)
}

getOsmoseJar = function(par) {
  jar.file = par$osmose
  jar = if(!is.null(jar.file)) jar.file else "osmose.jar"
  return(jar)
}


getJavaOptions = function(par) {
  options = par$java.options
  if(is.null(options)) options=""
  return(options)
}

