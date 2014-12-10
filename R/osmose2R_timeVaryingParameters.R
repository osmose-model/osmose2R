
seasonalCycle = function(A, a, d=1, T, dt) {
  
  if(d!=1 && d!=2) stop("d must be 1 (annual) or 2 (biannual)")
  a = ((a+0.5)%%12)/12 # torus for month parameter
  t = seq(from=1/(2*dt), by=1/dt, len=T*dt) 
  y = A*sin(2*pi*d*(t-a + 1/(4*d)))
  return(y)
  
}

splineTrend = function(x, T, dt) {
  
  time = seq(from=1/(2*dt), by=1/dt, len=T*dt) 
  timeSpline = seq(0, T, len=length(x))
  sf = splinefun(x=timeSpline, y=x)
  y = sf(time)
  return(y)
  
}

splinePeriodic = function(x, T, dt) {
  
  time = seq(from=1/(2*dt), by=1/dt, len=T*dt) 
  timeSpline = seq(0, 1, len=length(x))
  sf = splinefun(x=timeSpline, y=x)
  y = sf(time %% 1)
  return(y)
  
}

replicateParameter = function(x, T, dt) {
  
  n = T*dt/length(x)
  if(n%%1 != 0) stop("T*dt must be a multiple of length(x)")
  out = rep(x, each=n)
  return(out)
  
}

.logisticShift = function(time, L1, L2, shift, amplitude, prob=0.9) {
  if(prob > 1 | prob < 0) stop("prob must be between 0 and 1")
  if(is.na(L2)) return(L1)
  if(L1 == L2) return(L1)
  K = L2 - L1
  k = prob*L2 - prob*L1
  r = 2*log(k/(K-k))/amplitude
  time = time - as.numeric(shift)
  out = L1 + K/(1+exp(-r*time))
  return(out)
}

regimeShift = function(base, shift, T, dt, amplitude=2, log=FALSE) {
  times = seq(from=0, by=1/dt, len=T*dt)
  if(!isTRUE(log)) base = exp(base)
  out = .logisticShift(time=times, L1=base[1], L2=base[2], shift=shift, amplitude=amplitude) 
  if(!isTRUE(log)) out = log(out)
  return(out)
}

timeSeriesParameter = function(base, year=NULL, seasonal=NULL, type, T, dt, log=FALSE,
                               shift=NULL) {
  
  if(!is.null(year) | !is.null(seasonal)) {

    output = switch(type,
                    modelA = modelA(year=year, seasonal=seasonal, T=T, dt=dt),
                    modelB = modelB(year=year, seasonal=seasonal, T=T, dt=dt),
                    stop("model 'type' not recognized.")
    )
     
  } else {
    
    output = rep(0, length=T*dt)
  }
  
  shift = if(is.null(shift)) NA else shift
  base0 = if(is.na(shift)) base[1] else regimeShift(base=base, shift=shift, T=T, dt=dt, log=log)
  
  output = base0 + output
  
  output = if(isTRUE(log)) output else exp(output)
  
  return(output)
  
}

modelA = function(year, seasonal, T, dt) {
  
  if(length(year)!=(T+1) & !is.null(year)) 
    stop("length of 'year' parameters must be T+1.")
  
  if(length(seasonal)!=3 & !is.null(seasonal)) 
    stop("length of 'seasonal' parameters must be 3.")
  
  A = seasonal[1]
  a = seasonal[2]
  d = seasonal[3]
  
  x.year   = if(!is.null(year)) splineTrend(x=year, T=T, dt=dt) else 0
  x.season = if(!is.null(seasonal)) seasonalCycle(A=A, a=a, d=d, T=T, dt=dt) else 0
  
  output = x.year + x.season
  
  return(output)
     
}

modelB = function(year, seasonal, T, dt) {
  
  if(length(year)!=T & !is.null(year)) 
    stop("length of 'year' parameters must be T.")
  if(length(seasonal)!=(12*T) & !is.null(seasonal)) 
    stop("length of 'seasonal' parameters must be 12T.")
  
  x.year   = if(!is.null(year)) replicateParameter(x=year, T=T, dt=dt) else 0
  x.season = if(!is.null(seasonal)) replicateParameter(x=seasonal, T=T, dt=dt) else 0
  
  output = x.year + x.season
  
  return(output)
  
}

randomSeasonalVariates = function(T, min=0.01, max=0.05) {
  mF = NULL
  for(t in 1:T) {
    mf = min + diff(sort(runif(13, max=max))) 
    mf = mf/sum(mf) 
    mF = c(mF, mf)
  }
  return(log(12*mF))
}

calculateMigrationFlux = function(base, par, type, T, dt, proxy=NA) {
  
    output = switch(type,
                    gaussian    = gaussianFlux(base=base, par=par, T=T, dt=dt),
                    linear      = linearFlux(base=base, par=par, T=T, dt=dt),
                    linearproxy = linearProxy(base=base, par=par, T=T, dt=dt, proxy=proxy),
                    logproxy    = logProxy(base=base, par=par, T=T, dt=dt, proxy=proxy),
                    null        = rep(0, len=T*dt),
                    stop("model 'type' not recognized.")
                    )

  return(output)  
   
}



calculatePlanktonAccessibility = function(base, par, type, T, dt, proxy=NA) {
  
  output = switch(type,
                  proxy  = logisticProxy(base=base, par=par, T=T, dt=dt, proxy=proxy),
                  dummy  = dummyProxy(base=base, par=par, T=T, dt=dt, proxy=proxy),
                  null   = rep(base, len=T*dt),
                  stop("model 'type' not recognized.")
  )
  
  return(output)  
  
}
  
linearProxy = function(base, par, T, dt, proxy) {
  
  if(is.null(par)) return(rep(base, dt))
  
  if(any(is.na(proxy))) stop("You must provide a valid time series proxy (NA not allowed).")
  
  proxy = 2*(proxy - min(proxy))/(max(proxy) - min(proxy)) - 1
  
  b = par[1]
  
  output = base*(1 + b*proxy) # b in [-1,1]
                
  return(output)
  
}

logProxy = function(base, par, T, dt, proxy) {
  
  if(is.null(par)) return(rep(base, dt))
  
  if(any(is.na(proxy))) stop("You must provide a valid time series proxy (NA not allowed).")
  
  proxy = (proxy - mean(proxy))/sd(proxy)
  
  b = par[1]
  
  output = base*exp(-b*proxy)
  
  return(output)
  
}

dummyProxy = function(base, par, T, dt, proxy) {
  
  if(is.null(par)) return(rep(base, dt))
  
  if(any(is.na(proxy))) stop("You must provide a valid time series proxy (NA not allowed).")
  
  proxy = as.numeric(as.factor(proxy))

  if(length(unique(proxy))>3) stop("You must provide a proxy with at most 3 levels.")
  
  output = numeric(length(proxy))
  output[proxy==1] = base
  output[proxy==2] = par[1]
  output[proxy==3] = par[2]
  
  return(output)
  
}

logisticProxy = function(base, par, T, dt, proxy) {
  
  if(is.null(par)) return(rep(base, dt))
  
  if(any(is.na(proxy))) stop("You must provide a valid time series proxy (NA not allowed).")
  
  b = par[1]
  
  output = 1/(1 + (1/base - 1)*exp(-b*proxy))
  
  return(output)
  
}

gaussianFlux = function(base, par, T, dt) {
  
  t0 = par[1]
  sd = par[2]
  
  times = seq(from=0, by=1/dt, len=T*dt) + 0.5/dt
  
  output = round(dnorm(times, mean=t0, sd=sd),2)
  output = base*output/sum(output)
  
  
  return(output)
  
}

linearFlux = function(base, par, T, dt) {
  
  output = seq(from=base, to=par[1], len=T*dt)
  
  return(output)
  
}


getLarvalByYear = function(parameters, species, start, end, inter=TRUE) {
  larval = paste("L", seq(from=start, to=end+1), sep=".")
  output = if(isTRUE(inter)) parameters[larval, species] else matrix(0, nrow=end - start + 2, ncol=length(species))
  colnames(output) = species
  return(output)
}

getFishingByYear = function(parameters, species, start, end, inter=TRUE) {
  fishing = paste("F", seq(from=start, to=end), sep=".")
  output = if(isTRUE(inter)) parameters[fishing, species] else matrix(0, nrow=end - start + 1, ncol=length(species))
  colnames(output) = species
  return(output)
}

getFishingByMonth = function(parameters, species, start, end, inter=TRUE) {
  fishing = if(isTRUE(inter)) paste("F", seq(from=start, to=end), sep=".") else paste0("F.", start)
  fishing = paste(rep(fishing, each=12), tolower(month.abb), sep=".")
  output = parameters[fishing, species]
  if(!isTRUE(inter)) output = apply(output, 2, rep, times=end-start+1)
  return(output)
}
