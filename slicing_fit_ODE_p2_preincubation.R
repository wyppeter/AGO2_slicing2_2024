# Slicing assay fitting and graphing
# For preincubation datasets
# Requires helper file ODElib.R
# Bartel Lab
# Peter Y. Wang 2024

# LIBRARIES ----
# tidyverse
library(tidyverse)
# graphing
library(ggnewscale)
library(viridis)
library(ggpubr)
# NLS (guesser)
library(nlstools)
library(minpack.lm)
# ODEs
library(deSolve)
library(FME)

# SELECT MODEL MODE ----
modelmode = "CTS"

# SAVE OR NOT? ----
SAVE_PLOT = T
SAVE_TABLE = T

# DISPLAY SETUP ----
options(scipen = 5)
formatCoefs = function(x, dp = 3){ as.character(format(round(as.numeric(x), dp), nsmall = dp)) }

# CONSTANTS ----
difflim = log(6/60)  # given known miR-7 kon; nM-1 s-1 (10-10 M-1 s-1); diffusion limit; 6 in nM-1 min-1, 0.1 in nM-1 s-1
inf.conc = 1e6  # nM; a large number to graph limit behavior (no kon effects)
Rconc = 0.02  # nM  # only for plotting here
EXCHANGE_DIL = 1/2.75

# Set ggplot2 theme ----
theme0 = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black", linewidth = 0.75),
  axis.ticks = element_line(colour = "black", linewidth = 0.75),
  axis.text = element_text(color = "black", size = 12),
  axis.title = element_text(color = "black", size = 12),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12, colour = "black"),
  strip.text.y = element_text(size = 12, colour = "black", angle = 0)
)

# READ DATA ----
# [A] >= [R]*4
df.raw = read.csv("slicing_data_preincubation.csv",
  stringsAsFactors = T) %>%
  mutate(rxn = factor(paste0(miR, "_", preincubation)))

# PROCESS AND COMBINE DATA ----
df = df.raw %>%
  # Concentration dilution correction due to buffer exchange
  mutate(conc = conc * EXCHANGE_DIL) %>%
  ## MISC CLEANUP
  mutate(time = as.integer(round(time*60))) %>%  # rounding to the time resolution of ODE, and convert to seconds
  # filter(conc != 0, time != 0) %>%  # discard control points
  mutate(conc_d = factor(conc)) %>%  # factor-ize some variables
  arrange(rxn, conc)

# Note down the guide names used here
rxns = as.character(unname(levels(df$rxn)))

# Simplify df for fitting
df.dt = df %>%
  select(rxn, miR, preincubation, conc, Rconc, time, fraccleaved)

#############################################################
# ODE SETUP ----

# Asymptotic floating point protection
eps = 1e-25 # 1 molecule in 10 uL is 1/(Avo*10^-5) M ~ 10^-19 M = 10^-10 nM.
# Event triggered if state variable <= eps
rootfun = function (t, y, pars) y - eps
# Sets state variable to zero
eventfun = function (t, y, pars) {
  y[y <= eps] = 0
  return(y)
  }

# Generate time points for ODE
# First non-zero time point is one second
Tend = 10 * 60 * 60  # 10 hours
Tthres = c(10, 30, 60) * 60  # 10 min, 0.5 hr, then 1 hr
allTime = c(
  seq(0,         Tthres[1]-1, 1),   # 1 second
  seq(Tthres[1], Tthres[2]-1, 5),   # 5 seconds
  seq(Tthres[2], Tthres[3]-1, 60),  # 1 minute
  seq(Tthres[3], Tend,        1*60) # 1 minute
)
timeGen = function(Tmax) {
  if(Tmax>Tend){stop("time exceeds expected upper bound")}
  return(allTime[allTime <= Tmax])
}

############################################################
# ODE FUNCTIONS ----

initConc = function(d){
  # Generate conc setup as a named vector
  c(conc = d$conc[1], Rconc = d$Rconc[1])
}

ODEcalc = function(state.setup, times, parameters, params.static){
  # Given initial concs, time points, and current parameters,
  # provide the entire simulation state trajectory over time points

  # Define A_T and R_T
  A_T = state.setup["conc"]  %>% unname()
  R_T = state.setup["Rconc"] %>% unname()

  # Include static params
  paramsplus = c(parameters, params.static)

  # Generate the initial conditions
  state = ODEstateinit(modelmode, A_T, R_T, paramsplus)

  # Apply ODE
  out.model = ode(
    y = state,
    times = times,
    func = ODElib(modelmode),
    parms = paramsplus,
    method = "lsode",
    mf = 22,
    rtol = 1e-6,
    atol = 1e-10,
    # these two lines provide a floating point protection at asymptotes
    rootfun = rootfun,
    events = list(func = eventfun, root = T)
  ) %>%
    as.data.frame() %>%
    mutate(conc = A_T, Rconc = R_T,
           fraccleaved.ode = P/R_T)
  out.model # return value
}

ODEcompare = function(df.this, params, params.static){
  # Given the data of a miR at a particular conc set, calculate the ODE outcome,
  # and give a vector of residuals

  # Generate time sequence
  Tmax = max(df.this$time)
  t = timeGen(Tmax)

  # Calculate ODE simulation and compare
  overlaid = ODEcalc(
    state.setup = initConc(df.this),
    times = t,
    parameters = params,
    params.static = params.static
  ) %>%
    left_join(
      df.this, .,
      by = c("conc","Rconc","time"),
      keep = F
    ) %>%
    ungroup() %>%
    mutate(residuals = fraccleaved.ode - fraccleaved)

  if(any(is.na(overlaid$fraccleaved.ode))){
    print(overlaid)
    stop("Crit. error at ODEcompare!")
  }

  overlaid$residuals # return value
}

ODEcost = function(df.data, params, params.static) {
  # Cost function wrapper
  # For every conc-Rconc pair,
  # given parameters, generate simulations, and calculate how much it deviates from the data
  # Returns a compiled vector of residuals (per modFit requirements)

  df.cost = df.data %>%
    group_by(conc, Rconc) %>%
    do(res = ODEcompare(
      df.this = .,
      params = params,
      params.static = params.static
    ))

  unlist(df.cost$res) # return value
}

#####################################################################

# Generating the parameter bounds ----
# Fa.lo.lim = 0.60
# Fa.lo.lim.guess = 0.75
Fa.lo.lim = 0.2
Fa.lo.lim.guess = 0.3

Fa.default = 0.85
kph2.hi.lim =       log(0.2   /60) # s-1  # shouldn't be higher than kcat
kph2.hi.lim.guess = log(0.01  /60) # s-1  # greater constraint for guessing
kph2.lo.lim =       log(0.0001/60) # s-1  # detection limit
kph2.guess = log(0.0002/60) # s-1  # same guess for all. prevents starting at 0 which traps fit to Inf in ODE
uLim = 1 - 1e-2  # defense against limit placement of init ## changed
lLim = 1e-6      # defense against limit placement of init
kcat.init = 1/60  # lit lower end value
kcat.uLim = 60/60  # set ceiling to avoid runaway fits for NLS guessing

if(modelmode == "MAX"){
  #              kon            kcat         kph2                Fa
  lowerBound = c(kon = -Inf,    kcat = -Inf,                     Fa = Fa.lo.lim )
  upperBound = c(kon = difflim, kcat = Inf,                      Fa = 1         )
} else {
  lowerBound = c(kon = -Inf,    kcat = -Inf, kph2 = kph2.lo.lim, Fa = Fa.lo.lim )
  upperBound = c(kon = difflim, kcat = Inf,  kph2 = kph2.hi.lim, Fa = 1         )
}

# Generating the parameter init guesses ----
genGuess = function(df.here, modelmode, ignoredGuesses){
  # Generate guess values of the fitted parameters

  rxn.here = unique(df.here$rxn)
  
  # init for Fa is 95% of max fraccleaved
  Fa.guess0 = unname(max(df.here$fraccleaved) * 0.95)
  kon.init = 3/60 # s-1 ## changed # always fast for these ones that we tested
  Fa.guess0 = if_else(Fa.guess0 > 0.4, Fa.guess0, Fa.default)  # account for reactions with very low endpoint
  
  # Get naive guess for NLS to optimize
  naive.guess = c(kon = kon.init,
                  kcat = kcat.init,
                  kph2 = exp(kph2.guess),
                  Fa = Fa.guess0
                  )
  # We need to get a min-1 version to print
  naive.guess.p = naive.guess
  naive.guess.p[names(naive.guess.p)!="Fa"] = naive.guess.p[names(naive.guess.p)!="Fa"] * 60
  print("Naive guess:")
  print(as.numeric(formatCoefs(naive.guess.p)))

  # We now fit an NLS of combined exponential to optimize our init values
  # Append a 100% endpoint to prevent ones without long-time data failing
  #   (NOT done for ODE-fitting. Only for approximation.)
  df.here = df.here %>%
    bind_rows(data.frame(
      rxn = rxn.here,
      conc = inf.conc,
      Rconc = Rconc,
      time = 60*60*24,  # 24 hours
      fraccleaved = 1
    ))
  NLSresults = tryCatch(
    df.here %>%
      ungroup() %>%
      nlsLM(
        # Conc-dependent (series regime)
        fraccleaved ~ Fa  * (1-exp(-time / ( 1/kcat + 1/(kon*conc) ))) +  # using quasi-steady state assumption
                   (1-Fa) * (1-exp(-time * kph2)),
        algorithm = "LM",  # Levenberg-Marquardt
        start = naive.guess,
        upper = c(kon = exp(difflim)*uLim, kcat = kcat.uLim,
                  kph2 = exp(kph2.hi.lim.guess),    Fa = 1*uLim ),
        lower = c(kon = lLim,              kcat = lLim,
                  kph2 = exp(kph2.lo.lim)*(1+lLim), Fa = Fa.lo.lim.guess*(1+lLim) ),
        control = nls.lm.control(maxiter = 100),
        data = .
      ) %>%
        coef(),
    error = function(e){
      # Failure to fit? Fall back to these defaults
      eMsg = conditionMessage(e)
      print(eMsg)
      return(
        NLSresults = c(kon = kon.init,
                       kcat = kcat.init,
                       kph2 = exp(kph2.guess),
                       Fa = Fa.default)
      )
    }
  )

  # Adopt guesses
  kon.guess  = unname(NLSresults["kon"])
  kcat.guess = unname(NLSresults["kcat"])
  kph2.guess = unname(NLSresults["kph2"])
  Fa.corr = 1  # simulations show that Fa is overestimated if there is a lot of slow kon effects
  Fa.guess   = unname(NLSresults["Fa"]) * Fa.corr

  # If NLS switches up kcat and kph2, switch it back
  if(kph2.guess > kcat.guess){
    temp.hold.k = kcat.guess

    kcat.guess = kph2.guess
    kph2.guess = temp.hold.k

    # Fa would be switched too. Fall back to original naive guess
    Fa.guess = Fa.guess0
  }

  # Fa correction
  # Kinetic competitions impose Fa factor on relevant k constants
  if(modelmode %in% c(
    "BNS", "BNR"
  )){kon.guess = kon.guess * Fa.guess}
  if(modelmode %in% c(
    "UNL", "CFI", "CFR", "CFS"
  )){kcat.guess = kcat.guess * Fa.guess}
  if(modelmode %in% c(
    "CTS", "CFS", "BNS"
  )){kph2.guess = kph2.guess * Fa.guess}

  # Output
  guess = c(kon  = log(kon.guess),
            kcat = log(kcat.guess),
            kph2 = log(kph2.guess),
            Fa   = Fa.guess
            )

  outguess = guess[!names(guess) %in% ignoredGuesses]

  # Print guess for ref
  # Linear, and in min-1
  outguess.p = outguess
  outguess.p[names(outguess.p)!="Fa"] = exp(outguess.p[names(outguess.p)!="Fa"]) * 60
  print("NLS guess:")
  print(as.numeric(formatCoefs(outguess.p)))

  return(outguess)
}

# Coefs and goodness-of-fit ----
# CI is 95%, assuming normal distribution around model
CIcalc = function(param, fitpars, val){ if( !is.na(val) & any(!is.na(fitpars)) ) {
  unlist(fitpars[param,2])*1.96
  } else {NA_real_}
}
Pcalc  = function(param, fitpars, val){ if( !is.na(val) & any(!is.na(fitpars)) ) {
  fitpars[param,4]
  } else {NA_real_}
}
runaway = 1e3 # runaway error margin (when no data to fit)
ci_protect = function(ci, runaway.thres) {
  if_else(is.na(ci), NA_real_, if_else(ci > runaway.thres, Inf, ci))
}
extractCoef = function(df.here, fitO, firstpass = T){
  # Extract coefs and their CI and p values

  # Check which rxn this is first
  rxn.this = unique(as.character(df.here$rxn))

  # Get number of data points (sampleN)
  sampleN = nrow(df.here)

  print(
    paste0(
      "Now fitting: ",
      rxn.this,
      " (N = ",
      sampleN,
      ")"
    )
  )

  # print(fitO)

  coefs   = coef(fitO)
  fitpars = tryCatch(
    summary(fitO)$par,
    warning = function(w){
      warnMsg = conditionMessage(w)
      message(warnMsg)
      return(NA)
    }
  )
  # print(fitpars)

  kon  = coefs["kon"]
  kcat = coefs["kcat"]
  kph2 = coefs["kph2"]
  Fa   = coefs["Fa"]

  kon.CI  = ci_protect(CIcalc("kon",  fitpars, kon),  log(runaway))
  kcat.CI = ci_protect(CIcalc("kcat", fitpars, kcat), log(runaway))
  kph2.CI = ci_protect(CIcalc("kph2", fitpars, kph2), log(runaway))
  Fa.CI   = ci_protect(CIcalc("Fa",   fitpars, Fa),       runaway)

  print("Fit output:")
  print(suppressWarnings(as.numeric(formatCoefs(c(
    kon  = exp(kon)*60,
    kcat = exp(kcat)*60,
    kph2 = exp(kph2)*60,
    Fa = Fa
  )))))

  # return in unit of s-1, linear
  return(data.frame(
    kon     = exp(kon),
    kon.hi  = exp(kon  + kon.CI),
    kon.lo  = exp(kon  - kon.CI),
    kon.p   = Pcalc("kon",  fitpars, kon),

    kcat    = exp(kcat),
    kcat.hi = exp(kcat + kcat.CI),
    kcat.lo = exp(kcat - kcat.CI),
    kcat.p  = Pcalc("kcat", fitpars, kcat),

    kph2    = exp(kph2),
    kph2.hi = exp(kph2 + kph2.CI),
    kph2.lo = exp(kph2 - kph2.CI),
    kph2.p  = Pcalc("kph2", fitpars, kph2),

    Fa      =     Fa,
    Fa.hi   =     Fa   + Fa.CI,
    Fa.lo   =     Fa   - Fa.CI,
    Fa.p    = Pcalc("Fa",   fitpars, Fa),

    sampleN = sampleN,

    row.names = rxn.this
  ))
}

###############################################################

# FITTING TIME ----

MAXITER = 200  # maximum fitting iterations

# Compile ignored coefs for args for re-fitting
givenplat = Fa.default

get.ignore.here = function(y){
  lim.here = y %>% select(
    rapidon,  # no fit kon
    nobipha,  # no fit kph2
    fixplat   # no fit Fa
    ) %>% distinct()

  ignore.here = c()
  if(lim.here$rapidon){
    ignore.here = c(ignore.here, "kon")
  }
  if(lim.here$nobipha){
    ignore.here = c(ignore.here, "kph2")
  }
  if(lim.here$fixplat){
    ignore.here = c(ignore.here, "Fa")
  }

  if(modelmode == "MAX") ignore.here = c(ignore.here, "kph2") %>% unique()

  ignore.here
}
lim.guess.wrapper = function(y){
  ignore.here = get.ignore.here(y)
  genGuess(y, modelmode, ignoredGuesses = ignore.here)
}
lim.bound.wrapper = function(y){
  ignore.here = get.ignore.here(y)
  list(lowerBound[!names(lowerBound) %in% ignore.here],
       upperBound[!names(upperBound) %in% ignore.here])
}
lim.static.wrapper = function(y){
  ignore.here = get.ignore.here(y)
  full.static = c(kon = difflim,
                  kph2 = log(0),
                  Fa = givenplat
  )
  full.static[names(full.static) %in% ignore.here]
}

## First pass ----
ODEfit = df.dt %>%
  ## Remove special data points
  mutate(rapidon = FALSE,
         nobipha = FALSE,
         fixplat = FALSE) %>%
  group_by(rxn) %>%
  do(extractCoef(.,
                 modFit(
                   f = ODEcost,
                   p = lim.guess.wrapper(.),
                   df.data = .,
                   params.static = lim.static.wrapper(.),
                   lower = lim.bound.wrapper(.)[[1]],
                   upper = lim.bound.wrapper(.)[[2]],
                   method = "Marq",
                   jac = NULL,
                   control = list(
                     nprint = 0,
                     maxiter = MAXITER
                   ),
                   hessian = T
                 )
  )) %>%
  mutate(rapidon = F, nobipha = F, fixplat = F)

## Identify ones that went wrong ----
# Send to second pass with kon at diffusion limit

## Second pass (no kon effects) ----
secondPass = ODEfit %>%
  mutate(
    rapidon = any(
      is.na(kon.lo), 
      is.infinite(kon.hi), 
      kon.p > 0.5,  
      kon > exp(difflim) * 0.95 
    ),
    nobipha = F,
    fixplat = F,
    firstPassOk = !rapidon
  ) %>%
  filter(!firstPassOk) %>%
  select(rxn, rapidon, nobipha, fixplat)

secondPass.rxns = as.character(secondPass$rxn)

if(length(secondPass.rxns) > 0){
  print("Now fitting a second pass for the following:")
  cat(paste(c(secondPass.rxns, ""), collapse = "\n"))

  ODEfit.r = df.dt %>%
    filter(rxn %in% secondPass.rxns) %>%
    left_join(secondPass, by = "rxn") %>%
    group_by(rxn) %>%
    do(extractCoef(.,
                   modFit(
                     f = ODEcost,
                     p = lim.guess.wrapper(.),
                     df.data = .,
                     params.static = lim.static.wrapper(.),
                     lower = lim.bound.wrapper(.)[[1]],
                     upper = lim.bound.wrapper(.)[[2]],
                     method = "Marq",
                     jac = NULL,
                     control = list(
                       nprint = 0,
                       maxiter = MAXITER
                     ),
                     hessian = T
                   ),
                   firstpass = F
      ) %>%
      mutate(rapidon = is.na(kon), nobipha = is.na(kph2),
             fixplat = F)) %>%
    mutate(
      kon = if_else(rapidon, exp(difflim), kon),
      kph2 = if_else(nobipha, 0, kph2)
      )  # give the fixed values back to the panel

  # Combine ODEfit tables
  ODEfit.first.bu = ODEfit
  ODEfit = ODEfit %>%
    filter(!rxn %in% secondPass.rxns) %>%
    rbind(ODEfit.r)
}

## Identify ones that went wrong ----
# Send to third pass with kph2 no longer considered

## Third pass (no second phase) ----
thirdPass = ODEfit %>%
  mutate(
    rapidon = kon == exp(difflim),  # inherit rapidon designation
    nobipha = any(
      kph2.p > 0.5,
      kph2 < exp(kph2.lo.lim) * 1.05
    ),
    secondPassOk = !nobipha,
    fixplat = F
  ) %>%
  filter(!secondPassOk) %>%
  select(rxn, rapidon, nobipha,
         fixplat)

thirdPass.rxns = as.character(thirdPass$rxn)

if(length(thirdPass.rxns) > 0 & modelmode != "MAX"){
  print("Now fitting a third pass for the following:")
  cat(paste(c(thirdPass.rxns, ""), collapse = "\n"))

  ODEfit.r3 = df.dt %>%
    filter(rxn %in% thirdPass.rxns) %>%
    left_join(thirdPass, by = "rxn") %>%
    group_by(rxn) %>%
    do(extractCoef(.,
                   modFit(
                     f = ODEcost,
                     p = lim.guess.wrapper(.),
                     df.data = .,
                     params.static = lim.static.wrapper(.),
                     lower = lim.bound.wrapper(.)[[1]],
                     upper = lim.bound.wrapper(.)[[2]],
                     method = "Marq",
                     jac = NULL,
                     control = list(
                       nprint = 0,
                       maxiter = MAXITER
                     ),
                     hessian = T
                   ),
                   firstpass = F
      ) %>%
      mutate(rapidon = is.na(kon), nobipha = is.na(kph2),
             fixplat = F)) %>%
    mutate(
      kon = if_else(rapidon, exp(difflim), kon),
      kph2 = if_else(nobipha, 0, kph2)
    )  # give the fixed values back to the panel

  # Combine ODEfit tables
  ODEfit.second.bu = ODEfit
  ODEfit = ODEfit %>%
    filter(!rxn %in% thirdPass.rxns) %>%
    rbind(ODEfit.r3)
}

####################################################################

# SETUP FOR PLOTS ----

## Get fit graphs ----
expandTime = 1/1

ODEcalc.graph = function(d){
  # Given a one-row df of setup parameters, generate an ODE trajectory

  t.plot = timeGen(d$Tmax*(1+expandTime))

  params.static = c(
    kon =  log(d$kon),
    kcat = log(d$kcat),
    kph2 = log(d$kph2),
    Fa =   d$Fa
    )

  return(
    ODEcalc(
      state.setup = c(conc = d$conc, Rconc = d$Rconc),
      times = t.plot,
      parameters = c(),
      params.static = params.static
    ) %>%
      filter(time > 0)  # remove 0 time point from sim; we are plotting in log
  )
}

# Get Tmax for plot
df.times = df.dt %>%
  ungroup() %>%
  select(rxn, time) %>%
  group_by(rxn) %>%
  summarize_all(max) %>%
  rename(Tmax = time)

ODEfit.plot = ODEfit %>%
  inner_join(., df.times, by = c("rxn"))

# Expe concs
df.concs = df.dt %>%
  ungroup() %>%
  select(rxn, conc) %>%
  distinct() %>%
  mutate(Rconc = Rconc)
df.graph.setup = ODEfit.plot %>%
  select(-ends_with(match = c("\\.lo","\\.hi","\\.p")), -rapidon, -nobipha) %>%
  inner_join(df.concs, ., by = c("rxn"))
df.graph = df.graph.setup %>%
  group_by(rxn, conc, Rconc) %>%
  do(ODEcalc.graph(.)) %>%
  separate_wider_delim(rxn, "_", names = c("miR", "preincubation"))

# Limit line (no kon delay, ~ Inf conc)
n.rxns    = length(rxns)
df.concs.inf = data.frame(
  rxn   = rxns,
  conc  = rep(inf.conc,   times = n.rxns),
  Rconc = rep(Rconc, times = n.rxns)
)
df.graph.setup.inf = ODEfit.plot %>%
  select(-ends_with(match = c("\\.lo","\\.hi","\\.p")), -rapidon, -nobipha) %>%
  inner_join(df.concs.inf, ., by = c("rxn"))
df.graph.inf = df.graph.setup.inf %>%
  group_by(rxn, conc, Rconc) %>%
  do(ODEcalc.graph(.)) %>%
  mutate(conc = Inf) %>%
  separate_wider_delim(rxn, "_", names = c("miR", "preincubation"))

###########################################################################

# PLOT ----
FACETS = 4
plt = df.dt %>%
  group_by(miR) %>% mutate(concset = if_else(conc == max(conc), "high", "low")) %>% ungroup() %>%
  ggplot(aes(x = time/60, y = fraccleaved,
             # color = preincubation,
             color = paste(preincubation, concset),
             group = factor(paste(preincubation, conc))
             )) +
  coord_cartesian(xlim = c(0, 25), ylim = c(0,1)) +
  
  geom_line(data = df.graph %>%
              group_by(miR) %>% mutate(concset = if_else(conc == max(conc), "high", "low")) %>% ungroup(),
            aes(y = fraccleaved.ode), 
            linetype = "dashed",
            linewidth = 0.75) +
  geom_line(data = df.graph.inf %>%
              mutate(concset = "inf"),
            aes(y = fraccleaved.ode), 
            linewidth = 0.75) +
  geom_point(size = 1.2, stroke = 1.5, shape = 1) +
  
  geom_text(data = ODEfit.plot %>%
              separate_wider_delim(rxn, "_", names = c("miR", "preincubation")),
            inherit.aes = FALSE, show.legend = F,
            size = 4.5,
            x = Inf,
            y = -Inf,
            color = "black",
            hjust = 1, vjust = -0.2,
            aes(
              label = paste0(preincubation, ":n = ", sampleN)
            )) +
  
  facet_wrap("miR", ncol = FACETS, scales = "free") +
  
  # scale_color_manual(
  #   breaks = c("none", "noCa", "wCa"),
  #   labels = c("No preincubation", "0 mM Ca2+", "1 mM Ca2+"),
  #   values = c("darkslateblue", "darkgoldenrod4", "firebrick"),
  #   name = "Preincubation for 1 h without Mg2+"
  # ) +
  
  scale_color_manual(
    breaks = c("none inf", 
               "none high", 
               "none low", 
               "noCa inf", 
               "noCa high", 
               "noCa low", 
               "wCa inf", 
               "wCa high", 
               "wCa low"),
    values = c("steelblue4",
               "steelblue",
               "steelblue1",
               
               "darkgoldenrod4",
               "darkgoldenrod",
               "lightgoldenrod3",
               
               "red4",
               "orangered3",
               "lightcoral"
               ),
    name = "Preincubation for 1 h without Mg2+"
  ) +
  
  scale_x_continuous(
    breaks = seq(0, 200, 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0),
    breaks = seq(0, 1, 0.2)) +
  xlab("Time (min)") +
  ylab("Fraction sliced") +
  theme0
plt

##############################################################
# SAVING ----

## Save plot ----
if(SAVE_PLOT){
  ggsave(
    filename = paste0(
      "./",
      "slicing_preincubationing_plot",
      ".pdf"),
    plot = plt,

    width = 11, height = 3,

    units = "in"
  )
}

## Save output values ----
if(SAVE_TABLE){
  ODEfit.write = ODEfit %>%
    ungroup() %>%
    transmute(
      rxn       = rxn,

      kon       = kon*60,
      kon.lo    = kon.lo*60,
      kon.hi    = kon.hi*60,
      kon.pP    = -log10(kon.p),

      kcat      = kcat*60,
      kcat.lo   = kcat.lo*60,
      kcat.hi   = kcat.hi*60,
      kcat.pP   = -log10(kcat.p),

      kph2      = kph2*60,
      kph2.lo   = kph2.lo*60,
      kph2.hi   = kph2.hi*60,
      kph2.pP   = -log10(kph2.p),

      Fa        = Fa,
      Fa.lo     = Fa.lo,
      Fa.hi     = Fa.hi,
      Fa.pP     = -log10(Fa.p),

      sampleN = sampleN
    )
  write.csv(ODEfit.write, paste0(
    "./",
    "slicing_ks_preincubation_raw",
    ".csv"),
    quote = F, row.names = F
  )
}

#@#@#
print("Done!")
