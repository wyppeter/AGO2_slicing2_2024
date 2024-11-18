# Slicing assay re-graphing
# Requires helper file ODElib.R
# Bartel Lab
# Peter Y. Wang 2024

# LIBRARIES ----
# tidyverse
library(tidyverse)
# graphing
library(viridis)
library(ggpubr)
# ODEs
library(deSolve)
library(FME)

# SELECT MODEL MODE ----
modelmode = "CTS"

# DISPLAY SETUP ----
options(scipen = 5)

# CONSTANTS ----
difflim = log(6/60)  # given known miR-7 kon; nM-1 s-1 (10-10 M-1 s-1); diffusion limit; 6 in nM-1 min-1, 0.1 in nM-1 s-1
inf.conc = 1e6  # nM; a large number to graph limit behavior (no kon effects)
Rconc = 0.05  # nM

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
df.raw = read.csv("slicing_data_p2.csv",
  stringsAsFactors = T) %>%
  filter(targ == "perfect") %>%
  mutate(rxn = fct_inorder(factor(paste0(miR, "#", AGOmut, "_", targ)))) %>%
  select(-miR, -targ, -AGOmut) %>%
  mutate(Rconc = Rconc)

# PROCESS AND COMBINE DATA ----
df = df.raw %>%
  ## MISC CLEANUP
  mutate(time = as.integer(round(time*60))) %>%  # rounding to the time resolution of ODE, and convert to seconds
  filter(conc != 0, time != 0) %>%  # discard control points
  mutate(conc_d = factor(conc)) %>%  # factor-ize some variables
  arrange(rxn, conc)

# Note down the guide names used here
rxns = as.character(unname(levels(df$rxn)))

# Simplify df
df.dt = df %>%
  select(rxn, conc, Rconc, time, fraccleaved) %>%
  mutate(conc_f = factor(conc))
df.dt.lowconc = df.dt %>%
  mutate(conc = 2)

# Get ktable
df.k = read.csv(
  "slicing_ks_p2_raw.csv",
  stringsAsFactors = T 
)
ODEfit = df.k %>%
  transmute(
    rxn = rxn,
    kon = log(kon/60),
    kcat = log(kcat/60),
    kph2 = log(kph2/60),
    Fa = Fa,
    sampleN = sampleN
  ) %>%
  filter(as.character(rxn) != sub("_perfect", "", as.character(rxn)))

df.ft.lowconc.withk = df.dt.lowconc %>%
  left_join(ODEfit %>% select(
    -kon, -sampleN
  ), by = "rxn") %>%
  filter(time < (log(2)/exp(kcat)) * 2)  # only look at short time points relative to slicing rate

# Reference of WT's kon
WTkon = 4.03020745277037/60


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
  seq(Tthres[3], Tend,        60)   # 1 minute 
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
  # MODIFIED TO INCORPORATE KON LOSS SIM

  # Define A_T and R_T
  A_T = state.setup["conc"]  %>% unname()
  R_T = state.setup["Rconc"] %>% unname()

  # Include static params
  sim.kon = c(kon = log(exp(unname(parameters["konloss"])) * WTkon))
  paramsplus = c(sim.kon, params.static)
  
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

ODEgap = function(konloss, df.this){
  # Given the data of a miR at a particular conc set, calculate the ODE outcome,
  # and provide...
  # The least negative deviation, or the smallest positive deviation when you overshoot
  
  # Generate time sequence
  Tmax = max(df.this$time)
  t = timeGen(Tmax)
  
  # Get precalc ks
  params.static = df.this %>%
    select(kcat, kph2, Fa) %>%
    distinct() %>%
    unlist()
  
  # Calculate ODE simulation and compare
  overlaid = ODEcalc(
    state.setup = initConc(df.this),
    times = t,
    parameters = konloss,
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
    stop("Crit. error at ODEgap!")
  }
  
  resids = overlaid$residuals
  
  # return
  if(any(resids > 0)){
    min(resids[resids > 0])^2
  } else {
    min(-resids[resids <= 0])^2
  }
}

lossfit = df.ft.lowconc.withk %>%
  group_by(rxn) %>%
  do(., {
    d.here = .
    optim.out = optim(
    par = c(konloss = log(0.1)),
    fn = ODEgap,
    df.this = d.here,
    method = "Nelder-Mead",
    )
    if(optim.out["message"] != "NULL") print(optim.out["message"])
    data.frame(
      konloss = optim.out["par"] %>% unlist() %>% unname(),
      resid = optim.out["value"] %>% unlist() %>% unname()
    )
  }) %>%
  mutate(kon.foldloss.limit = exp(konloss))

write.csv(lossfit,
          "./kon-loss-lower-lim.csv",
          row.names = F)

