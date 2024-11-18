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

# Get ktable
df.k = read.csv(
  "slicing_ks_p2_raw.csv",
  stringsAsFactors = T 
)
ODEfit = df.k %>%
  transmute(
    rxn = rxn,
    kon = kon/60,
    kcat = kcat/60,
    kph2 = kph2/60,
    Fa = Fa,
    sampleN = sampleN
  )

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

# SETUP FOR PLOTS ----

## Get fit graphs ----
expandTime = 1/3

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
  inner_join(df.concs, ., by = c("rxn"))
df.graph = df.graph.setup %>%
  group_by(rxn, conc, Rconc) %>%
  do(ODEcalc.graph(.)) %>%
  mutate(conc_f = factor(conc))

# Limit line (no kon delay, ~ Inf conc)
n.rxns    = length(rxns)
df.concs.inf = data.frame(
  rxn   = rxns,
  conc  = rep(inf.conc,   times = n.rxns),
  Rconc = rep(Rconc, times = n.rxns)
)
df.graph.setup.inf = ODEfit.plot %>%
  inner_join(df.concs.inf, ., by = c("rxn"))
df.graph.inf = df.graph.setup.inf %>%
  group_by(rxn, conc, Rconc) %>%
  do(ODEcalc.graph(.)) %>%
  mutate(conc = Inf,
         conc_f = factor(Inf))

# Slow kon simulations
WTkon = 4.03020745277037/60
rxns.to.sim = c("miR-7-L22#H56A_perfect",
                "miR-7-L22#R68A_perfect",
                "miR-7-L22#R97A_perfect",
                "miR-7-L22#K98A_perfect",
                "miR-7-L22#H56A.K98A_perfect",
                "miR-7-L22#R97A.K98A_perfect",
                "miR-7-L22#quadA_perfect",
                "miR-7-L22#R97E.K98E_perfect",
                
                "miR-7-L22#R710A_perfect",
                "miR-7-L22#H712A_perfect",
                "miR-7-L22#R710A.H712A_perfect"
                )
slower.kons = data.frame(slow.fold = c(10)) %>%
  mutate(kon.wt = WTkon,
         kon = kon.wt/slow.fold) %>%
  crossing(data.frame(rxn = rxns.to.sim))
df.graph.setup.slowkonsim = df.graph.setup %>%
  filter(rxn %in% rxns.to.sim) %>%
  select(-kon) %>%
  left_join(slower.kons, by = "rxn", relationship = "many-to-many")
df.graph.slowkonsim = df.graph.setup.slowkonsim %>%
  group_by(rxn, conc, Rconc, slow.fold) %>%
  do(ODEcalc.graph(.)) %>%
  mutate(conc_f = factor(conc))

###########################################################################

# PLOT ----
# FACETS = as.integer(round(sqrt(n.rxns*3/2)))
# FACETS = 2
FACETS = 3
rxn.plot.sublist = c(
  "miR-7-L22#WT_perfect",
  
  "miR-7-L22#H56A_perfect",
  "miR-7-L22#R68A_perfect",
  "miR-7-L22#R97A_perfect",
  "miR-7-L22#K98A_perfect",
  "miR-7-L22#H56A.K98A_perfect",
  "miR-7-L22#R97A.K98A_perfect",
  "miR-7-L22#quadA_perfect",
  "miR-7-L22#R97E.K98E_perfect"
  
  # "miR-7-L22#WT_16bp",
  # "miR-7-L22#H56A_16bp",
  # "miR-7-L22#R68A_16bp",
  # "miR-7-L22#R97A.K98A_16bp",
  # "miR-7-L22#R97E.K98E_16bp",
  # "miR-7-L22#quadA_16bp"

  # "miR-7-L22#WT_perfect",
  # "miR-7-L22#H712A_perfect",
  # "miR-7-L22#R710A_perfect",
  # "miR-7-L22#R710A.H712A_perfect"
  
  # "miR-7-L22#WT_mm10",
  # "miR-7-L22#H712A_mm10",
  # "miR-7-L22#R710A_mm10",
  # "miR-7-L22#R710A.H712A_mm10"
  
  # "miR-7-L22#WT_mm11",
  # "miR-7-L22#H712A_mm11",
  # "miR-7-L22#R710A_mm11",
  # "miR-7-L22#R710A.H712A_mm11"
  
  # "miR-7-L22#WT(Sept)_perfect", "miR-7-L22#R710A(Sept)_perfect", "miR-7-L22#H712A(Sept)_perfect", "miR-7-L22#R635A(Sept)_perfect",
  # "miR-7-X3#WT(Sept)_perfect", "miR-7-X3#R710A(Sept)_perfect", "miR-7-X3#H712A(Sept)_perfect", "miR-7-X3#R635A(Sept)_perfect"
)
df.dt %>%
  filter(rxn %in% rxn.plot.sublist) %>% mutate(rxn = factor(rxn, levels = rxn.plot.sublist)) %>%
  ggplot(aes(x = time/60, y = fraccleaved,
             
             color = conc_f, group = conc_f
             # color = factor(paste(rxn, conc)), group = factor(paste(rxn, conc))
             
             )) +
  coord_cartesian(xlim = c(
    0, 4
    # 0, 70
    # 0, 11
    # 0, 250
    )) +
  geom_line(data = df.graph %>%
              filter(rxn %in% rxn.plot.sublist) %>% mutate(rxn = factor(rxn, levels = rxn.plot.sublist)),
            aes(y = fraccleaved.ode),
            linewidth = 0.75) +
  geom_line(data = df.graph.inf %>%
              filter(rxn %in% rxn.plot.sublist) %>% mutate(rxn = factor(rxn, levels = rxn.plot.sublist)),
            aes(y = fraccleaved.ode),
            linewidth = 0.75) +
  geom_point(size = 1.2, stroke = 1.5, shape = 1) +
  geom_text(data = ODEfit.plot %>%
              filter(rxn %in% rxn.plot.sublist) %>% mutate(rxn = factor(rxn, levels = rxn.plot.sublist)),
            inherit.aes = FALSE,
            size = 4.5,
            color = "black",
            x = Inf,
            y = -Inf,
            hjust = 1, vjust = -0.2,
            aes(
              label = paste0(rxn, ":n = ", sampleN),
              # y = 0+as.numeric(rxn)*0.02
            )) +
  
  facet_wrap("rxn", ncol = FACETS, scales = "free") +

  scale_color_manual(breaks = c(2, 5, 10, Inf),
                     # values = c("gray75", "gray50", "gray25", "black"),
                     values = c("gray65", "gray45", "gray25", "black"),
                     name = "[RISC] (nM)") +
  
  # ###/\/\/\/\/
  geom_line(data = df.graph.slowkonsim %>%
              filter(rxn %in% rxn.plot.sublist) %>% mutate(rxn = factor(rxn, levels = rxn.plot.sublist)),
            aes(y = fraccleaved.ode),
            linewidth = 0.75, linetype = "dashed") +
  # ###/\/\/\/\/
  
  # ###/\/\/\/\/
  # geom_line(data = df.graph.slowkonsim %>%
  #             filter(rxn %in% rxn.plot.sublist) %>% mutate(rxn = factor(rxn, levels = rxn.plot.sublist)),
  #           aes(y = fraccleaved.ode,
  #               group = interaction(conc_f, slow.fold),
  #               color = factor(paste(rxn, conc))),
  #           linewidth = 0.75, linetype = "dashed") +
  # scale_alpha_manual(breaks = c(2, 5, 10),
  #                    values = c(0.2, 0.5, 0.8),
  #                    name = "[RISC] (nM)") +
  # ###/\/\/\/\/

  # scale_color_manual(breaks = c("miR-7-L22#WT_perfect Inf",          "miR-7-L22#WT_perfect 10",          "miR-7-L22#WT_perfect 5",          "miR-7-L22#WT_perfect 2",
  #                               "miR-7-L22#R710A_perfect Inf",       "miR-7-L22#R710A_perfect 10",       "miR-7-L22#R710A_perfect 5",       "miR-7-L22#R710A_perfect 2",
  #                               "miR-7-L22#H712A_perfect Inf",       "miR-7-L22#H712A_perfect 10",       "miR-7-L22#H712A_perfect 5",       "miR-7-L22#H712A_perfect 2",
  #                               "miR-7-L22#R710A.H712A_perfect Inf", "miR-7-L22#R710A.H712A_perfect 10", "miR-7-L22#R710A.H712A_perfect 5", "miR-7-L22#R710A.H712A_perfect 2"
  #                               ),
  #                    values = c(
  #                      "black", "gray30", "gray50", "gray70",
  #                      "steelblue4", "steelblue", "steelblue2", "skyblue1",
  #                      "violetred4", "palevioletred3", "palevioletred1", "#FFBBCC",
  #                      "purple4", "purple3", "mediumpurple3", "mediumpurple1"
  #                      ),
  #                    name = "[RISC] (nM)") +

  # scale_color_manual(breaks = c(
  #   "miR-7-L22#WT_mm10 Inf",          "miR-7-L22#WT_mm10 10",          "miR-7-L22#WT_mm10 5",          "miR-7-L22#WT_mm10 2",
  #   "miR-7-L22#R710A_mm10 Inf",       "miR-7-L22#R710A_mm10 10",       "miR-7-L22#R710A_mm10 5",       "miR-7-L22#R710A_mm10 2",
  #   "miR-7-L22#H712A_mm10 Inf",       "miR-7-L22#H712A_mm10 10",       "miR-7-L22#H712A_mm10 5",       "miR-7-L22#H712A_mm10 2",
  #   "miR-7-L22#R710A.H712A_mm10 Inf", "miR-7-L22#R710A.H712A_mm10 10", "miR-7-L22#R710A.H712A_mm10 5", "miR-7-L22#R710A.H712A_mm10 2",
  #   "miR-7-L22#WT_mm11 Inf",          "miR-7-L22#WT_mm11 10",          "miR-7-L22#WT_mm11 5",          "miR-7-L22#WT_mm11 2",
  #   "miR-7-L22#R710A_mm11 Inf",       "miR-7-L22#R710A_mm11 10",       "miR-7-L22#R710A_mm11 5",       "miR-7-L22#R710A_mm11 2",
  #   "miR-7-L22#H712A_mm11 Inf",       "miR-7-L22#H712A_mm11 10",       "miR-7-L22#H712A_mm11 5",       "miR-7-L22#H712A_mm11 2",
  #   "miR-7-L22#R710A.H712A_mm11 Inf", "miR-7-L22#R710A.H712A_mm11 10", "miR-7-L22#R710A.H712A_mm11 5", "miR-7-L22#R710A.H712A_mm11 2"
  #   ),
  #   values = c(
  #     "black", "gray30", "gray50", "gray70",
  #     "steelblue4", "steelblue", "steelblue2", "skyblue1",
  #     "violetred4", "palevioletred3", "palevioletred1", "#FFBBCC",
  #     "purple4", "purple3", "mediumpurple3", "mediumpurple1",
  #     "black", "gray30", "gray50", "gray70",
  #     "steelblue4", "steelblue", "steelblue2", "skyblue1",
  #     "violetred4", "palevioletred3", "palevioletred1", "#FFBBCC",
  #     "purple4", "purple3", "mediumpurple3", "mediumpurple1"
  #   ),
  #   name = "[RISC] (nM)") +

  # scale_color_manual(breaks = c(
  #   "miR-7-L22#WT(Sept)_perfect Inf",    "miR-7-L22#WT(Sept)_perfect 10",    "miR-7-L22#WT(Sept)_perfect 5",    "miR-7-L22#WT(Sept)_perfect 2",   
  #   "miR-7-L22#R710A(Sept)_perfect Inf", "miR-7-L22#R710A(Sept)_perfect 10", "miR-7-L22#R710A(Sept)_perfect 5", "miR-7-L22#R710A(Sept)_perfect 2",
  #   "miR-7-L22#H712A(Sept)_perfect Inf", "miR-7-L22#H712A(Sept)_perfect 10", "miR-7-L22#H712A(Sept)_perfect 5", "miR-7-L22#H712A(Sept)_perfect 2",
  #   "miR-7-L22#R635A(Sept)_perfect Inf", "miR-7-L22#R635A(Sept)_perfect 10", "miR-7-L22#R635A(Sept)_perfect 5", "miR-7-L22#R635A(Sept)_perfect 2",
  #   "miR-7-X3#WT(Sept)_perfect Inf",     "miR-7-X3#WT(Sept)_perfect 10",     "miR-7-X3#WT(Sept)_perfect 5",     "miR-7-X3#WT(Sept)_perfect 2",    
  #   "miR-7-X3#R710A(Sept)_perfect Inf",  "miR-7-X3#R710A(Sept)_perfect 10",  "miR-7-X3#R710A(Sept)_perfect 5",  "miR-7-X3#R710A(Sept)_perfect 2", 
  #   "miR-7-X3#H712A(Sept)_perfect Inf",  "miR-7-X3#H712A(Sept)_perfect 10",  "miR-7-X3#H712A(Sept)_perfect 5",  "miR-7-X3#H712A(Sept)_perfect 2", 
  #   "miR-7-X3#R635A(Sept)_perfect Inf",  "miR-7-X3#R635A(Sept)_perfect 10",  "miR-7-X3#R635A(Sept)_perfect 5",  "miR-7-X3#R635A(Sept)_perfect 2"
  #   ),
  # values = c(
  #   "black", "gray30", "gray50", "gray70",
  #   "steelblue4", "steelblue", "steelblue2", "skyblue1",
  #   "violetred4", "palevioletred3", "palevioletred1", "#FFBBCC",
  #   "darkgreen", "palegreen4", "palegreen3", "darkolivegreen2",
  #   "black", "gray30", "gray50", "gray70",
  #   "steelblue4", "steelblue", "steelblue2", "skyblue1",
  #   "violetred4", "palevioletred3", "palevioletred1", "#FFBBCC",
  #   "darkgreen", "palegreen4", "palegreen3", "darkolivegreen2"
  # ),
  # name = "[RISC] (nM)") +

  scale_x_continuous(
    breaks = seq(0, 200, 1),
    # breaks = seq(0, 200, 20),
    # breaks = seq(0, 200, 5),
    # breaks = seq(0, 200, 50),
    expand = c(0, 0)
    ) +
  scale_y_continuous(
    limits = c(0, 1),
    # limits = c(0, .45),
    
    expand = c(0, 0),
    breaks = seq(0, 1, 0.2)) +
  xlab("Time (min)") +
  ylab("Fraction sliced") +
  theme0
  
