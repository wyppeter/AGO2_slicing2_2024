# For RISC quantification by fitting to quadratic titration curves
# Bartel Lab
# Peter Y. Wang 2024

# LIBRARIES ----
library(tidyverse)
library(minpack.lm)
library(ggpubr)

# Set ggplot2 theme ----
theme0 = theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(colour = "black"),
  axis.text = element_text(color = "black", size = 12),
  axis.title = element_text(color = "black", size = 12),
  legend.background = element_blank(),
  legend.key = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_text(size = 12, colour = "black"),
  strip.text.y = element_text(size = 12, colour = "black", angle = 0)
)

# READ DATA ----

# quantoligoconc = 1000  # pM; miR-200b
quantoligoconc = 2000  # pM; let-7a

df = read.csv("titration_RISC_quant_data.csv", stringsAsFactors = T) %>%
  # sort by max fracbound for visualization
  filter(
    # miR == "miR-200b"
    miR == "let-7a"
      ) %>%
  select(-miR) %>%
  group_by(AGOprep) %>%
  arrange(desc(max(fracbound)), dil) %>%
  mutate(AGOprep = fct_inorder(AGOprep))

AGOpreps = df$AGOprep %>% unname() %>% unique() %>% factor(levels = c("WT","2A","5E","2A5E"))

# Equation for fraction bound ----
calcFB = function(At, dil, Rt, Kd, Fmax = 1){
  (((At*dil+Rt+10^Kd) - sqrt((At*dil+Rt+10^Kd)^2 - 4*At*dil*Rt))/(2*Rt)) * Fmax
}

# NLS modeling ----

## NLS fitting ----
nlsFit = df %>%
  mutate(Rt = quantoligoconc) %>%
  do(data.frame(model = broom::tidy(
    nlsLM(
      
      ## miR-200b
      # formula = fracbound ~ (((At*dil+Rt+10^Kd) - sqrt((At*dil+Rt+10^Kd)^2 - 4*At*dil*Rt))/(2*Rt)),
      # start   = c(At = 100 * 1000,   Kd = 1.5),
      # upper   = c(At = 1000 * 1000,  Kd = 2.5),
      # lower   = c(At = 0,            Kd = 0.5),
      
      # ## let-7a
      formula = fracbound ~ (((At*dil+Rt+10^Kd) - sqrt((At*dil+Rt+10^Kd)^2 - 4*At*dil*Rt))/(2*Rt)) * Fmax,
      start   = c(At = 100 * 1000,   Kd = 1.5, Fmax = 0.9),
      upper   = c(At = 2000 * 1000,  Kd = 10, Fmax = 1),
      lower   = c(At = 0,            Kd = -10, Fmax = 0),
      
      control = c(maxiter = 1000),
      data = .
    )
  ))) %>%
  ungroup() %>%
  transmute(AGOprep = AGOprep,
            coeff = model.term,
            v  = model.estimate,
            hi = model.estimate+model.std.error*1.96,
            lo = model.estimate-model.std.error*1.96,
            p = model.p.value) %>%
  pivot_wider(id_cols = AGOprep, names_from = coeff, names_sep = "_", values_from = c(v, hi, lo, p)) %>%
  mutate(AGOprep = factor(AGOprep, levels = AGOpreps)) %>%
  arrange(AGOprep)

## Create NLS fit model graph ----
df.maxconc = df %>%
  group_by(AGOprep) %>%
  summarize_at("dil", max) %>%
  rename(max.dil = dil)
xmax = 0.5
xstep = 0.0001
nlsFit.graph.raw = data.frame(
  dil     =         rep( seq(0, xmax, xstep), times = length(AGOpreps)),
  AGOprep = factor( rep( AGOpreps,            each  = xmax/xstep+1),      levels = AGOpreps)
) %>%
  left_join(df.maxconc, by = "AGOprep")
df.coefs = nlsFit %>%
  transmute(AGOprep = AGOprep,
            At = v_At,
            Kd = v_Kd,
            Rt = quantoligoconc,
            Fmax = v_Fmax
            )
nlsFit.graph = nlsFit.graph.raw %>%
  left_join(df.coefs, by = "AGOprep") %>%
  mutate(fracbound = calcFB(
    
    # At, dil, Rt, Kd
    At, dil, Rt, Kd, Fmax
    
    ))

# Plotting time
df %>%
  ggplot(aes(x = dil, y = fracbound, color = AGOprep, group = AGOprep)) +
  geom_hline(yintercept = c(0,1)) +
  geom_vline(xintercept = c(0)) +
  geom_line(data = nlsFit.graph,
            aes(group = AGOprep),
            linetype = "dashed",
            linewidth = 0.5, alpha = 0.75,
            show.legend = T) +

  geom_point(shape = 1, size = 1.5, stroke = 1.25, alpha = 0.8,
             show.legend = F) +

  scale_x_continuous(labels = scales::percent,
                     expand = expansion(mult = c(0.05, 0.2), add = 0)) +

  coord_cartesian(
    
    # xlim = c(0, 0.006),
    # ylim = c(0, 0.75)
    
    xlim = c(0, 0.025),
    ylim = c(0, 1)

    ) +
  labs(x = "Dilution from stock", y = "Fraction bound of target RNA") +

  theme0

