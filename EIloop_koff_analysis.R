# Dissociation kinetics analysis
# Bartel Lab
# Peter Y. Wang 2024

library(tidyverse)
library(minpack.lm)

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

muts_scale_breaks = c("WT", "EI_5xE")
muts_scale_labels = c("WT", "EI 5xE")
muts_scale_values = c("gray25", "indianred")

# Data processing ----
df.fb = read.csv("EIloop_koff_data.csv") %>%
  mutate(AGOmut = factor(AGOmut, levels = c("WT", "EI_5xE")),
         expe = factor(expe))

# NLS fit ----
df.fit = df.fb %>%
  filter(time > 0) %>%
  group_by(miR, AGOmut) %>%
  do(data.frame(model = broom::tidy(
    nlsLM(
      formula = fracbound ~ A * exp(-time * exp(koff)) + void,
      start   = c(A = 0.75, koff = 0.2, void = 0.05),
      upper   = c(A = 1, koff = Inf, void = 1),
      lower   = c(A = 0, koff = -Inf, void = 0),
      control = c(maxiter = 1000),
      data = .
    )
  )) %>%
    mutate(miR = .$miR %>% unique())) %>%
  arrange(model.term, miR, AGOmut)

# Make graph ----
df.graph = df.fit %>%
  select(miR, AGOmut, model.term, model.estimate) %>%
  pivot_wider(names_from = model.term, values_from = model.estimate) %>%
  crossing(data.frame(time = 10^seq(-4, 3, 0.01))) %>%
  mutate(fracbound = A * exp(-time * exp(koff)) + void)

# Plot time course ----
# MIRHERE = "miR-200b"
MIRHERE = "let-7a"
df.fb %>%
  filter(time > 0,
         miR == MIRHERE
         ) %>%
  ggplot(aes(x = time, 
             y = fracbound,
             color = AGOmut)) +
  geom_line(data = df.graph %>%
              filter(miR == MIRHERE),
            linewidth = 0.75,
            alpha = 0.75) +
  geom_point(
    aes(shape = expe),
    size = 2,
    alpha = 0.75) +
  facet_wrap("miR") +
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)) +
  scale_x_continuous(
    # breaks = seq(0, 1000, 5),
    breaks = seq(0, 1000, 60),
    expand = c(0, 0)) +
  coord_cartesian(
    ylim = c(0.0, 1.0),
    # xlim = c(0, 16)
    xlim = c(0, 190)
    ) +
  scale_color_manual(
    breaks = muts_scale_breaks,
    labels = muts_scale_labels,
    values = muts_scale_values
  ) +
  scale_shape_manual(
    breaks = c(1, 2, 3),
    values = c(15, 16, 17)
  ) +
  labs(x = "Time (min)", y = "Fraction bound",
       color = "AGO2 mutant",
       shape = "Experiment replicate") +
  theme0

df.fb %>%
  filter(time > 0) %>%
  group_by(miR, AGOmut) %>%
  count()

# Plot koff ----
df.koff = df.fit %>%
  filter(model.term == "koff") %>%
  transmute(miR = miR,
            AGOmut = AGOmut,
            miR_AGOmut = paste0(miR, " ", AGOmut),
            koff = exp(model.estimate),
            koff.log = model.estimate,
            koff.log.se = model.std.error,
            koff.lo = exp(model.estimate-1.96*model.std.error),
            koff.hi = exp(model.estimate+1.96*model.std.error),
            Toff = log(2)/koff,
            Toff.s = Toff * 60
  ) %>%
  mutate(miR_AGOmut = factor(miR_AGOmut, levels = c(
    "miR-200b WT", "miR-200b EI_5xE",
    "let-7a WT", "let-7a EI_5xE"
  )))

df.koff %>%
  mutate(miR_AGOmut = fct_rev(miR_AGOmut)) %>%
  ggplot(aes(x = koff, y = miR_AGOmut, color = AGOmut)) +
  geom_errorbarh(aes(xmin = koff.lo, xmax = koff.hi),
                 height = 0.4, linewidth = 0.75, color = "black") +
  geom_point(shape = 18, size = 4, show.legend = F) +
  scale_color_manual(
    breaks = muts_scale_breaks,
    labels = muts_scale_labels,
    values = muts_scale_values
  ) +
  scale_x_continuous(trans = "log10", expand = c(0, 0)) +
  labs(x = "koff (min-1)", y = "AGO2 mutant") +
  coord_cartesian(xlim = c(
    0.002, 4
    )) +
  theme0
