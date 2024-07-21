# Slicing kinetics plots
# Bartel Lab
# Peter Y. Wang 2024

library(tidyverse)

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
theme.ks = theme(
  axis.text.y = element_text(color = "black", size = 12, hjust = 0)
)

df.k = read.csv("slicing_ks_p2_raw.csv")

#### CHOOSE ONE
rxnLEVELS = c("miR-7-L22#WT_perfect",
              "miR-7-L22#H56A_perfect",
              "miR-7-L22#R68A_perfect",
              "miR-7-L22#R97A_perfect",
              "miR-7-L22#R97A.K98A_perfect",
              "miR-7-L22#R97E.K98E_perfect",
              "miR-7-L22#quadA_perfect",
              "miR-7-L22#WT_16bp",
              "miR-7-L22#H56A_16bp",
              "miR-7-L22#R68A_16bp",
              "miR-7-L22#R97A.K98A_16bp",
              "miR-7-L22#R97E.K98E_16bp",
              "miR-7-L22#quadA_16bp")
# rxnLEVELS = c("miR-7-L22#WT_perfect", "miR-7-L22#R710A_perfect", "miR-7-L22#H712A_perfect", "miR-7-L22#R710A.H712A_perfect")
# rxnLEVELS = c("miR-7-L22#WT_mm10", "miR-7-L22#R710A_mm10", "miR-7-L22#H712A_mm10", "miR-7-L22#R710A.H712A_mm10",
#               "miR-7-L22#WT_mm11", "miR-7-L22#R710A_mm11", "miR-7-L22#H712A_mm11", "miR-7-L22#R710A.H712A_mm11")

# kslice (ranges) ----
df.k %>%
  filter(rxn %in% rxnLEVELS) %>% mutate(rxn = factor(rxn, levels = rev(rxnLEVELS))) %>%
  ggplot(aes(x = kcat, y = rxn)) +
  geom_errorbarh(aes(xmin = kcat.lo, xmax = kcat.hi),
                 color = "black", linewidth = 0.75, height = 0.4) +
  geom_point(shape = 18, size = 4, show.legend = F, color = "black") +
  scale_x_continuous(trans = "log10", 
                     breaks = 10^seq(-10,10,1),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(
    0.001, 10
    # 0.012, 6
    # 0.00005, 0.025
    )) +
  labs(x = "k_slice (min-1)", y = "") +
  theme0 +
  theme.ks

#################

# mm10/11 effects
cent.mm.expes = c(
  "miR-7-L22#WT_perfect", "miR-7-L22#R710A_perfect", "miR-7-L22#H712A_perfect", "miR-7-L22#R710A.H712A_perfect",
  "miR-7-L22#WT_mm10", "miR-7-L22#R710A_mm10", "miR-7-L22#H712A_mm10", "miR-7-L22#R710A.H712A_mm10",
  "miR-7-L22#WT_mm11", "miR-7-L22#R710A_mm11", "miR-7-L22#H712A_mm11", "miR-7-L22#R710A.H712A_mm11"
)
df.centmm = df.k %>%
  filter(rxn %in% cent.mm.expes) %>%
  separate(col = "rxn", sep = "#|_", into = c("miR", "AGOmut", "targ")) %>%
  mutate(kcat.se = (log(kcat.hi)-log(kcat))/1.96)

df.centmm.perfect = df.centmm %>%
  filter(targ == "perfect") %>%
  transmute(AGOmut = AGOmut, kcat.perfect = kcat, kcat.se.perfect = kcat.se)

df.centmm.delta = df.centmm %>%
  filter(targ != "perfect") %>%
  left_join(df.centmm.perfect, by = "AGOmut") %>%
  mutate(
    kcat.FC = kcat.perfect/kcat,
    kcat.se.comp = sqrt(kcat.se^2+kcat.se.perfect^2),
    kcat.FC.hi = kcat.FC * exp(kcat.se.comp*1.96),
    kcat.FC.lo = kcat.FC / exp(kcat.se.comp*1.96)
  )

# df.centmm.delta %>%
#   mutate(AGOmut = factor(AGOmut, levels = c("WT", "R710A", "H712A", "R710A.H712A"))) %>%
#   ggplot(aes(x = AGOmut, y = kcat.FC)) +
#   geom_col(width = 0.75, fill = "gray50") +
#   geom_errorbar(aes(ymin = kcat.FC.lo, ymax = kcat.FC.hi),
#                 linewidth = 0.75, width = 0.3) +
#   geom_hline(yintercept = 1, color = "black", linewidth = 0.75) +
#   scale_y_continuous(trans = "log10") +
#   coord_cartesian(ylim = c(1, 10000)) +
#   facet_wrap("targ") +
#   labs(x = "AGO variant",
#        y = "Fold reduction in kslice upon target mismatch") +
#   theme0 +
#   theme(axis.text.x = element_text(color = "black", size = 12, hjust = 0, angle = -45))
#   

options(scipen = 5)
df.centmm.delta %>%
  mutate(AGOmut = factor(AGOmut, levels = rev(c("WT", "R710A", "H712A", "R710A.H712A")))) %>%
  ggplot(aes(y = AGOmut, x = 1/kcat.FC)) +
  geom_errorbarh(aes(xmin = 1/kcat.FC.lo, xmax = 1/kcat.FC.hi),
                 color = "black", linewidth = 0.75, height = 0.4) +
  geom_point(shape = 18, size = 4, show.legend = F, color = "black") +
  geom_hline(yintercept = -Inf, color = "black", linewidth = 0.75) +
  scale_x_continuous(trans = "log10") +
  coord_cartesian(xlim = 10^c(-3.8, -1.5)) +
  facet_wrap("targ", ncol = 1) +
  labs(y = "AGO variant",
       x = "Fold change in kslice upon target mismatch, log10") +
  theme0
