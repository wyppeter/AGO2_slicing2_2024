# Code to plot major-groove widths
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

Aform = 13.4
vdW = 5.8

grooves = read.csv("refined_major_grooves_Oct.csv") %>%
  pivot_longer(cols = -position, names_to = "model.name", values_to = "groove.width") %>%
  filter(!is.na(groove.width)) %>%
  mutate(model.name = sub("^x", "", model.name),
         groove.width = groove.width - vdW)  # subtract vdW

grooves %>%
  ggplot(aes(x = position, y = groove.width-(Aform-vdW), color = model.name)) +
  geom_hline(yintercept = 0, linewidth = 0.75, color = "black", linetype = "dashed") +
  geom_line(linewidth = 0.75) +
  geom_point(shape = 16, size = 2) + 
  scale_x_continuous(breaks = seq(0, 24, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-20, 20, 2), expand = c(0, 0)) +
  scale_color_manual(
    breaks = c(
      "7SWF",
      "9CMP"
    ),
    values = c(
      "hotpink1",
      "purple4"
  )) +
  coord_cartesian(xlim = c(4.2, 19.8), ylim = c(-7, 7)) +
  labs(x = "Guide/target position", y = "Major groove width (A), relative to standard A-form") +
  theme0
