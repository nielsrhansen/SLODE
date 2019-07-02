## Theme for ggplot ##

library(ggplot2)

angle_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 17, face = "bold"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top")

ggplot_theme <- theme_bw() + 
  theme(axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 17, face = "bold"),
    strip.text.x = element_text(size = 15),
    strip.text.y = element_text(size = 15),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.title = element_text(hjust = 0.5),
    legend.position = "top")


## Color scheme
col_scheme <- c(AIM = "#0072B2", EGM = "#CC79A7", IM = "#009E73", SCAD = "#E69F00")
