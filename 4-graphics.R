# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/08/2026
# Purpose: Creating nice spatial overlay graphics.
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)


inputdir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"

# reading in my dataframes
df_iso7IF <- read_csv(file.path(inputdir, "iso7_coords_clean.csv"))
df_iso7visium <- read_csv(file.path(inputdir, "iso7spatial_sig.csv"))

df_pd1IF <- read_csv(file.path(inputdir, "pd1-9_coords_final.csv"))
df_pd1visium <- read_csv(file.path(inputdir, "pd1spatial_sig.csv"))


# -----------------------------------------------------------------------------
draw <- function(df_IF, df_visium, trans) { 
  df_visium_long <- df_visium %>% # have to flip the format so it recognizes gene cols
    pivot_longer(
      cols = c(
        # X1,
        #X2
        X3
        # X4,
        # X5,
        #X6
        #X7,
        #X8
        # X9,
        # X10
        #cdc1, cdc2
      ),
      names_to = "factor",
      values_to = "value") %>%
    group_by(factor) %>%
    filter(value > quantile(value, 0.8)) %>%        # isolate highest expressed areas
    ungroup()
    
  ggplot() +
  
    stat_density_2d(
      data = filter(df_IF, cell_type == "tdtomato"), # density avoiding overplotting
      aes(x = x, y = y),
      fill        = "red",
      geom        = "density_2d_filled",
      contour_var = "ndensity",
      alpha       = 0.10,
      n           = 400,
      na.rm       = TRUE
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "cd8"),    # points for clarity
      aes(x = x, y = y),
      color = "#02819e",
      alpha = 0.20,
      size  = 0.01
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "gc3ai"),
      aes(x = x, y = y),
      color = "green4",
      size  = 0.1, 
      alpha = 0.25
    ) +  
    geom_point(
      data  = filter(df_IF, cell_type == "lectin"),    
      aes(x = x, y = y),
      color = "#D9C20B",
      alpha = 0.1,
      size  = 0.005) +
    
    geom_point(
      data = filter(df_visium_long),                # overlaying
      aes(
        x = cx,                                     # renamed vars in jupyter
        y = cy,
        color = factor
      ),
      size = 0.1
    ) +
    
# formatting
    scale_x_continuous(expand = c(0.05, 0.05)) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    coord_equal(clip = "off") +
    theme(
      panel.background = element_rect(fill = "black"),
      plot.background  = element_rect(fill = "black"),
      panel.grid       = element_blank(),
      
      axis.text.x  = element_text(color = "white", size = 12),
      axis.text.y  = element_text(color = "white", size = 12),
      
      plot.title = element_text(color = "white", size = 15, face = "bold"),
      
      axis.line = element_line(color = "white")
    )  +
    
    scale_color_manual(
      values =  c(
        # "X1" = "",
        #"X2" = color2
        "X3" = color2
        # "X4" = "",
        # "X5" = "",
        #"X6" = color1
        # "X7" = "#9678B6",                         # purple mountain majesty, courtesy of DB
       #"X8" = color1           
        # "X9" = "",
        # "X10" = ""  
        #cdc1 = color1, cdc2 = color2
      ),
      name = "gene signature",
      guide = guide_legend(
        override.aes = list(shape = 15, size = 5, alpha = 1)))
}

# -----------------------------------------------------------------------------
# Visualization & exporting for more graphics
test <- draw(df_iso7IF, df_iso7visium, 0.2)
test

setwd(
  "I:/Hu Lab/Sophie/Visium/visium image manual spot selection/20260413_final_merge/current data/graphics creation")

# automating frame creation, join to make movie in Adobe
for (i in 0:30) { 
  alpha <- 0.005*i 
  p <- draw(df_IF, df_visium, alpha)               # by transparency intervals
    
  ggsave(
    filename <- sprintf("frames/frame_%02d_alpha_%0.1f.png", i, alpha),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}
