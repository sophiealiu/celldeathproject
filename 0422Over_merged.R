# this file is visualization to confirm data. Can generate a nice graphic

library(readr)
library(dplyr)
library(ggplot2)

inputdir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"


df_IF <- read_csv(file.path(inputdir, "iso7_coords_clean.csv"))
df_visium <- read_csv(file.path(inputdir, "iso7spatial_sig.csv"))


# *****************************************************************************
# creating plots **************************************************************
check <- function(df_IF, df_visium, trans) { 
  ggplot() +
    # IF Imaris DATA FIRST ********************************************************
    stat_density_2d(
      data = filter(df_IF, cell_type == "tdtomato"), # used density bc too many tdtomato
      aes(x = x, y = y),
      fill        = "red",
      geom        = "density_2d_filled",
      contour_var = "ndensity",
      alpha       = 0.10,
      n           = 400,
      na.rm       = TRUE
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "cd8"),    # but points ok for others
      aes(x = x, y = y),
      color = "#02819e",
      alpha = 0.20,
      size  = 0.01
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "gc3ai"),
      aes(x = x, y = y),
      color = "green4",
      size  = 0.0005, 
      alpha = 0.15
    ) +  
    
    # overlay plotting ************************************************************
    geom_point(
      data = filter(df_visium, gene %in% c(cd8_effector, 
                                           monocytes, 
                                           tregs,
                                           cdc1,
                                           fibroblasts)),
      aes(x = x, y = y),
      alpha = trans,
      size = 0.01
    ) +
    
    # formatting so it looks better ***********************************************
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
      values = c(
        "cd8_effector" = "white",
        "monocytes" = "purple",
        "tregs" = "yellow",
        "cdc1" = "pink",
        "fibroblast" = "orange"  # add in others but for color/visualization purposes
      ),
      name = "gene signature",
      guide = guide_legend(
        override.aes = list(shape = 15, size = 5, alpha = 1)
    )
  
}


# *****************************************************************************
# automating frame creation, join to make movie in Adobe **********************

see <- draw(df_iso7IF, score_over, 0.05)
see # make sure it looks okay before running below


setwd(
  "I:/Hu Lab/Sophie/Visium/visium image manual spot selection/20260413_final_merge/current data/graphics creation")
                                      # print frames up to alpha = 0.15, should be bright enough
for (i in 0:30) { 
  alpha <- 0.005*i 
  p <- draw(df_IF, df_visium, alpha) # by transparency intervals
    
  ggsave(
    filename <- sprintf("frames/frame_%02d_alpha_%0.1f.png", i, alpha),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}
