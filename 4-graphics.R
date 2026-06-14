# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/29/2026 11:34am CDT
# Purpose: Nice lil spatial overlay graphics.
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)

# importing files
df_isoIF <- read_csv("iso_coordsNT.csv")
df_iso_vis <- read_csv("9NMF_isoC.csv")

df_pd1IF <- read_csv("pd1_coordsNT.csv")
df_pd1_vis <- read_csv("9NMF_Apd1.csv")


# -----------------------------------------------------------------------------
# 1. regular overlay factors on top of IF
draw <- function(df_IF, df_visium,
                 factor_cols,
                 colors, size, trans,
                 threshold) { 
  
  df_visium_long <- df_visium %>%                         # flipping to orient
    pivot_longer(
      cols = factor_cols,
      names_to = "factor",
      values_to = "value") %>%
    group_by(factor) %>%
    filter(value > quantile(value, threshold)) %>%        # isolate highest expressed areas
    ungroup()
  
  ggplot() +
    
    stat_density_2d(
      data = filter(df_IF, cell_type == "tdtomato"), 
      aes(x = x, y = y),
      fill        = "red",
      geom        = "density_2d_filled",
      contour_var = "ndensity",
      alpha       = 0.10,
      n           = 400,
      na.rm       = TRUE
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "cd8"),    
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
      data = df_visium_long,     
      aes(
        x = x,                                           
        y = y,
        color = factor
      ),
      size = size,
      alpha = trans
    ) +
    
    # formatting
    scale_x_continuous(expand = c(0.05, 0.05)) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    coord_equal(clip = "off") +
    theme(
      panel.background = element_rect(fill = "black"),
      plot.background  = element_rect(fill = "black"),
      panel.grid       = element_blank(),                 # change later if I want gridlines 
      
      axis.text.x  = element_text(color = "white", size = 12),
      axis.text.y  = element_text(color = "white", size = 12),
      
      plot.title = element_text(color = "white", size = 15, face = "bold"),
      
      axis.line = element_line(color = "white")
    )  +
    
    scale_color_manual(
      values = colors,
      name = "factor",
      guide = guide_legend(
        override.aes = list(shape = 15, size = 5, alpha = 1)))
}


# -----------------------------------------------------------------------------
# 2. dead or cd8 spots on top of heatmap of factors
heat <- function(df_IF, df_visium, 
                 factor_num, type, 
                 color, pt_size, trans, 
                 bin_num, thresh_fact) {
  
   df_vis_filt <- df_visium %>% 
     filter(.data[[factor_num]] > quantile(.data[[factor_num]], thresh_fact))
   
   ggplot() +
  
   geom_hex(data = df_vis_filt, aes(x = x, y = y),
            bins = bin_num,
            alpha = trans) +
   scale_fill_viridis_c(option = "A") +
   
   geom_point(
     data  = filter(df_IF, cell_type == type),
     aes(x = x, y = y),
      color = color,
    size  = pt_size
  ) + 
     
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
     )
}


# -----------------------------------------------------------------------------
# 3. actually happening now
colors <- c("#9678B6")                                   # purple mountain majesty, courtesy of DB
factor_cols <- c("X11")

test <- draw(df_isoIF, df_isovisium, factor_cols, colors, 1, 1, 0.6)
test

#D9C20B
test2 <- heat(df_isoIF, df_isovisium, "X17", "lectin",
              "#D9C20B", 3, 0.6, 30, 0)                 
test2


# ------------------------------------------------------------------------------
# being *extra* and automating frame creation, joining to make movie in Adobe
setwd(
  "I:/Hu Lab/Sophie/Visium/visium image manual spot selection/20260413_final_merge/current data/graphics creation")

for (i in 0:30) { 
  alpha <- 0.005*i 
  p <- draw(df_IF, df_visium, alpha)                      # by transparency intervals
    
  ggsave(
    filename <- sprintf("frames/frame_%02d_alpha_%0.1f.png", i, alpha),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}

