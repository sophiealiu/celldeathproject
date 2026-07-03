# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/29/2026 11:34am CDT
# Purpose: Nice lil spatial overlay graphics.
# -----------------------------------------------------------------------------

library(dplyr)
library(gghalves)
library(ggplot2)
library(pheatmap)
library(tidyr)
library(viridis)

# importing files. same µm coordinate space
datadir <- "path/to/your/working/directory"
df_isoIF <- read_csv(file.path(datadir,"iso_coords.csv"))
df_iso_vis <- read_csv(file.path(datadir,"NMF_isoC.csv"))

df_pd1IF <- read_csv(file.path(datadir,"pd1_coords.csv"))
df_pd1_vis <- read_csv(file.path(datadir,"NMF_Apd1.csv"))


# -----------------------------------------------------------------------------
# 1. regular overlay factors on top of IF
draw <- function(df_IF, df_visium,
                 factorcols,
                 colors, size, trans,
                 threshold) { 
  
  df_visium_long <- df_visium %>%                # flipping to orient
    pivot_longer(
      cols = factorcols,
      names_to = "factor",
      values_to = "value") %>%    
      group_by(factor) %>%                       # isolate highest expressed areas
      filter(value > quantile(value, threshold, na.rm = TRUE)) %>%
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
      panel.grid       = element_blank(),                                  # later gridlines
      
      axis.text.x  = element_text(color = "white", size = 12),
      axis.text.y  = element_text(color = "white", size = 12),
      axis.line = element_line(color = "white")
      
      plot.title = element_text(color = "white", size = 15, face = "bold"),
    )  +
    
    scale_color_manual(
      values = colors,
      name = "factor",
      guide = guide_legend(
        override.aes = list(shape = 15, size = 5, alpha = 1))) 
}

# now employing
colors <- c("#9678B6")                                            # purple mountain majesty, courtesy of DB
factorcols <- c("factor1")

test <- draw(df_isoIF, df_iso_vis, factorcols, colors, 1, 1, 0.6)
test


# -----------------------------------------------------------------------------
# 2. summary graphs for experimental condition. repeat for factors of interest
# gene heats
df_weights <- read.csv(file.path(datadir, "NMF_W.csv"))

list_weights <- list()
fact_cols <- colnames(df_weights)[-1]

for (i in fact_cols) {
  top10 <-df_weights[
    order(df_weights[[i]], decreasing = TRUE),
    1][1:10]
  
  list_weights <- c(list_weights, top10)
}

clean_weights <- unique(list_weights)

subset_mat <- df_weights[df_weights[,1] %in% clean_weights, -1]
rownames(subset_mat) <- df_weights[df_weights[,1] %in% clean_weights, 1]

pheatmap(subset_mat,
         scale = "row", color = viridis(50),
         main = "Top gene expression in each factor (z-scored)",
         cluster_rows = FALSE,   
         cluster_cols = FALSE,  
         cellwidth = 20,         
         cellheight = 10)


# -----------------------------------------------------------------------------
# factor weights
ggplot(df_pd1_vis, aes(x = x, y = y, z = factor1)) +
  stat_summary_2d(fun = mean, binwidth = c(300, 300),
                  alpha = 0.9) +      
  scale_fill_viridis_c(option = "magma") +
  
  # geom_point(
  #   data  = filter(df_pd1IF, cell_type == "gc3ai"),
  #   aes(x = x, y = y),
  #   inherit.aes = FALSE,
  #   color = "green4",
  #   size  = 0.1,
  #   alpha = 1
  # ) +
  geom_point(
    data  = filter(df_pd1IF, cell_type == "lectin"),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "#D9C20B",
    alpha = 0.1,
    size  = 0.005) +
  
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "black"),
    axis.title       = element_blank(),              
    axis.line        = element_blank(),
    axis.text        = element_blank()
  ) +
  
  scale_x_continuous(breaks = seq(0, 10000, by = 2000)) +
  scale_y_continuous(breaks = seq(0, 10000, by = 2000)) +
  
  annotate("segment", x = 1000, xend = 2000, y = 2000, yend = 2000,        
           colour = "white", linewidth = 3)


# ------------------------------------------------------------------------------
# 3. being *extra* and automating frame creation, joining to make movie in Adobe
for (i in 0:30) { 
  alpha <- 0.005*i 
  p <- draw(df_pd1IF, df_pd1_vis, alpha)                      # by transparency intervals
    
  ggsave(
    filename <- sprintf("frames/frame_%02d_alpha_%0.1f.png", i, alpha),
    plot = p,
    width = 6,
    height = 6,
    dpi = 300
  )
}

