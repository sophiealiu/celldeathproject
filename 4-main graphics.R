# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Purpose: Nice lil spatial overlay graphics.
# -----------------------------------------------------------------------------

library(dplyr)
library(gghalves)
library(ggplot2)
library(pheatmap)
library(scales)
library(tidyr)
library(viridis)

# importing files
# datadir <- "path/to/your/working/directory"
datadir <- "I:/Hu Lab/Sophie/1. Cell death/all final data"

df_isoIF <- read.csv(file.path(datadir,"iso_coords.csv"))
df_iso_vis <- read.csv(file.path(datadir,"NMF_isoC.csv"))

df_pd1IF <- read.csv(file.path(datadir,"apd1_coords.csv"))
df_pd1_vis <- read.csv(file.path(datadir,"NMF_Apd1.csv"))

df_pd1_trans <- df_pd1_vis
df_pd1_trans$x <- df_pd1_vis$x / 1.5454            # scaling factor from Loupe browser 
df_pd1_trans$y <- df_pd1_vis$y / 1.5454            # to align in the same µm coordinate space

df_iso_trans <- df_iso_vis                         # avoiding overwriting
df_iso_trans$x <- df_iso_vis$x / 1.5454
df_iso_trans$y <- df_iso_vis$y / 1.5454


# -----------------------------------------------------------------------------
# 1. function for regular overlay factors on top of IF
# Figure 5, panel B. replace with clusters when calling the function
draw <- function(df_IF, df_visium,
                 factorcols,
                 colors, size, trans,
                 threshold) { 
  
  df_visium_long <- df_visium %>%            # flipping format for R to understand/process
    pivot_longer(                            # not necessary for clusters
      cols = factorcols,
      names_to = "factor",
      values_to = "value") %>%    
      group_by(factor) %>%                   # isolate highest expressed areas
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
    )  +
    
    scale_color_manual(
      values = colors,
      name = "factor",
      guide = guide_legend(
        override.aes = list(shape = 15, size = 5, alpha = 1))) 
}

# now employing. remove factors to "" to get panel with no overlay or comment out from above function
colors <- c("#9678B6", "pink")                    # hex code for purple mountain majesty, courtesy of DB
facts <- c("factor1", "factor2")                  # etc. sub in factors of interest

# recall function takes size, transparency, and threshold percentile (1 is 100th) for top expression
# df_clus <- read.csv(file.path(datadir,"apd1_clusters.csv"))
test <- draw(df_isoIF, df_iso_trans, facts, colors, 1, 1, 0.6)
test


# -----------------------------------------------------------------------------
# 2. summary graphs for experimental condition. repeat for factors of interest
# gene heatmap using internally generated file with gene list for each factor (from block 1)
df_weights <- read.csv(file.path(datadir, "NMF_W.csv"))

list_weights <- list()
fact_cols <- colnames(df_weights)[-1]

for (i in fact_cols) {
  top10 <-df_weights[
    order(df_weights[[i]], decreasing = TRUE),
    1][1:10]          
  
  list_weights <- c(list_weights, top10)
}

# accounting for overlap in genes across factors (Figure 5, panel D)
clean_weights <- unique(list_weights)

subset_mat <- df_weights[df_weights[,1] %in% clean_weights, -1]
rownames(subset_mat) <- df_weights[df_weights[,1] %in% clean_weights, 1]

pheatmap(subset_mat,
         scale = "row", color = viridis(50),
         main = "Top gene expression in each factor (z-scored)",
         cluster_rows = TRUE,   
         cluster_cols = FALSE,  
         cellwidth = 20,         
         cellheight = 10)


# -----------------------------------------------------------------------------
# 3. summary graphs: Figure 5, panel G
total_umisI <- iso_raw$nCount_Spatial.008um
total_umisP <- pd1_raw$nCount_Spatial.008um

names(total_umisI) <- paste0("iso_", names(total_umisI))
names(total_umisP) <- paste0("pd1_", names(total_umisP))

df_iso_vis$total_umis <- total_umisI[df_iso_vis$X]
df_pd1_vis$total_umis <- total_umisP[df_pd1_vis$X]

# change for each factor
pd1_exp9 <- df_pd1_vis$factor9/ df_pd1_vis$total_umis
iso_exp9 <- df_iso_vis$factor9/ df_iso_vis$total_umis

df_temp <- cbind(pd1_exp9, iso_exp9)   # for density
df_plot <- as.data.frame(df_temp)

colnames(df_plot) <- c("PD1", "ISO")
df_long <- pivot_longer(
  df_plot,
  cols = everything(),
  names_to = "condition",
  values_to = "value"
)

 ggplot(df_long, aes(x = "", y = value, fill = condition)) +
    geom_half_violin(
      data = subset(df_long, condition == "ISO"),
      side = "l",
      trim = FALSE
    ) +
    geom_half_violin(
      data = subset(df_long, condition == "PD1"),
      side = "r",
      trim = FALSE
    ) +
    geom_boxplot(
      data = subset(df_long, condition == "ISO"),
      width = 0.05, color = "black", fill = NA, outlier.shape = NA,
      position = position_nudge(x = -0.05) 
    ) + 
    geom_boxplot(
      data = subset(df_long, condition == "PD1"),
      width = 0.05, color = "black", fill = NA, outlier.shape = NA,
      position = position_nudge(x = 0.05) 
    ) + 
    coord_cartesian(ylim = c(0, 3e-8)) + 
    theme_classic()                    
    labs(title = "factor 4 / umi", x = "", y = "expression")


# -----------------------------------------------------------------------------
# 4. plotting location and intensity of factors, accounting for sequencing depth
# Figure 5, panel H. 
# Panel C is from an older dataframe (block 1.1 giving cell type signatures)
# to plot panel C, just replace z with the cell type and the dataframe in the function to df_sig.
df_sig <- read.csv(file.path(datadir, "apd1_sig.csv"))
df_sig_trans <- df_sig
df_sig_trans$x <- df_sig$x / 1.5454            # scaling factor from Loupe browser again
df_sig_trans$y <- df_sig$y / 1.5454            

ggplot(df_pd1_vis, aes(x = x, y = y, z = factor9/total_umis)) +
  stat_summary_2d(binwidth = c(48, 48),
                  alpha = 1) +      
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1e-7)           # alter upper limit for comparison btwn conditions but not along factors
  ) +
  
  geom_point(
    data  = filter(df_pd1IF, cell_type == "gc3ai"),
    aes(x = x, y = y),
    inherit.aes = FALSE,
    color = "green4",
    size  = 0.5,
    alpha = 0.3
  ) +
  # geom_point(
  #   data  = filter(df_pd1IF, cell_type == "lectin"),
  #   aes(x = x, y = y),
  #   inherit.aes = FALSE,
  #   color = "#D9C20B",
  #   alpha = 0.1,
  #   size  = 0.005) +
  
  theme(
    panel.background = element_rect(fill = "black"),
    plot.background  = element_rect(fill = "black"),
    axis.title       = element_blank(),              
    axis.line        = element_blank(),
    axis.text        = element_blank()
  ) + 
  
  scale_x_continuous(breaks = seq(0, 10000, by = 2*1000)) +
  scale_y_continuous(breaks = seq(0, 10000, by = 2*1000)) +
  
  annotate("segment", x = 1000, xend = 2000, y = 500, yend = 500,        
           colour = "white", linewidth = 3)                   # scaled in our data 1000 = 1mm


# ------------------------------------------------------------------------------
# 5. being *extra* and automating frame creation, joining to make movie in Adobe
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

