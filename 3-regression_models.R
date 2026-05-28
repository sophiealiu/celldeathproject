# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/25/2026
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(dplyr)
library(glmnet)
library(pheatmap)
library(tidyr)

datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"
iso <- read.csv(file = file.path(datadir, "0525_NMF_iso_calc.csv"))
iso_sub <- read.csv(file = file.path(datadir, "0519_NMF_iso_calc_smallRad.csv"))
                                              # 10 micron vs 40 micron radius. verification

pd <- read.csv(file = file.path(datadir, "0525_NMF_pd1_calc.csv"))     
pd_sub <- read.csv(file = file.path(datadir, "0519_NMF_pd_calc_smallRad.csv"))


# -----------------------------------------------------------------------------
# 0. primary checks of known cell type and my NMF factors
fact_cols <- names(pd)[12:28]                         # indices, NMF columns begin @12
cell_types <- unique(pd$cell_type)
factor_sums <- colSums(pd[, fact_cols])

# Distribution of cell across factors
percent_type <- list()
for (type in cell_types) {
  subset_df <- pd[cell_types == type, ]
  
  by_factor <- colSums(subset_df[, fact_cols])
  total <- sum(by_factor)
  
  percent <- (by_factor/total) * 100
  percent_type[[type]] <- percent
}
percent_type <- do.call(rbind, percent_type)

# visualization
mat <- as.matrix(percent_type)

pheatmap(mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# LISTING ASSOCIATED GENES WITH FACTORS
sorted_genes <- (sort(W[,4], decreasing = TRUE))
head(names(sorted_genes), 40)

# ----------------------------------------------------------------------------_
# BEGINNING MODELING

# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights
# - attempted to restore pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. Cleaning, removing sparsity effects/edge effects
pd <- pd %>% drop_na()
iso <- iso %>% drop_na()

x <- as.matrix(iso[, fact_cols])
x_scaled <- scale(x)

# y_binar <- pd$exist_dying 
y_disc <- iso$n_dying
y_cont <- iso$prop_dying   # regions with no dying tumor are still useful

# -----------------------------------------------------------------------------
# 2. Determining the best model
# looking at distribution of data to inform family choice
hist(y_cont, breaks = 30)

# basic log-norm. bad residuals
ln <- glm(y_cont ~ x_scaled,
          family = "gaussian")

r_sq <- 1 - (ln$deviance / ln$null.deviance)

# negative binomial
library(MASS)
nb <- glm.nb(y_disc ~ x_scaled)

# zero-inflated NB
library(glmmTMB)
zero <- glmmTMB(y_disc ~ x_scaled, 
                ziformula = ~1,
                data = pd, 
                family = nbinom2)    # does variance increase linearly? nbinom1

# using dx plots to alter my 10 vars. will be hard to do non-manually later
library(splines)
zero2 <- glmmTMB(y_disc ~ X1 + 
                  log1p(X2 + 1) +    # adding transform for upwards slope
                  ns(X3, df = 3) +   # adding a spline
                  X4 + X5 + X6 + X7 + X8 + X9 +
                  X10 + I(X10^2),    # negative parabolic curve
                
                ziformula = ~1,
                data = pd, 
                family = nbinom2)

# -----------------------------------------------------------------------------
# 4. summary graphics
library(performance)
check_zeroinflation(zero)

library(DHARMa)
sim_res <- simulateResiduals(nb)
plot(sim_res)               # QQ and residual plots

# finding the culprit throwing off my data
for(i in 1:10) {
  col_name <- paste0("X", i)
  plotResiduals(sim_res, form = pd[[col_name]])
  
  # visuals are slow.
  readline(prompt = paste("factor", col_name, "- [Enter] for next"))
}


# -----------------------------------------------------------------------------
# 5. permutation test. Want to see the model break
y_messUp <- sample(y_disc)

nb_mess <- glm.nb(y_messUp ~ x_scaled)

zero_mess <- zeroinfl(y_messUp ~ x_scaled | x_scaled, 
                      data = pd, dist = "negbin")

