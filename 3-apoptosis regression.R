# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/15/2026 4:02pm MDT
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(DHARMa)       # some packages will have to be forced in using remotes::install_github
library(dplyr)
library(ggplot2)
library(MASS)
library(mgcv)
library(tidyr)

datadir <- "path/to/your/working/directory"
iso <- read.csv(file.path(datadir,"NMF_iso40.csv"))       # repeat with varying radii for sensitivity analysis,
pd1 <- read.csv(file.path(datadir,"NMF_pd140.csv"))       # which involves regenerating from file 2
df <- rbind(iso, pd1)


# -----------------------------------------------------------------------------
# BEGINNING MODELING
# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - NMF additive non-orthogonal nature. raw counts, per-gene weights
# - pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
fact_cols <- grep("^factor", names(df), value = TRUE)   
x <- as.matrix(df[, fact_cols])         
x_scaled <- scale(x)                    # visible representation of RNA counts. z-scores
                
trt <- ifelse(df$sample == "apd1", 1, 0)


# -----------------------------------------------------------------------------
# 2. using counts considers sequencing depth. 
# distribution explains family choice
hist(df$n_dying, 
     breaks = 30,
     main = "Distribution of count dying",
     xlab = "number of dying cells in 40 micron vicinity")

# a. negative binomial.
nb <- glm.nb(
  n_dying ~ n_immune + trt + x_scaled + offset(log(n_tumor)),
  data = df
)

# b. geographical trends, Moran's I. Not combined w/ zero. 
# excess overfit possible with smoothing taking credit.
spatial <- gam(
  n_dying ~ n_immune + trt + x_scaled +
    s(x, y, bs = "gp", k = 30) +        # depends on power of gaussian process smoothing
    offset(log(n_tumor)),
  data = df,
  family = nb()
)

# parametric p-vals adjusted
p_para <- summary(spatial)$p.pv
p_adj <- p.adjust(p_para, method = "BH")

# c. treatment interactions
nb_int <- glm.nb(
  n_dying ~ n_immune + trt * x_scaled + offset(log(n_tumor)),
  data = df
)


# -----------------------------------------------------------------------------
# 3. diagnostic stats
sim_res <- simulateResiduals(nb)
plot(sim_res)

testZeroInflation(sim_res)
testSpatialAutocorrelation(simulationOutput = sim_res, 
                           x = df$x, 
                           y = df$y)


# -----------------------------------------------------------------------------
# 4. permutation test, random noise effects
set.seed(42)
df_perm <- df
df_perm$n_dying <- sample(df_perm$n_dying)

nb_mess <- glm.nb(
  n_dying ~ n_immune + trt + x_scaled + offset(log(n_tumor)),
  data = df_perm
)

