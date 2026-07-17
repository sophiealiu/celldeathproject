# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(DHARMa)       # some packages will have to be forced in using remotes::install_github
library(dplyr)
library(ggplot2)
library(MASS)
library(mgcv)
library(tidyr)

datadir <- "path/to/your/working/directory"
df_iso <- read.csv(file.path(datadir,"NMF_iso40.csv"))         # repeat with varying radii for sensitivity analysis,
df_pd1 <- read.csv(file.path(datadir,"NMF_pd140.csv"))         # which involves regenerating from file 2
df_merged <- rbind(df_iso, df_pd1)


# -----------------------------------------------------------------------------
# BEGINNING MODELING
# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - NMF additive non-orthogonal nature. raw counts, per-gene weights
# - pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
fact_cols <- grep("^factor", names(df_merged), value = TRUE)   
fact <- as.matrix(df_merged[, fact_cols])         
fact_sc <- scale(fact)                      # visible representation of RNA counts. z-scores
                
trt <- ifelse(df_merged$sample == "apd1", 1, 0)


# -----------------------------------------------------------------------------
# 2. using counts considers sequencing depth. 
# distribution explains family choice
hist(df_merged$n_dying, 
     breaks = 30,
     main = "Distribution of count dying",
     xlab = "number of dying cells in 40 micron vicinity")

# a. negative binomial family model, number dying was our response variable.
nb <- glm.nb(
  n_dying ~ n_immune + trt + fact_sc + offset(log(n_tumor)),        # offset accounting for density effects
  data = df_merged
)
summary(nb)

# b. geographical trends, Moran's I to address spatial autocorrelation, cautious overfitting. 
# excess overfit possible with smoothing taking credit.
spatial <- gam(
  n_dying ~ n_immune + trt + fact_sc +
    s(x, y, bs = "gp", k = 30) +            # smoothing using coordinates depends on power of gaussian process
    offset(log(n_tumor)),
  data = df_merged,
  family = nb()
)

# parametric p-vals adjusted
p_para <- summary(spatial)$p.pv                # extracts p-vals from spatial model
p_adj <- p.adjust(p_para, method = "BH")       # benjamini-hochberg correction

# c. treatment interactions
nb_int <- glm.nb(
  n_dying ~ n_immune + trt * fact_sc + offset(log(n_tumor)),
  data = df_merged
)


# -----------------------------------------------------------------------------
# 3. diagnostic stats
sim_res <- simulateResiduals(nb)
plot(sim_res)

testZeroInflation(sim_res)
testSpatialAutocorrelation(simulationOutput = sim_res, 
                           x = df_merged$x, 
                           y = df_merged$y)


# -----------------------------------------------------------------------------
# 4. permutation test, random noise effects
set.seed(42)
df_merged_perm <- df_merged
df_merged_perm$n_dying <- sample(df_merged_perm$n_dying)

nb_mess <- glm.nb(
  n_dying ~ n_immune + trt + fact_sc + offset(log(n_tumor)),
  data = df_merged_perm
)

