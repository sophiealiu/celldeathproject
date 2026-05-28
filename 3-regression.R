# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/18/2026
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(dplyr)
library(glmnet)
library(pheatmap)

datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"
iso <- read.csv(file = file.path(datadir, "0518_NMF_iso_calc.csv"))     
pd1 <- read.csv(file = file.path(datadir, "0518_NMF_pd1_calc.csv"))     


# -----------------------------------------------------------------------------
# 0. primary sanity checks of known cell type and my NMF factors
# nothing specific found with cd8s. but yes for tdtomato.
fact_cols <- names(iso)[12:21]                         # indices, NMF columns begin @12
cell_types <- unique(iso$cell_type)
factor_sums <- colSums(iso[, fact_cols])

# Distribution of cell across factors
percent_type <- list()
for (type in cell_types) {
  subset_df <- iso[iso$cell_type == type, ]
  
  by_factor <- colSums(subset_df[, fact_cols])
  total <- sum(by_factor)
  
  percent <- (by_factor/total) * 100
  
  percent_type[[type]] <- percent
}
percent_type <- do.call(rbind, percent_type)
percent_type

# visualization
mat <- as.matrix(percent_type)

pheatmap(mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

# top 20 associated genes with factors
head(sort(W[,7], decreasing = TRUE), 20)


# ----------------------------------------------------------------------------_
# BEGINNING MODELING

# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights


# -----------------------------------------------------------------------------
# 1. Cleaning, removing sparsity effects (0,1) instead of [0,1] inclusive
df_nonempty <- pd1 %>%
  filter(prop_dying != 0 & prop_dying != 1)

x <- as.matrix(pd1[, fact_cols])
x_scaled <- scale(as.matrix(df_nonempty[, fact_cols]))


# -----------------------------------------------------------------------------
# 2. Beta regression. determined from 0515 that this was the best model
# determining significance from experimental condition because that's where I trained NMF
library(betareg)

# y <- pd1$exist_dying                  # binomial regression [N] from inspection
y_cont <- df_nonempty$prop_dying        # linear gaussian     [N]
                                        # beta regression     [possible]
y_disc <- df_nonempty$n_dying           # poisson             [N], but
                                        # negative binomial   [possible]
# - may require transformation of response variable addressing variance
# - poisson dispersion ~ 2. residuals fan out


beta <- betareg(y_cont ~ x_scaled)


# -----------------------------------------------------------------------------
# 2. Negative binomial regression for counts discrete data, 
# maintains more uncertainty
library(MASS)

# requires more incorporation of polynomial or spline terms. 
# however those still maintain fanned residuals
nb <- glm.nb(y_disc ~ x_scaled)
nb_poly <- glm.nb(y_disc ~ poly(x_scaled,2))

library(mgcv)
nb_spline <- gam(y_disc ~ s(x_scaled),  # larger k = more "curve" in model
                 family = nb(), 
                 data = overall)


# -----------------------------------------------------------------------------
# 4. summary stats, NB shows worse residuals, therefore maybe not the best model
summary(nb)
summary(beta)

plot(beta)
qqnorm(residuals(beta, type = "quantile"))
qqline(residuals(beta, type = "quantile"))

# displays different default plots
plot(nb)
plot(cooks.distance(nb), type = "h", 
     ylab = "Cook's Distance", xlab = "Obs. Number")


# -----------------------------------------------------------------------------
# 4. isolating high-death chunk vs low-death (1mm^2), dying pocket (0.25mm^2)
# issue is too sparse!
df_high_death <- df_nonempty %>%
  filter(cx >= 2500 & cx <= 3500,      # recall I renamed the variables in jupyter
         cy >= 4500 & cy <= 5500)

df_low_death <- df_nonempty %>%
  filter(cx >= 5000 & cx <= 6000,
         cy >= 4500 & cy <= 5500)

df_death_pocket <- df_nonempty %>%
  filter(cx >= 3500 & cx <= 4500,
         cy >= 3000 & cy <= 4000)

x_spec <- as.matrix(df_death_pocket[, fact_cols])
y_spec <- df_death_pocket$n_dying
nb_spec <- glm.nb(y_spec ~ x_spec)

