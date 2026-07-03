# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/15/2026 4:02pm MDT
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(DHARMa)       # some libraries will have to be forced in using remotes::install_github
library(dplyr)
library(ggplot2)
library(MASS)
library(mgcv)
library(tidyr)

datadir <- "path/to/your/working/directory"

iso <- read.csv(file.path(datadir,"NMF_iso40.csv"))   # repeat with varying radii for sensitivity analysis,
pd1 <- read.csv(file.path(datadir,"NMF_pd140.csv      # involves regenerating from file 2
df <- rbind(iso, pd1)


# -----------------------------------------------------------------------------
# preliminary checks of NMF factors
fact_cols <- names(df)[10:18]                         # indices, check head first

# summary of top genes
W <- read.csv(file.path(datadir, "NMF_W.csv"),row.names = 1)

top_idx <- order(W$V8, decreasing = TRUE)
rownames(W)[top_idx][1:10]


# -----------------------------------------------------------------------------
# BEGINNING MODELING

# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - NMF additive non-orthogonal nature. raw counts, per-gene weights
# - pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
x <- as.matrix(df[, fact_cols])         # replace with cell cols or gene cols based on which block 1
x_scaled <- scale(x)                    # visible representation of RNA counts. z-scores

y_binar <- df$exist_dying               # binomial  [X]
y_disc <- df$n_dying                    # overdispersion -> negative binomial
                                        # zero inflation overfit from performance stats
trt <- ifelse(df$sample == "pd1-9", 1, 0)


# -----------------------------------------------------------------------------
# 3. basic comparative univariate odds-ratios after moving to block 1.2
OR <- glm(
  reformulate(df$n_immune + trt + x_scaled + offset(log(df$n_tumor)), 
              response = y_binar),
  family = "binomial"
)

OR_clean <- tidy(OR, conf.int = TRUE, exponentiate = TRUE)

ggplot(OR_clean %>% 
         aes(x = reorder(term, estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) 


# -----------------------------------------------------------------------------
# 4. using counts considers sequencing depth. iterate after diagnosis in 5
# distribution explains family choice
hist(y_disc, 
     breaks = 30,
     main = "Distribution of count dying",
     xlab = "number of dying cells in 40 micron vicinity")

# a. negative binomial.
nb <- glm.nb(y_disc ~ df$n_immune + trt + x_scaled + offset(log(df$n_tumor)))

# b. geographical trends, Moran's I. Not combined w/ zero. AIC = 4340.9, BIC increase 100
# excess overfit possible with smoothing taking credit.
spatial <- gam(y_disc ~ n_immune + trt + x_scaled + 
                 s(x, y, bs = "gp", k = 10),      # gaussian process term
                 offset = log(df$n_tumor),
                 data = df,                         
                 family = nb())

# parametric p-vals adjusted
p_para <- summary(spatial)$p.pv
p_adj <- p.adjust(p_para, method = "BH")


# -----------------------------------------------------------------------------
# 5. diagnostic stats
sim_res <- simulateResiduals(nb)
plot(sim_res)

# testing if I need an altered NB model (spatial or zero-inflated)
testZeroInflation(sim_res)
testSpatialAutocorrelation(simulationOutput = sim_res, 
                           x = df$x, 
                           y = df$y)


# -----------------------------------------------------------------------------
# 6. permutation test, model loses significance, good. dismissing artifact
set.seed(42)
y_messUp <- sample(y_disc)

nb_mess <- glm.nb(y_messUp ~ df$n_immune + trt + x_scaled + offset(log(df$n_tumor)))


