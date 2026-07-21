# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(DHARMa)       # some packages will have to be forced in using remotes::install_github
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(MASS)
library(mgcv)
library(tidyr)

# datadir <- "path/to/your/working/directory"
datadir <- "I:/Hu Lab/Sophie/1. Cell death/all final data"

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
# distribution explains family choice. Supplemental figure
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

# b. spatial correction for NB model to address autocorrelation, cautious of overfitting. 
# excess overfit possible with smoothing taking credit.
spatial <- gam(
  n_dying ~ n_immune + trt + fact_sc +
    s(x, y, bs = "gp", k = 30) +            # smoothing using coordinates depends on power of gaussian process
    offset(log(n_tumor)),
  data = df_merged,
  family = nb()
)
summary(spatial)

# parametric p-vals adjusted
p_para <- summary(spatial)$p.pv                   # extracts p-vals from spatial model
p_adj <- p.adjust(p_para, method = "BH")          # benjamini-hochberg correction
print(p_para)
print(p_adj)

# c. incorporating treatment interactions into NB model
nb_int <- glm.nb(
  n_dying ~ n_immune + trt * fact_sc + offset(log(n_tumor)),
  data = df_merged
)
summary(nb_int)


# -----------------------------------------------------------------------------
# 3. diagnostic stats to check model. Supplemental figure
sim_res <- simulateResiduals(spatial)
plot(sim_res)              

testZeroInflation(sim_res)
testSpatialAutocorrelation(simulationOutput = sim_res,         # Moran's I
                           x = df_merged$x, 
                           y = df_merged$y)

nice_labels <- c(
  "fact_scfactor1"    = "1",
  "fact_scfactor2"    = "2",
  "fact_scfactorX3"    = "3",
  "fact_scfactor4"    = "4",
  "fact_scfactor5"    = "5",
  "fact_scfactor6"    = "6",
  "fact_scfactor7"    = "7",
  "fact_scfactor8"    = "8",
  "fact_scfactor9"    = "9"
)

nb_nice <- as.data.frame(summary(nb_test)$coeff)
nb_nice <- nb_nice %>%
  rename(
    `effect size` = "Estimate",
    p_value = `Pr(>|z|)`
  )

nb_nice$"effect size" <- as.numeric(nb_nice$"effect size")
nb_nice$BH <- p.adjust(nb_nice$p_value, method = "BH")

nb_nice <- nb_nice[grepl("^xsc", rownames(nb_nice)), ]  
rownames(nb_nice) <- nice_labels[rownames(nb_nice)]

# Figure 5, panel E. visualization of BH
EnhancedVolcano(nb_nice,
                lab = rownames(nb_nice), 
                selectLab = rownames(nb_nice),
                x = 'effect size',
                y = 'BH',  pCutoff = 0.05, FCcutoff = 0.5,
                xlim = c(-0.6,0.6),
                ylim = c(0,30),
                xlab = "effect size (regression coefficient)",
                title = "Benjamini-Hochberg adjusted p-value vs effect size",
                
                labSize = 7,
                pointSize = 4,
                
                maxoverlapsConnectors = Inf, 
                drawConnectors = TRUE
)


# -----------------------------------------------------------------------------
# 4. permutation to test random noise effects
# would like to see permutation model lose significance. checking diagnostics as well
set.seed(42)
df_perm <- df_merged
df_perm$n_dying <- sample(df_merged_perm$n_dying)

nb_mess <- glm.nb(
  n_dying ~ n_immune + trt + fact_sc + offset(log(n_tumor)),
  data = df_merged_perm
)
summary(nb_mess)
