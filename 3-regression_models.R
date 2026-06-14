# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/14/2026 10:01pm MDT
# Purpose: model iterations, regression dying tumor ~ factors
# -----------------------------------------------------------------------------

library(dplyr)
library(glmnet)
library(tidyr)

# importing files
iso <- read.csv("9NMF_iso40.csv")
pd1 <- read.csv("9NMF_pd140.csv")
df <- rbind(iso, pd1)


# -----------------------------------------------------------------------------
# preliminary checks of NMF factors
fact_cols <- names(df)[10:18]                         # indices, check head first

sorted_genes <- (sort(W[,9], decreasing = TRUE))      # summary of top genes
head(names(sorted_genes), 10)

sorted_factor <- list()
for (i in 1:length(fact_cols)) {
  ord <- order(W[, i], decreasing = TRUE)             # re-import gene factors
  
  sorted_factor[[i]] <- data.frame(
    gene   = rownames(W)[ord],
    weight = W[ord, i]
  )
  
  colnames(sorted_factor[[i]]) <- c(
    paste0("factor_", fact_cols[i], "_gene"),
    paste0("factor_", fact_cols[i], "_weight")
  )
}

df_weighted <- do.call(cbind, sorted_factor)
write.csv(head(df_weighted, 40), file.path(localdir, "factor_weights40.csv"),
          rownames = FALSE)

# ----------------------------------------------------------------------------_
# BEGINNING MODELING

# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights
# - pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
x <- as.matrix(df[, fact_cols])         
x_scaled <- scale(x)                    # visible representation of RNA counts. z-scores

y_binar <- df$exist_dying               
y_disc <- df$n_dying                    # overdispersion -> negative binomial
                                        # zero inflation overfit from performance stats
trt <- ifelse(df$sample == "apd1", 1, 0)
tum_off <- offset(log(df$n_tumor+1))    # offset term, don't treat same as coeffs


# -----------------------------------------------------------------------------
# 2. basic comparative univariate odds-ratios
OR <- glm(
  reformulate(df$n_immune + trt + x_scaled + tum_off, response = y_binar),
  family = "binomial"
)

OR_clean <- tidy(OR, conf.int = TRUE, exponentiate = TRUE)

ggplot(OR_clean %>% 
         aes(x = reorder(term, estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) 


# -----------------------------------------------------------------------------
# 3. using counts considers sequencing depth. iterate after diagnosis in 5
hist(y_disc, 
     breaks = 30,
     main = "Distribution of count dying",
     xlab = "number of dying cells in 40 micron vicinity")


# a. negative binomial. AIC = 4390.5
library(MASS)
nb <- glm.nb(y_disc ~ df$n_immune + trt + x_scaled + tum_off)

p_vals <- sort(                                   # first visualization
  coef(summary(nb))[, "Pr(>|z|)"], 
  decreasing = FALSE)
p_vals_adj <- p.adjust(p_vals, method = "BH")     # after benjamini still optimistic

# b. zero-inflation. AIC = 4459.6
library(glmmTMB)
zero <- glmmTMB(y_disc ~ n_immune + trt + x_scaled, 
                data = df, 
                family = nbinom2)

# c. geographical trends, Moran's I. Not combined w/ zero. AIC = 4358.5
library(mgcv)
spatial <- gam(y_disc ~ n_immune + trt + x_scaled + tum_off + 
                 s(cx, cy, bs = "gp"),     # spatial smoothing
                data = df,                         
                family = nb())           

# d. reduced model, only "significant" predictors. AIC = 4356.9
red_col <- fact_cols[c(1, 4)]
x_red <- as.matrix(df[, red_col])
x_red_sc <- scale(x_red)

spatial_red <- gam(y_disc ~ n_immune + trt + x_red_sc + tum_off +
                     s(cx, cy, bs = "gp"),
                   data = df,
                   family = nb())

# e. elastic net re-visualization. consistent effect size
library(glmnet)
# finding optimal penalization term
crval <- cv.glmnet(df$n_immune + trt + x_scaled + tum_off, y_disc, 
                   alpha = 0.5, nfolds = 10)
optim <- crval$lambda.min

elastic <- glmnet(df$n_immune + trt + x_scaled + tum_off, y_disc, 
                  alpha = 0.5, lambda = optim)
coef(elastic) 


# -----------------------------------------------------------------------------
# 5. diagnostic stats
library(DHARMa)
sim_resR <- simulateResiduals(spatial_red)
plot(sim_res)

# ratio of zeroes, okay need to do this.
testZeroInflation(sim_res)

# finding the culprit throwing off my data
full_vars <- cbind(x_red_sc, df$n_immune, trt)
for(i in 1:4) {                                   # number of reduced factors
  plotResiduals(sim_res, full_vars[, i],
                main = colnames(full_vars)[i])
}

testSpatialAutocorrelation(simulationOutput = sim_res0, 
                           x = df$cx, 
                           y = df$cy)

# c. visualization of treatment interaction
ggplot(df, aes(x = scale(X1), y = y_disc,
               color = factor(trt,
                              levels = c(0,1),
                              labels = c("control","treated")))) +
  geom_smooth(method = "loess") +
  scale_color_manual(
    values = c("control" = "steelblue",
               "treated" = "pink"),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, 1.645)) + 
  labs(
    x = "factor 1 expression (standardized)",
    y = "counts dying in 40 micron radius",
    color = NULL
  )


# -----------------------------------------------------------------------------
# 6. permutation test. Want to see the model lose
set.seed(42)
y_messUp <- sample(y_disc)

nb_mess <- glm.nb(y_messUp ~ df$n_immune + trt + x_scaled + tum_off)


# -----------------------------------------------------------------------------
# 7. radii sensitivity. base 40
df10 <- read.csv("NMF_rad10.csv")
df20 <- read.csv("NMF_rad20.csv") 
df80 <- read.csv("NMF_rad80.csv")

x10 <- as.matrix(df10[, fact_cols])
x10sc <- scale(x10)  

x20 <- as.matrix(df20[, fact_cols])
x20sc <- scale(x20)  

x80 <- as.matrix(df80[, fact_cols])
x80sc <- scale(x80)  

y10d <- df10$n_dying
y20d <- df20$n_dying
y80d <- df80$n_dying

nb10 <- glm.nb(y10d ~ x10sc)
nb20 <- glm.nb(y20d ~ x20sc)
nb80 <- glm.nb(y80d ~ x80sc)

yd_test <- df$n_lectin            # vasculature for verification. cell type proximity
nb_test <- glm.nb(yd_test ~ x_scaled)

