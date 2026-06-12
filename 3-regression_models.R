# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/11/2026 11:30am
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(dplyr)
library(glmnet)
library(pheatmap)
library(tidyr)

datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"
df <- read.csv(file = file.path(datadir, "0528_22NMF_iso_40.csv"))


# -----------------------------------------------------------------------------
# preliminary checks of NMF factors
fact_cols <- names(df)[13:32]                         # indices, check head first

sorted_genes <- (sort(W[,17], decreasing = TRUE))     # summary of top genes
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
write.csv(head(df_weighted, 40), file.path(datadir, "factor_weights40.csv"),
          rownames = FALSE)

# ----------------------------------------------------------------------------_
# BEGINNING MODELING

# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights
# - attempted to restore pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
x <- as.matrix(df[, fact_cols])         # 2x repeated control condition & treated
x_scaled <- scale(x)                    # visible representation of RNA counts. z-scores

y_binar <- df$exist_dying               # binomial  [X]
y_disc <- df$n_dying                    # overdispersion -> negative binomial
                                        # zero inflation overfit from performance stats

y_cont <- df$prop_dying                 # non-gaussian
dens_correct <- offset(log(df$tumor))   # incorporating tumor density effects, 
                                        # more tumor likely more dying


# -----------------------------------------------------------------------------
# 2. basic comparative univariate odds-ratios
OR <- glm(
  reformulate(n_immune + trt + x_scaled + dens_corect, response = y_binar),
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


# a. negative binomial. AIC = 3994.1
library(MASS)
nb <- glm.nb(y_disc ~ n_immune + trt + x_scaled + dens_corect)

p_vals <- sort(                                   # first visualization
  coef(summary(nb))[, "Pr(>|z|)"], 
  decreasing = FALSE)
p_vals_adj <- p.adjust(p_vals, method = "BH")     # after benjamini still optimistic

# b. to control for geographical trends from Moran's I. AIC = 3915.2
library(mgcv)
spatial <- gam(y_disc ~ n_immune + trt + x_scaled + dens_corect + 
                        s(cx, cy, bs = "gp"),     # spatial smoothing
               data = df,                         
               family = nb())

# c. reduced model, only "significant" predictors. AIC = 3959.2, not worrisome
red_col <- fact_cols[c(1, 3, 7, 11, 17, 22)]
x_red <- as.matrix(df[, red_col])
x_red_sc <- scale(x_red)

spatial_red <- gam(y_disc ~ n_immune + trt + x_red_sc + dens_corect +
                            s(cx, cy, bs = "gp"),
                   data = df,
                   family = nb())

# d. elastic net re-visualization of collinearity and significant factors. 
library(glmnet)
# finding optimal penalization term that minimizes MSE
crval <- cv.glmnet(n_immune + trt + x_scaled + dens_corect, y_disc, 
                   alpha = 0.5, nfolds = 10)
optim <- crval$lambda.min

elastic <- glmnet(n_immune + trt + x_scaled + dens_corect, y_disc, 
                  alpha = 0.5, lambda = optim)
coef(elastic) 

# e. adding interaction terms. AIC = 3955.4
inter24 <- x_red_sc[,2] * x_red_sc[,4]             # slope interaction
spatial_int <- gam(y_disc ~ n_immune + trt + x_red_sc + inter24 + dens_corect + 
                            s(cx, cy, bs = "gp"),
                   data = df,
                   family = nb())

# f. addressing the curved residuals. parsimony!! AIC = 3945.9
x_red_rf <- x_red_sc
x_red_rf[,2] <- x_red_sc[,2] + I(x_red_sc[,2]^2)   # polynomial term
spatial_intC <- gam(y_disc ~ n_immune + trt + x_red_rf + inter24 + dens_corect
                             s(cx, cy, bs = "gp"),
                   data = df,
                   family = nb())


# -----------------------------------------------------------------------------
# 4. Examining continuous models, however, proportion is oversmoothed
library(betareg)
beta <- betareg(y_cont ~ n_immune + trt + x_scaled + dens_corect) 


# -----------------------------------------------------------------------------
# 5. diagnostic stats
# a. discrete
library(DHARMa)
sim_res <- simulateResiduals(spatial_red)
plot(sim_res)       

# finding the culprit throwing off my data
for(i in 1:6) {                                   # number of reduced factors
  plotResiduals(sim_res, x_red_rf[, i],
                main = colnames(x_red_rf)[i])
}

# b. continuous
qqnorm(residuals(beta, type = "quantile"))
qqline(residuals(beta, type = "quantile"))

# c. checking spatial autocorrelation
library(sdfep)
coords <- as.matrix(df[, c("cx", "cy")])
nb <- knn2nb(knearneigh(coords, k = 4))           # baseline querying 4 nearest neighbors
weights <- nb2listw(nb, style = "W")
testSpatialAutocorrelation(sim_res, x = df$cx, y = df$cy)


# -----------------------------------------------------------------------------
# 6. permutation test. Want to see the model fail
set.seed(42)
y_messUp <- sample(y_disc)

nb_mess <- glm.nb(y_messUp ~ x_scaled)


# -----------------------------------------------------------------------------
# 7. radii sensitivity
df10 <- read.csv(file = file.path(datadir, "0602_22NMF_df1_10incl.csv"))
df20 <- read.csv(file = file.path(datadir, "0602_22NMF_df1_20incl.csv"))  
df80 <- read.csv(file = file.path(datadir, "0602_22NMF_df1_80incl.csv"))  

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

yd_test <- df10$n_lectin            # vasculature for verification. cell type proximity
nb_test <- glm.nb(yd_test ~ x10sc)

