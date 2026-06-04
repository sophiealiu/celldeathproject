# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/02/2026 4:59pm
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(dplyr)
library(glmnet)
library(pheatmap)
library(tidyr)

datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"
iso <- read.csv(file = file.path(datadir, "0528_22NMF_iso_40.csv"))

pd <- read.csv(file = file.path(datadir, "0602_22NMF_pd1_40incl.csv"))    


# -----------------------------------------------------------------------------
# preliminary checks of NMF factors
fact_cols <- names(pd)[12:33]                         # indices, NMF columns begin @12

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
# ASSUMPTIONS AND LIMITATIONS
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights
# - attempted to restore pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
x <- as.matrix(iso[, fact_cols])       # 2x control cond. & exp
x_scaled <- scale(x)                  # visible representation of RNA counts

y_binar <- pd$exist_dying             # binomial  [X]
y_disc <- iso$n_dying                  # overdispersion -> negative binomial
                                      # zero inflation overfit from performance stats
yd_test <- pd$n_lectin                # known vasculature for verification

y_cont <- pd$prop_dying               # beta, non-gaussian


# -----------------------------------------------------------------------------
# 2. basic comparative univariate odds-ratios
OR <- glm(
  reformulate(x_scaled, response = y_binar),
  family = "binomial"
)

OR_clean <- tidy(OR, conf.int = TRUE, exponentiate = TRUE)

ggplot(OR_clean %>% 
         filter(term != "(Intercept)"),          # only displaying the actual factors
       aes(x = reorder(term, estimate), y = estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) 
)


# -----------------------------------------------------------------------------
# 3. using counts considers sequencing depth
hist(y_disc, 
     breaks = 30,
     main = "Distribution of count dying",
     xlab = "number of dying cells in 40 micron vicinity")

# a. negative binomial AIC = 3994.1
library(MASS)
nb <- glm.nb(y_disc ~ x_scaled)

p_vals <- sort(
  coef(summary(nb_test))[, "Pr(>|z|)"], 
  decreasing = FALSE)
p_vals_adj <- p.adjust(p_vals, method = "BH")     # prefer benjamini but still optimistic

# b. elastic net is another way to see what factors are most predictive
library(glmnet)
# finding optimal penalization term using cross-validation to minimize MSE
crval <- cv.glmnet(x_scaled, y_disc, alpha = 0.5, nfolds = 10)
optim <- crval$lambda.min

elastic <- glmnet(x_scaled, y_disc, alpha = 0.5, lambda = optim)
coef(elastic)                       

# c. reduced model, only "significant" predictors. AIC = 4008.8
red_col <- fact_cols[c(1, 3, 11, 13, 22)]
x_red <- as.matrix(pd[, red_col])
x_red_sc <- scale(x_red)

nb_red <- glm.nb(y_disc ~ x_red_sc)

# d. refining the model nonlinearity. parsimony
x_red_rf <- x_red_sc
x_red_rf[,3] <- x_red_sc[,3] + I(x_red_sc[,3]^2)   # polynomial term

nb_rf <- glm.nb(y_disc ~ x_red_rf)

intthir <- x_red_sc[,3] * x_red_sc[,4]             # slope interaction
nb_int <- glm.nb(y_disc ~ x_red_sc + 
                          intthir)


# -----------------------------------------------------------------------------
# 4. Examining continuous models, however, proportion is oversmoothed
library(betareg)
pd_val <- pd %>%
  filter(pd$n_immune != 0)

xv <- as.matrix(pd_val[, fact_cols])
xv_scaled <- scale(xv)

yc2 <- pd_val$efficacy                             # depends on definition
                                                   # gamma if using dying counts per instead
beta <- betareg(yc2 ~ xv_scaled)  


# -----------------------------------------------------------------------------
# 5. diagnostic stats
# a. discrete
library(DHARMa)
sim_res <- simulateResiduals(nb_red)
plot(sim_res)       

# finding the culprit throwing off my data
for(i in 1:5) {                     # k defined upstream, number of predictors
  plotResiduals(sim_res, x_red_sc[, i],
                main = colnames(x_red_sc)[i])
}

# b. continuous
qqnorm(residuals(beta, type = "quantile"))
qqline(residuals(beta, type = "quantile"))


# -----------------------------------------------------------------------------
# 6. permutation test. Want to see the model fail
set.seed(42)
y_messUp <- sample(y_disc)

nb_mess <- glm.nb(y_messUp ~ x_scaled)


# -----------------------------------------------------------------------------
# 7. radii sensitivity
pd10 <- read.csv(file = file.path(datadir, "0602_22NMF_pd1_10incl.csv"))
pd20 <- read.csv(file = file.path(datadir, "0602_22NMF_pd1_20incl.csv"))  
pd80 <- read.csv(file = file.path(datadir, "0602_22NMF_pd1_80incl.csv"))  

x10 <- as.matrix(pd10[, fact_cols])
x10sc <- scale(x10)  

x20 <- as.matrix(pd20[, fact_cols])
x20sc <- scale(x20)  

x80 <- as.matrix(pd80[, fact_cols])
x80sc <- scale(x80)  

y10d <- pd10$n_dying
y20d <- pd20$n_dying
y80d <- pd80$n_dying

nb10 <- glm.nb(y10d ~ x10sc)
nb20 <- glm.nb(y20d ~ x20sc)
nb80 <- glm.nb(y80d ~ x80sc)


