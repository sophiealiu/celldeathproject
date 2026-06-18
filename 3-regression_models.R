# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/15/2026 4:02pm MDT
# Purpose: model iterations, regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)

# data from block 2 using 1.1
# iso <- read.csv("sig_iso40.csv"))
# pd1 <- read.csv("sig_pd140.csv"))
# df <- rbind(iso, pd1)

# iso <- read.csv("all_genes_iso40.csv"))
# pd1 <- read.csv("all_genes_pd140.csv"))
# df <- rbind(iso, pd1)

# iso <- read.csv("DE_iso40.csv"))
# pd1 <- read.csv("DE_pd140.csv"))
# df <- rbind(iso, pd1)

# 1.2 MAIN
iso <- read.csv("9NMF_iso40.csv"))
pd1 <- read.csv("9NMF_pd140.csv"))
df <- rbind(iso, pd1)

# 1.3
# iso <- read.csv("9NSF_iso40.csv"))
# pd1 <- read.csv("9NSF_pd140.csv"))
# df <- rbind(iso, pd1)


# -----------------------------------------------------------------------------
# preliminary checks of NMF factors
fact_cols <- names(df)[10:18]                         # indices, check head first

# summary of top genes
W <- read.csv(file.path(localdir, "0613_9NMF_W.csv"),row.names = 1)

top_idx <- order(W$V8, decreasing = TRUE)
rownames(W)[top_idx][1:10]


# ----------------------------------------------------------------------------_
# BEGINNING MODELING

# What are our assumptions?
# - know data is nonlinear and non-normal
# - means smoothing need to incorporate density covariate
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights
# - pseudo-independence by sampling non-overlapping disks


# -----------------------------------------------------------------------------
# 1. setting vars and family
x <- as.matrix(df[, fact_cols])         # replace with cell cols or gene cols based on which block 1
x_scaled <- scale(x)                    # visible representation of RNA counts. z-scores

# distribution explains family choice
hist(y_disc, 
     breaks = 30,
     main = "Distribution of count dying",
     xlab = "number of dying cells in 40 micron vicinity")

y_binar <- df$exist_dying               # binomial  [X]
y_disc <- df$n_dying                    # overdispersion -> negative binomial
                                        # zero inflation overfit from performance stats
trt <- ifelse(df$sample == "pd1-9", 1, 0)
tum_off <- offset(log(df$n_tumor+1))    # offset term, don't treat same as coeffs


# -----------------------------------------------------------------------------
# 2. weak signals from individual genes and signatures (DE) block 1.1
# simple monotonic relationship strength
spearman <- cor(x_scaled, y_disc, method = "spearman")

names(spearman) <- cell_cols
spearman_ordered <- sort(spearman, decreasing = TRUE)

# high collinearity
cor_mat <- cor(df[cell_cols], use = "pairwise.complete.obs")
cor_table <- round(cor_mat, 2)
pheatmap(cor_mat)


# -----------------------------------------------------------------------------
# 3. basic comparative univariate odds-ratios after moving to block 1.2
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
# 4. using counts considers sequencing depth. iterate after diagnosis in 5
# a. negative binomial. AIC = 4410.3
library(MASS)
nb <- glm.nb(y_disc ~ df$n_immune + trt + x_scaled + tum_off)

p_vals <- sort(                                   # first visualization
  coef(summary(nb))[, "Pr(>|z|)"], 
  decreasing = FALSE)
p_vals_adj <- p.adjust(p_vals, method = "BH")     # after benjamini still optimistic

# b. zero-inflation. AIC = 4428.1
library(glmmTMB)
zero <- glmmTMB(y_disc ~ n_immune + trt + x_scaled, 
                data = df, 
                family = nbinom2)

# c. geographical trends, Moran's I. Not combined w/ zero. AIC = 4342.3
library(mgcv)
spatial <- gam(y_disc ~ n_immune + trt + x_scaled + tum_off + 
                 s(x, y, bs = "gp"),     # spatial smoothing
                data = df,                         
                family = nb())           

# d. reduced model, only "significant" predictors. AIC = 4354.5
red_col <- fact_cols[c(1, 4)]
x_red <- as.matrix(df[, red_col])
x_red_sc <- scale(x_red)

spatial_red <- gam(y_disc ~ n_immune + trt + x_red_sc + tum_off +
                     s(x, y, bs = "gp"),
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

# f. interaction (careful overfitting). AIC = 4314.8
inter1 <- scale(df$X1)*trt
inter4 <- scale(df$X4)*trt
sr_int <- gam(y_disc ~ n_immune + trt + x_red_sc + tum_off + inter1 + inter4
                     + s(x, y, bs = "gp"),
                   data = df,
                   family = nb())

# -----------------------------------------------------------------------------
# 5. diagnostic stats
library(DHARMa)
sim_resR <- simulateResiduals(sr_int)
plot(sim_resR)

# a. ratio of zeroes and autocorrelation
sim_res0 <- simulateResiduals(nb)
testZeroInflation(sim_res0)
testSpatialAutocorrelation(simulationOutput = sim_res0, 
                           x = df$x, 
                           y = df$y)

# finding the culprit throwing off my data
full_vars <- cbind(x_red_sc, trt)
for(i in 1:4) {                                   # number of reduced factors
  plotResiduals(sim_resR, full_vars[, i],
                main = colnames(full_vars)[i])
}

# b. visualization of treatment interaction. not that informative but potential
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
  coord_cartesian(xlim = c(0, 2)) + 
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
# 7. vasculature for verification. cell type proximity
iso10 <- read.csv("9NMF_iso10.csv")
pd110 <- read.csv("9NMF_pd110.csv")
df10 <- rbind(iso10, pd110)

x10 <- as.matrix(df10[, fact_cols])
x10sc <- scale(x10)
trt10 <- ifelse(df10$sample == "pd1-9", 1, 0)
tum10_off <- offset(log(df10$n_tumor+1)) 

yd_test <- df10$n_lectin            
nb_test <- glm.nb(yd_test ~ x10sc)

y10d <- df10$n_dying
nb10 <- glm.nb(y10d ~ df10$n_immune + trt + tum10_off+ x10sc)


# -----------------------------------------------------------------------------
# 8. more radii sensitivity
iso20 <- read.csv("9NMF_iso20.csv")
pd120 <- read.csv("9NMF_pd120.csv")
df20 <- rbind(iso20, pd120)

iso80 <- read.csv("9NMF_iso80.csv")
pd180 <- read.csv("9NMF_pd180.csv")
df80 <- rbind(iso20, pd120)  

x20 <- as.matrix(df20[, fact_cols])
x20sc <- scale(x20)  
trt20 <- ifelse(df20$sample == "pd1-9", 1, 0)
tum20_off <- offset(log(df20$n_tumor+1)) 

x80 <- as.matrix(df80[, fact_cols])
x80sc <- scale(x80)  
trt80 <- ifelse(df80$sample == "pd1-9", 1, 0)
tum80_off <- offset(log(df80$n_tumor+1)) 

y20d <- df20$n_dying
y80d <- df80$n_dying

nb20 <- glm.nb(y20d ~ df20$n_immune + trt20 + tum20_off+ x20sc)
nb80 <- glm.nb(y80d ~ df80$n_immune + trt80 + tum80_off+ x80sc)

