# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/15/2026
# Purpose: separate file, attempt to fit regression dying tumor ~ NMF factors
# -----------------------------------------------------------------------------

library(dplyr)
library(glmnet)
library(MASS)
library(pheatmap)
datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"
overall <- read.csv(file = file.path(datadir, "0514_NMF20calc.csv"))     


# What are our assumptions?
# - know data is nonlinear and non-normal
# - noisy experimental
# - better but still non-orthogonal when using NMF. raw counts, per-gene weights
# - many samples but independence violated because of spatial overlap

# -----------------------------------------------------------------------------
# 1. Model iteration
fact_cols <- names(overall)[12:31]                         # indices, NMF columns begin @12
x <- as.matrix(overall[, fact_cols])


# y <- overall$exist_dying              # binomial regression [N] from inspection
# y <- overall$prop_dying               # linear gaussian     [N]
y <- overall$n_dying                    # poisson             [N] OR
                                        # negative binomial   [possible]
# - may require transformation of response variable addressing variance
# - poisson dispersion ~ 2. residuals fan out


# a. general NB model
nb <- glm.nb(y ~ x)
summary(nb)

# b. okay this makes it worse
x_red <- as.matrix(overall[, fact_cols[c(1, 2, 7, 9, 12)]]) # reduced model using first sig

# c. adding mild nonlinearity (splines join multiple polynomial segments)
nb <- glm.nb(y ~ poly(x,2))


# plotting
plot(nb$fitted.values, residuals(nb, type = "pearson"))
abline(h = 0, col = "red")

library(DHARMa)
qq <- simulateResiduals(fittedModel = nb)
plot(qq)


# -----------------------------------------------------------------------------
# 2a. Distribution of signal across cell type
cell_types <- unique(overall$cell_type)
factor_sums <- colSums(overall[, fact_cols])

percent_sig <- list()
for (type in cell_types) {
  subset_df <- overall[overall$cell_type == type, ]
  
  type_sums <- colSums(subset_df[, fact_cols])
  percent <- (type_sums / factor_sums) * 100
  
  percent_sig[[type]] <- percent
}

percent_sig <- do.call(rbind, percent_sig)


# 2b. Distribution of cell across factors
percent_type <- list()
for (type in cell_types) {
  subset_df <- overall[overall$cell_type == type, ]
  
  by_factor <- colSums(subset_df[, fact_cols])
  total <- sum(by_factor)
  
  percent <- (by_factor/total) * 100
  
  percent_type[[type]] <- percent
}
percent_type <- do.call(rbind, percent_type)
percent_type


# -----------------------------------------------------------------------------
# 3. Visualization
mat <- as.matrix(percent_type)

pheatmap(mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE)

