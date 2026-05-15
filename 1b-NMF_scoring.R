# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/15 10:39am
# Purpose: Creating NMF factors from Visium data and determining the best scaling method
# -----------------------------------------------------------------------------
# X ≈ WH (W is defining genes, H is spatial location)

library(Matrix)
library(ggplot2)
library(glmnet)
library(patchwork)
library(RcppML)
library(Seurat)

        
# importing files
datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"

pd1_raw <- readRDS(file = file.path(datadir, "Seurat objects", "pd1_raw.rds"))
pd1_n16 <- readRDS(file = file.path(datadir, "Seurat objects", "pd1_n16.rds"))
                          # larger-bin, normalized for testing w/o computational lag


# -----------------------------------------------------------------------------
# 1. Normalization
emat <- GetAssayData(pd1_raw, 
                     assay = "Spatial.008um",
                     layer = "counts")

# per-gene normalization
emat_sub <- emat[rowSums(emat) >0, ]          # some genes not present at all
gene_max <- apply(emat_sub, 1, max)

emat_by_gene <- sweep(emat, 1, gene_max, "/") # sweep 1 for rows

# per-bin normalization
bin_sums <- colSums(emat)
emat_by_bin <- sweep(emat, 2, bin_sums, "/")  # sweep 2 for columns


# -----------------------------------------------------------------------------
# 2. Solving expression matrix
# emat16 <- GetAssayData(pd1_n16, 
#                        assay = "Spatial.016um",
#                        layer = "data")  


# factorization stops when X ≈ WH is good enough 
ks <- c(5, 10, 20, 50)                   # varying number of programs for optimal model.
NMF_fits <- list()

# iterate/re-seed for validation
set.seed(42)                             # the answer to to the ultimate question of life, 
                                         # the universe, and everything.
for (k in ks) {
    NMF_fits[[paste0("k", k)]] <- nmf(emat_scaled,
                                      k = k,
                                      verbose = FALSE
  )
  rownames(NMF_fits[[paste0("k", k)]]$w) <- rownames(emat_scaled)
  colnames(NMF_fits[[paste0("k", k)]]$h) <- colnames(emat_scaled)
}                                                    

fit20 <- NMF_fits[["k20"]]
head(sort(fit20$w[, 12], decreasing = TRUE), 20)


# -----------------------------------------------------------------------------
# 3. preliminary visualization
fit <- NMF_fits[["k20"]]

plots <- list()
for (i in c(1,2,7,9,12)) {
  pd1_raw[[paste0("Factor", i)]] <- fit$h[i, ]

  p <- SpatialFeaturePlot(
    pd1_raw,
    features = paste0("Factor", i))
    #alpha = 0.5)
    
  plots[[length(plots) + 1]] <- p
  }
  
wrap_plots(plots, ncol = 3)


