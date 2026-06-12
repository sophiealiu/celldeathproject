# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/11/2026 5:58pm
# Purpose: Creating NMF factors from Visium data and optimizing for downstream
# -----------------------------------------------------------------------------

library(Matrix)
library(ggplot2)
library(RcppML)                # RcppML (and Gene NMF) for sparse matrices
library(Seurat)

# importing files
datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"

iso_raw <- readRDS(file = file.path(datadir, "Seurat objects", "iso_raw.rds"))
pd1_raw <- readRDS(file = file.path(datadir, "Seurat objects", "pd1_raw.rds"))


# -----------------------------------------------------------------------------
# 1. expression matrix X before decomposing. merging here because raw unable
mat_iso <- GetAssayData(
  iso_raw,
  assay = "Spatial.008um",
  layer = "counts"
)

mat_pd1 <- GetAssayData(
  pd1_raw,
  assay = "Spatial.008um",
  layer = "counts"
)

# need to have different. double check which ones it is
colnames(mat_iso)  <- paste0("iso_", colnames(mat_iso))
colnames(mat_pd1) <- paste0("pd1_", colnames(mat_pd1))

merged_mat <- cbind(mat_iso, mat_pd1)


# -----------------------------------------------------------------------------
# 2. testing k-ranks
# a. avg loss across runs gives us optimal rank at curvature
library(inflection)
ks <- 2:20
n <- 10                                              # 10 for estimation

loss_matrix <- matrix(nrow = length(ks), ncol = n)

for (i in seq_along(ks)) {
  for (r in 1:n) {                                   # 7 is my favorite number
    model <- nmf(merged_mat, k = i, seed = 7*n, verbose = FALSE)      
    loss_matrix[i, r] <- mse(merged_mat, model$w, 
                                         model$d, 
                                         model$h)
  }
}

mean_loss <- rowMeans(loss_matrix)
k_opt <- uik(x = ks, y = mean_loss)

# b. stability score at optimal rank
n3 <- 30                                             # here, CLT
W_list <- vector("list", n3)
for (i in 1:n3) {
  mod_stab <- nmf(merged_mat, k = k_opt, seed = 7*i, verbose = FALSE)
  W_list[[i]] <- mod_stab$w
}

pair_stab <- c()

for (i in 1:(n3 - 1)) {
  for (j in (i + 1):n3) {
    # pearson correlation between runs i& j. want close to 1
    cor_mat <- cor(W_list[[i]], W_list[[j]])
    
    matched_corrs <- numeric(ks)
    temp_mat <- cor_mat
    
    # greedy matching computationally more feasible. 
    for (f in 1:ks) {
      max <- which(temp_mat == max(temp_mat), arr.ind = TRUE)[1, ]
      matched_corrs[f] <- temp_mat[max[1], max[2]]
      temp_mat[max[1], ] <- -1
      temp_mat[, max[2]] <- -1
    }
    
    pair_stab <- c(pair_stab, mean(matched_corrs))
  }
}

final_stab <- mean(pair_stab)


# -----------------------------------------------------------------------------
# 4. final model, optimal rank was 10
best_mod <- nmf(merged_mat, k = 10, seed = 42)

# re-aligning row and column names (lost before)
rownames(best_mod$w) <- rownames(merged_mat)
colnames(best_mod$h) <- colnames(merged_mat)

# renaming for ease of reference
W <- best_mod$w
H <- best_mod$h  


# -----------------------------------------------------------------------------
# 5. factors are orthogonality. using cosine similarity btwn factors
library(proxyC)
factor_sim <- simil(
  t(W),                                            # required to transpose to get factor x gene
  method = "cosine"
)

library(pheatmap)
pheatmap(factor_sim, 
         clustering_distance_roW_list = "correlation",
         clustering_distance_cols = "correlation")

         
# -----------------------------------------------------------------------------
# 6. finalizing, binding spatial data
pd1_vis_coords <- readRDS(file.path(datadir, "pd1_vis_coords.rds"))
iso_vis_coords <- readRDS(datadir, "iso_vis_coords.rds")

rownames(pd1_vis_coords)  <- paste0("pd1_", rownames(pd1_vis_coords))
rownames(iso_vis_coords) <- paste0("iso_", rownames(iso_vis_coords))

tH <- t(H)   # required to transpose to line up. colnamesH holds both

common_pd1 <- rownames(pd1_vis_coords)[
  rownames(pd1_vis_coords) %in% colnames(H)
]

pd1_joined <- cbind(
  pd1_vis_coords[common_pd1, , drop = FALSE],
  tH[common_pd1, , drop = FALSE]
)

# scale from pixels to microns for euclidean spatial alignment
pd1_joined$x <- pd1_joined$x / 1.5454              # scaling factor from imaging metadata
pd1_joined$y <- pd1_joined$y / 1.5454              

# export, repeating above for isotype control
write.csv(pd1_joined, file.path(datadir, "0611_NMFs.csv"))

         
