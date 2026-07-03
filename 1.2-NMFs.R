# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 06/11/2026 5:58pm MDT
# Purpose: Creating NMF factors from Visium data and optimizing for downstream
# -----------------------------------------------------------------------------

library(clue)
library(ggplot2)
library(inflection)
library(progress)
library(Matrix)
library(RcppML)        # our data is too sparse to use regular NMF package
library(Seurat)

# importing files. no assumptions made, hence raw data. 10x recommended 8 micron bins
datadir <- "path/to/your/working/directory"
iso_raw <- readRDS(file.path(datadir,"iso_raw.rds"))
pd1_raw <- readRDS(file.path(datadir,"pd1_raw.rds"))

# -----------------------------------------------------------------------------
# 1. expression matrix before decomposing. 
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

# prefixing required so we can separate downstream
colnames(mat_iso)  <- paste0("iso_", colnames(mat_iso))
colnames(mat_pd1) <- paste0("pd1_", colnames(mat_pd1))

merged_mat <- cbind(mat_iso, mat_pd1)


# -----------------------------------------------------------------------------
# 2. testing k-ranks. 
# a. avg MSE dim. unit invariant knee, optimal rank = 9
ks <- 2:20                                               # from manual sweeping by 5s, k < 20
nruns_opt <- 10
loss_matrix <- matrix(nrow = length(ks), ncol = nruns_opt)

pb <- progress_bar$new(                                  # it takes ~ 3hrs, sorry
  format = "[:bar] :percent, approx. time left: :eta",
  total = length(ks)*nruns_opt,                          
  clear = FALSE, 
  width = 60
)

for (i in seq_along(ks)) {
  for (r in 1:nruns_opt) {
    model <- nmf(merged_mat, k = ks[i], seed = r,     
                 verbose = FALSE)
    loss_matrix[i, r] <- mse(merged_mat, model$w,        
                                         model$d,
                                         model$h)
    pb$tick()
  }
}

mean_loss <- rowMeans(loss_matrix)
k_opt <- uik(x = ks, y = mean_loss)

# b. stability score using Hungarian algorithm, gene weights across stochastic regeneration. 
# result: 0.9498465 yay close to 1
nruns_stab <- 30                                   # here CLT
pair_stab <- c()
W_list <- vector("list", nruns_stab)               # features

for (r in 1:nruns_stab) {
  save <- nmf(merged_mat, k = k_opt, seed = r)
  W_list[[r]] <- save$w
}

for (i in 1:(nruns_stab - 1)) {
  for (j in (i + 1):nruns_stab) {
    cor_mat <- abs(cor(W_list[[i]], W_list[[j]]))        # magnitude pref.
    assign <- solve_LSAP(cor_mat, maximum = TRUE)
    
    matched <- cor_mat[cbind(1:k_opt, assign)]
    pair_stab <- c(pair_stab, mean(matched))
  }
}


# -----------------------------------------------------------------------------
# 3. re-aligning row and column names (lost before). 
# just re-using the last save, but re-running doesn't take too long.
rownames(save$w) <- rownames(merged_mat)
colnames(save$h) <- colnames(merged_mat)

# renaming for ease of reference, file creation for GSEA block 5
W <- save$w
H <- save$h  
write.csv(W, "NMF_W.csv")


# -----------------------------------------------------------------------------
# 4. finalizing/ binding spatial data. 
pd1_vis_coords <- GetTissueCoordinates(pd1_raw, image = NULL)
iso_vis_coords <- GetTissueCoordinates(iso_raw, image = NULL)

# proper prefixing
rownames(pd1_vis_coords)  <- paste0("pd1_", rownames(pd1_vis_coords))
rownames(iso_vis_coords) <- paste0("iso_", rownames(iso_vis_coords))

tH <- t(H)                                  # transpose to line up weights

common_pd1 <- rownames(pd1_vis_coords)[
  rownames(pd1_vis_coords) %in% colnames(H)
]
common_iso <- rownames(iso_vis_coords)[
  rownames(iso_vis_coords) %in% colnames(H)
]

pd1_joined <- cbind(
  pd1_vis_coords[common_pd1, , drop = FALSE],
  tH[common_pd1, , drop = FALSE]
)

iso_joined <- cbind(
  iso_vis_coords[common_iso, , drop = FALSE],
  tH[common_iso, , drop = FALSE]
)

# export both individually, bind at regression step. keep in pixels
write.csv(pd1_joined, file.path(datadir, "NMF_Apd1.csv"))
write.csv(iso_joined, file.path(datadir, "NMF_isoC.csv"))

