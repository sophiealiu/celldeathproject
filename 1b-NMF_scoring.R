# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/21/2026 3:58pm
# Purpose: Creating NMF factors from Visium data and optimizing for downstream
# -----------------------------------------------------------------------------
# X ≈ WH (W is feature matrix/# gene components, H is coefficient matrix/weights)

library(Matrix)
library(ggplot2)
library(RcppML)        # faster computationally but treats it as sparse matrix
library(Seurat)

        
# importing files
datadir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"

iso_raw <- readRDS(file = file.path(datadir, "Seurat objects", "iso_raw.rds"))
pd1_raw <- readRDS(file = file.path(datadir, "Seurat objects", "pd1_raw.rds"))


# -----------------------------------------------------------------------------
# 1. Expression matrix X before decomposing. must merge here instead of upstream
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

# removing columns I don't want
mat_iso <- mat_iso[, Matrix::colSums(mat_iso) > 0]
mat_pd1 <- mat_pd1[, Matrix::colSums(mat_pd1) > 0]
mat_iso <- na.omit(mat_iso)
mat_pd1 <- na.omit(mat_pd1)

# merging genes
all_genes <- union(rownames(mat_iso), rownames(mat_pd1))

mat_iso <- mat_iso[match(all_genes, rownames(mat_iso)), ]
mat_pd1 <- mat_pd1[match(all_genes, rownames(mat_pd1)), ]

merged_mat <- cbind(mat_iso, mat_pd1)


# -----------------------------------------------------------------------------
# 2. testing different k-ranks to optimize least error. when MSE improvement flattens 
# ideally this also is stable with 30 runs
library(inflection)

ks <- 2:50
fit <- list()
merged_mse <- list()

for (i in seq_along(ks)) {
  k <- ks[i]
  fit[[i]] <- nmf(merged_mat, k = k, verbose = FALSE)
  
  merged_mse[[i]] <- mse(merged_mat, 
                         fit[[i]]$w, fit[[i]]$d, fit[[i]]$h, 
                         mask_zeros = FALSE)
}

names(merged_mse) <- ks

plot(2:50, delta, type = "b",
     xlab = "k",
     ylab = "change in MSE")

merged_mseNum <- as.numeric(merged_mse)
delta <- diff(merged_mseNum)


# -----------------------------------------------------------------------------
# 3. manually checking stability. re-seeding. from optimization, k = 16

n <- 30                                # CLT, number of runs
iterate_fit16 <- vector("list", n)

for (i in 1: n) {
  set.seed(7*i)                        # my lucky number is 7. 
  iterate_fit16[[i]] <- nmf(merged_mat, k = 17, verbose = FALSE)
  
  rownames(iterate_fit16[[i]]$w) <- rownames(merged_mat)
  colnames(iterate_fit16[[i]]$h) <- colnames(merged_mat)
}

# pair-wise comparison against reference
W_list <- lapply(iterate_fit16, function(x) x$w)
W_ref <- W_list[[1]]        

cor_matrix <- list()
for (i in 1:n) {
  W <- W_list[[i]]
  cor_matrix[[i]] <- cor(W, W_ref)
}

# getting stability score (how it compares)
stability <- sapply(cor_matrix, function(avg) {
  mean(apply(avg, 1, max))            # greedy algorithm matching
})

# bar graph for similarity. 0.75 as threshold 
run_num <- 1:length(stability)
bar_colors <- c("#9678B6", rep("lightgray", 29))

barplot(stability,
        names.arg = run_num,
        col = bar_colors,
        ylim = c(0, 1),
        xlab = "NMF run",
        ylab = "proportion similar",
        main = "NMF similarity across 30 runs (k = 16), reference column 1")
par(cex.main = 0.75)
par(cex.axis = 0.5)
abline(h = 0.75, col = "blue", lwd = 2, lty =2)


# -----------------------------------------------------------------------------
# 4. viewing if the factors are orthogonal. using cosine similarity btwn vectors
library(proxyC)
factor_sim <- simil(
  t(W),                        # required to transpose to get factor x gene
  method = "cosine"
)


library(corrplot)
corrplot(
  as.matrix(factor_sim),
  method = "color",
  type = "upper"
)


# -----------------------------------------------------------------------------
# 5. NMF on optimal k rank (16) on combined samples 
fit16 <- nmf(merged_mat, k = 16, verbose = FALSE)

# re-aligning row and column names (lost before)
rownames(fit16$w) <- rownames(merged_mat)
colnames(fit16$h) <- colnames(merged_mat)

# renaming for ease of reference
W <- fit16$w
H <- fit16$h  # do I need to scale this?


# -----------------------------------------------------------------------------
# 6. finalizing, binding spatial data
pd1_vis_coords <- readRDS(datadir, "pd1_vis_coords.rds")
iso_vis_coords <- readRDS(datadir, "iso_vis_coords.rds")

rownames(pd1_vis_coords)  <- paste0("pd1_", rownames(pd1_vis_coords))
rownames(iso_vis_coords) <- paste0("iso_", rownames(iso_vis_coords))

tH <- t(H)   # required to transpose to line up

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

# scale from pixels to microns for euclidean spatial alignment
pd1_joined$x <- pd1_joined$x / 1.5454                   # careful with overwriting.
pd1_joined$y <- pd1_joined$y / 1.5454

# export
write.csv(iso_joined, file.path(datadir, "0525_NMF_iso.csv"))

ds_merged <- merge(pd1_joined, iso_joined)              # does it make sense? separate model
                                                        # or together?
