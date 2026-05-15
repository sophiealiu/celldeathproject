# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/06/2026
# Purpose: Generation of gene scores based on known genes and signatures.
# -----------------------------------------------------------------------------

# Loading in necessary libraries & directories 
library(arrow)
library(dplyr)
library(ggplot2)
library(limma)
library(presto)      # Seurat V5 incompatible, forced in using remotes + Github
library(Seurat)
library(SeuratData)  # same as above


inputdir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"

pd1_norm <- readRDS(file.path(inputdir, "Seurat objects", "pd1_norm.rds"))
            # generated using Load10x function then log normalization
iso7_norm <- readRDS(file.path(inputdir, "Seurat objects", "iso7_norm.rds"))


# -----------------------------------------------------------------------------
# List of known signatures
cd8_effector <- c("Cd8a","Cd8b1","Gzmb","Gzmk","Prf1","Ifng",
                  "Nkg7","Cx3cr1","Klrk1","Ccl5","Xcl1")    # Cx3cr1 missing from our dataset
cd8_effector <- setdiff(cd8_effector, "Cx3cr1")

cd8_exhausted <- c("Pdcd1","Ctla4","Lag3","Havcr2","Tigit",
                   "Tox","Eomes","Cxcl13","Batf","Nr4a1")
cd4_effector <- c("Cd4","Cd40lg","Il2","Ifng", "Tbx21","Cxcr3","Ccl5")

treg <- c("Foxp3","Il2ra","Ctla4","Ikzf2","Tnfrsf18","Tnfrsf4","Entpd1")

nk <- c("Ncr1","Klrk1","Klrd1","Nkg7","Gzmb","Prf1","Ifng","Xcl1")

folr2_mac <- c( "Folr2","Lyve1","Mrc1","Sepp1",
                "Gas6","Timd4")                             # Sepp1 missing ""
folr2_mac <- setdiff(folr2_mac, "Sepp1")

monocytes <- c("Ly6c2","Ccr2","Sell","S100a8","S100a9",
               "Plac8","Lyz2")                              # Ly6c2 missing ""
monocytes <- setdiff(monocytes, "Ly6c2")

cdc1 <- c("Xcr1","Clec9a","Batf3","Irf8",
          "Cadm1","Itgae")
cdc2 <- c("Cd209a","Sirpa","Itgam","Irf4", "Clec10a","Fcgr2b")

endothelial <- c("Pecam1","Cdh5","Kdr","Flt1","Vwf","Emcn")

fibroblast <- c("Col1a1","Col1a2","Col3a1","Pdgfra",
                "Dcn","Lum","Col5a1","Col6a1")

# aggregating
signatures <- list(
  cd8_effector = cd8_effector,
  cd8_exhausted = cd8_exhausted,
  cd4_effector = cd4_effector,
  treg = treg,
  nk = nk,
  folr2_mac = folr2_mac,
  monocytes = monocytes,
  cdc1 = cdc1,
  cdc2 = cdc2,
  endothelial = endothelial,
  fibroblast = fibroblast
)


# -----------------------------------------------------------------------------
# 1. Binding expression per bin with bin coordinates
iso7_coords <- readRDS(file.path(inputdir, "Seurat objects", "iso7bin_coords.rds"))
    # from GetTissueCoordinates function, image = NULL

# expression matrix (all the genes by expression)
expr_matIso <- GetAssayData(iso7_norm,
                         assay = "Spatial.008um",
                         layer = "data")

# optional explicitly making sure there's no missing genes
all_genes <- rownames(expr_matIso)

signatures <- lapply(signatures, function(genes) {
  present <- intersect(genes, all_genes)
  if (length(present) == 0) {
    warning("No genes found for signature: ", paste(genes, collapse=", "))
  }
  present
})

# assigning the gene expression to each bin
sig_Iso <- do.call(rbind, lapply(signatures, function(genes) {
  colMeans(expr_matIso[genes, , drop = FALSE], na.rm = TRUE)
}))

# required to transpose to match row and column names
sig_Iso <- t(sig_Iso) %>%
  as.data.frame()

# assigning coordinates to bins
scoresIso <- cbind(iso7_coords, sig_Iso)

# transforming from pixels to microns. easier to do here/upstream
scoresIso$x <- scoresIso$x / 1.5454 
scoresIso$y <- scoresIso$y / 1.5454


# -----------------------------------------------------------------------------
# 2. Verify visually proper spatial location
df_iso7IF <- read.csv(file.path(inputdir, "iso7_coords_clean.csv"))
                # generated from Imaris

check <- function(df_IF, df_visium, trans) { 
  ggplot() +
    stat_density_2d(
      data = filter(df_IF, cell_type == "tdtomato"), # density, too many tdtomato
      aes(x = x, y = y),
      fill        = "red",
      geom        = "density_2d_filled",
      contour_var = "ndensity",
      alpha       = 0.10,
      n           = 400,
      na.rm       = TRUE
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "cd8"),    # points, clarity
      aes(x = x, y = y),
      color = "#02819e",
      alpha = 0.20,
      size  = 0.01
    ) +
    geom_point(
      data  = filter(df_IF, cell_type == "gc3ai"),
      aes(x = x, y = y),
      color = "green4",
      size  = 0.0005, 
      alpha = 0.15
    ) +  
    
    geom_point(
      data = df_visium,                             # overlay
      aes(x = x, y = y),
      color = "white",
      alpha = trans,
      size = 0.01
    ) +
    
# formatting
    scale_x_continuous(expand = c(0.05, 0.05)) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    coord_equal(clip = "off") +
    theme(
      panel.background = element_rect(fill = "black"),
      plot.background  = element_rect(fill = "black"),
      panel.grid       = element_blank(),
      
      axis.text.x  = element_text(color = "white", size = 12),
      axis.text.y  = element_text(color = "white", size = 12),
      
      axis.line = element_line(color = "white")
    )  
}

stopifnot(all(c("x", "y") %in% colnames(scoresIso)))
stopifnot(all(c("x", "y", "cell_type") %in% colnames(df_iso7IF)))

test <- check(df_iso7IF, scoresIso, 0.05)
test # make sure it looks okay before running below


# export
write.csv(scoresIso, file.path(inputdir, "iso7spatial_sig.csv"))


# -----------------------------------------------------------------------------
# 3. Using all genes instead of known signatures. Same logic.
DEgenes <- read.csv(file.path(inputdir, "topDE.csv"))
genelist <- DEgenes$gene

# expression to bin
sig_Iso <- do.call(rbind, lapply(genelist, function(genes) {
  colMeans(expr_matIso[genes, , drop = FALSE], na.rm = TRUE)
}))

# transpose
sig_Iso <- t(sig_Iso) %>%
  as.data.frame()

# coords to bins
scoresIso <- cbind(iso7_coords, sig_Iso)

# transform
scoresIso$x <- scoresIso$x / 1.5454 
scoresIso$y <- scoresIso$y / 1.5454

# export
write.csv(scoresIso, file.path(inputdir, "iso7spatial_mark.csv"))
