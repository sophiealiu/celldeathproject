# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/11/2026
# Purpose: Generation of gene scores based on known genes and signatures.
# -----------------------------------------------------------------------------

# Loading in necessary libraries & directories 
library(arrow)
library(dplyr)
library(ggplot2)
library(limma)
library(presto)      # Seurat V5 incompatible, forced in using remotes
library(Seurat)
library(SeuratData)  # same as above


inputdir <- "I:/Hu Lab/Sophie/1. Cell death/visium image manual spot selection/20260413_final_merge/data"

pd1_norm <- readRDS(file.path(inputdir, "Seurat objects", "pd1_norm.rds"))
            # generated using Load10x function then log normalization
iso7_norm <- readRDS(file.path(inputdir, "Seurat objects", "iso7_norm.rds"))


# -----------------------------------------------------------------------------
# 1. List of known signatures
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
  cd8_effector,
  cd8_exhausted,
  cd4_effector,
  treg,
  nk,
  folr2_mac,
  monocytes,
  cdc1,
  cdc2,
  endothelial,
  fibroblast
)

# -----------------------------------------------------------------------------
# 2. FindMarkers on manually annotated merged MC38s
celltypes <- levels(Idents(merged_annotated))       # inspecting metadata first

de.mark_known <- list()

for (ct in celltypes) {
  de.mark_known[[ct]] <- FindMarkers(
    merged_annotated,
    ident.1 = ct
  )
}

saveRDS(de.mark_known, file.path(inputdir, "Seurat objects", "mark_known.rds"))
head(rownames(de.mark_known[["Fibroblasts"]]), 10)   


# -----------------------------------------------------------------------------
# 3. binding spatially 
iso7_coords <- readRDS(file.path(inputdir, "Seurat objects", "iso7bin_coords.rds"))
    # from GetTissueCoordinates function, image = NULL

# expression matrix (all the genes by expression)
expr_matIso <- GetAssayData(iso7_norm,
                         assay = "Spatial.008um",
                         layer = "data")

# optional explicitly making sure there's no missing genes
all_genes <- rownames(expr_matIso)


# 1a. Obtaining known signature expression per bin
signatures <- lapply(signatures, function(genes) {   # for 2, sub mark_known for signatures
  present <- intersect(genes, all_genes)
  if (length(present) == 0) {
    warning("No genes found for signature: ", paste(genes, collapse=", "))
  }
  present
})

# 1b. Using all genes instead of known signatures. Same logic.
DEgenes <- read.csv(file.path(inputdir, "topDE.csv"))
genelist <- DEgenes$gene


# -----------------------------------------------------------------------------
# 2. Verify proper spatial location 
plot(scoresIso, x = x, y = y)


