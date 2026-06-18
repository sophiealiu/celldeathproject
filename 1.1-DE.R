# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Date : 05/28/2026 11:31am CDT
# Purpose: preliminary looking at individual genes to known signatures differential expression
# -----------------------------------------------------------------------------

# necessary libraries & directories 
library(dplyr)
library(ggplot2)
library(Seurat)

# recommended 8 micron bins
pd1 <- Load10X_Spatial("Iso-7_realign/outs", bin.size = 8)
iso <- Load10X_Spatial("PD1-9_realign/outs", bin.size = 8)

# log normalized with unique IDs
iso_norm <- NormalizeData(iso)
pd1_norm <- NormalizeData(pd1)
colnames(iso_norm) <- paste0("iso_", colnames(iso_norm))
colnames(pd1_norm) <- paste0("pd1_", colnames(pd1_norm))

merged_norm <- merge(iso_norm, pd1_norm)


# -----------------------------------------------------------------------------
# 2. signature finding. repeat with experimental replacing object
pd1_vis_coords <- GetTissueCoordinates(pd1_norm, image = NULL)
iso_vis_coords <- GetTissueCoordinates(iso_norm, image = NULL)

# expression matrix (all the genes by expression) 
merged_norm <- JoinLayers(
  merged_norm
)

merged_matN <- GetAssayData(
  merged_norm,
  layer = "data"
)

tmerged_matN <- t(merged_matN)
rownames(tmerged_matN) <- colnames(GetAssayData(merged_norm,  
                                        layer = "data"))

common_pd1 <- rownames(pd1_vis_coords)[
  rownames(pd1_vis_coords) %in% colnames(merged_matN)
]
common_iso <- rownames(iso_vis_coords)[
  rownames(iso_vis_coords) %in% colnames(merged_matN)
]

# must have large RAM to allocate vector size 58.8 Gb
pd1_joined <- cbind(
  pd1_vis_coords[common_pd1, , drop = FALSE],
  tmerged_matN[common_pd1, , drop = FALSE]
)

iso_joined <- cbind(
  iso_vis_coords[common_iso, , drop = FALSE],
  tmerged_matN[common_iso, , drop = FALSE]
)

# writing all genes to csv takes forever. repeat with both conditions
library(arrow) 
write_feather(iso_joined, "iso_all_genes")
write_feather(pd1_joined, "pd1_all_genes")


# -----------------------------------------------------------------------------
# 2. cell-type known signatures 
cd8_effector <- c("Cd8a","Cd8b1","Gzmb","Gzmk","Prf1","Ifng",
                  "Nkg7","Cx3cr1","Klrk1","Ccl5","Xcl1")     # Cx3cr1 missing
cd8_effector <- setdiff(cd8_effector, "Cx3cr1")

cd8_exhausted <- c("Pdcd1","Ctla4","Lag3","Havcr2","Tigit",
                   "Tox","Eomes","Cxcl13","Batf","Nr4a1")
cd4_effector <- c("Cd4","Cd40lg","Il2","Ifng", "Tbx21","Cxcr3","Ccl5")

treg <- c("Foxp3","Il2ra","Ctla4","Ikzf2","Tnfrsf18","Tnfrsf4","Entpd1")

nk <- c("Ncr1","Klrk1","Klrd1","Nkg7","Gzmb","Prf1","Ifng","Xcl1")

folr2_mac <- c( "Folr2","Lyve1","Mrc1","Sepp1",
                "Gas6","Timd4")                             # Sepp1 missing
folr2_mac <- setdiff(folr2_mac, "Sepp1")

monocytes <- c("Ly6c2","Ccr2","Sell","S100a8","S100a9",
               "Plac8","Lyz2")                              # Ly6c2 missing
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

# assigning the signature expression to each bin
sig_all <- do.call(rbind, lapply(signatures, function(genes) {
  colMeans(merged_matN[genes, , drop = FALSE], na.rm = TRUE)
}))

# transpose to line up
tsig_all <- t(sig_all) %>%
  as.data.frame()

common_pd12 <- rownames(pd1_vis_coords)[
  rownames(pd1_vis_coords) %in% colnames(sig_all)
]
common_iso2 <- rownames(iso_vis_coords)[
  rownames(iso_vis_coords) %in% colnames(sig_all)
]

pd1_joined2 <- cbind(
  pd1_vis_coords[common_pd12, , drop = FALSE],
  tsig_all[common_pd12, , drop = FALSE]
)

iso_joined2 <- cbind(
  iso_vis_coords[common_iso2, , drop = FALSE],
  tsig_all[common_iso2, , drop = FALSE]
)

# exporting for jupyter block 2
write.csv(iso_joined2, "iso_sig.csv")
write.csv(pd1_joined2, "pd1_sig.csv")


# -----------------------------------------------------------------------------
# 4. DE with FindMarkers on annotated cells. viewing metadata first
celltypes <- levels(Idents(merged_annotated))     

de.mark_known <- list()

for (ct in celltypes) {
  de.mark_known[[ct]] <- FindMarkers(
    merged_annotated,
    ident.1 = ct
  )
}

# exporting
head(rownames(de.mark_known[["Fibroblasts"]]), 10) # inspection
write.csv(de.mark_known, "mark_known.csv")



