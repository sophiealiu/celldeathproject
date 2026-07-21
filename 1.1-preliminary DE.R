# -----------------------------------------------------------------------------
# Author: Sophie A. Liu
# Purpose: Exploring individual genes then moving to differential expression
#        - similar to existing Seurat vignettes
# -----------------------------------------------------------------------------

# necessary libraries & directories 
library(dplyr)
library(ggplot2)
library(Seurat)

# recommended 8 micron bins
# datadir <- "path/to/your/working/directory"
datadir <- "I:/Hu Lab/Sophie/1. Cell death/all final data"

iso_raw <- readRDS(file.path(datadir,"iso_raw.rds"))
pd1_raw <- readRDS(file.path(datadir,"pd1_raw.rds"))

# log normalized with unique IDs. following default Seurat settings
iso_norm <- NormalizeData(iso_raw)
pd1_norm <- NormalizeData(pd1_raw)
colnames(iso_norm) <- paste0("iso_", colnames(iso_norm))
colnames(pd1_norm) <- paste0("pd1_", colnames(pd1_norm))

merged_norm <- merge(iso_norm, pd1_norm)


# -----------------------------------------------------------------------------
# 1. individual gene expression, 19060 total
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

# writing all genes to csv takes forever so use feather.
library(arrow) 
write_feather(iso_joined, "iso_all_genes")
write_feather(pd1_joined, "pd1_all_genes")


# -----------------------------------------------------------------------------
# 2. cell-type known signatures from literature
cd8_effector <- c("Cd8a","Cd8b1","Gzmb","Gzmk","Prf1","Ifng",
                  "Nkg7","Cx3cr1","Klrk1","Ccl5","Xcl1")     # Cx3cr1 missing in our data.
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
write.csv(pd1_joined2, "apd1_sig.csv")


# -----------------------------------------------------------------------------
# 3. DE with FindMarkers on annotated cells. viewing metadata first
celltypes <- levels(Idents(merged_annotated))     

de.mark_known <- list()

for (ct in celltypes) {
  de.mark_known[[ct]] <- FindMarkers(
    merged_annotated,
    ident.1 = ct
  )
}

markers_all <- bind_rows(
  lapply(names(de.mark_known), function(ct) {
    df <- de.mark_known[[ct]]
    df$gene <- rownames(df)
    df$celltype <- ct
    df
  })
)

# UMAP
de.mark_known <- RunUMAP(de.mark_known, dims = 1:20)

DimPlot(object = de.mark_known, 
        reduction = "umap", 
        group.by = "finercelltype", 
        label = TRUE, 
        pt.size = 0.5)


# -----------------------------------------------------------------------------
# 4. clustering
Assays(pd1_raw)
DefaultAssay(pd1_raw) <- "Spatial.008um" # do we want higher or lower resolution?
pd1_norm <- NormalizeData(pd1_raw)

# unsupervised: ML, finding similar gene expression profiles w/o knowing cell type
DefaultAssay(pd1_norm) <- "Spatial.008um"
pd1_norm <- FindVariableFeatures(pd1_norm)
pd1_norm <- ScaleData(pd1_norm)

pd1_norm <- SketchData(
  object = pd1_norm,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# beginning adding more metadata
pd1_norm <- RunPCA(pd1_norm, assay = "sketch", reduction.name = "pca.sketch")
pd1_norm <- FindNeighbors(pd1_norm, assay = "sketch", reduction = "pca.sketch", dims = 1:10)
pd1_norm <- FindClusters(pd1_norm, cluster.name = "seurat_cluster.sketched", resolution = 1) 
        #Shared Nearest Neighbor (SNN) graph, cells are nodes and edges weighted 
        #based on shared local neighborhoods.
pd1_norm <- RunUMAP(pd1_norm, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

# ProjectData maps results from the sketch subset to the entire dataset
options(future.globals.maxSize = 2000 * 1024^2)  # 2 GB Because otherwise too big, whoopsie!
pd1_norm <- ProjectData(
  object = pd1_norm, # FUNCTION DEFINITION
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch", # WHAT IS THIS
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

# clusters by similarity & density. Figure 5, panel B
Idents(pd1_norm) <- "seurat_cluster.sketched"
p008 <- DimPlot(pd1_norm, reduction = "umap.sketch", label = F) 
ggtitle("Sketched clustering (50000 cells)") 
theme(legend.position = "bottom") 
coord_fixed()
p008

Idents(pd1_norm) <- "seurat_cluster.projected"
cells <- CellsByIdentities(pd1_norm, idents = c(0:15))
      # it appears that only clusters 2&13 have much that's interesting tbh
pCluster <- SpatialDimPlot(
  pd1_norm,
  cells.highlight = cells[setdiff(names(cells), "NA")],
  cols.highlight = c("#FF0000", "transparent"), 
  image.alpha = 0,
  facet.highlight = TRUE, 
  combine = TRUE # get them in the same panel for easier visualization...
)
pCluster


# appending coordinates
coords2 <- GetTissueCoordinates(pd1_norm)
coords2$cluster <- Idents(pd1_norm)[rownames(coords2)]
head(coords2)

write.csv(coords2, file.path(datadir,"apd1_clusters.csv"), row.names = FALSE)

# supplemental figure
unique <- unique(top7$gene)
DotPlot(
  pd1_norm,
  features = unique,
  cols = c("lightgrey", "blue"),
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA
)
