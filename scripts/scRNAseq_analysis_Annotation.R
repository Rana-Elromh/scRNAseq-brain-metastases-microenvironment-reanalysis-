# =========================================================
# 1. Load packages
# =========================================================
library(ggplot2)
library(Seurat)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
library(SingleR)
library(celldex)

set.seed(42)

# =========================================================
# 2. Define samples and data directory
# =========================================================
samples <- list(
  BRBMET2  = "GSM7475324_BRBMET2",
  BRBMET3  = "GSM7475325_BRBMET3",
  BRBMET87 = "GSM7475326_BRBMET87",
  LUBMET7  = "GSM7475327_LUBMET7",
  LUBMET1  = "GSM7475328_LUBMET1"
)

data_dir <- "/Users/dell/scRNAseq-brain-metastases-microenvironment-reanalysis-/data/"


# =========================================================
# 3. Load the 5 samples
# =========================================================
seurat_list <- lapply(names(samples), function(s) {
  
  sample_dir <- file.path(data_dir, s)
  
  mat <- ReadMtx(
    mtx = file.path(sample_dir, "matrix.mtx.gz"),
    features = file.path(sample_dir, "features.tsv.gz"),
    cells = file.path(sample_dir, "barcodes.tsv.gz"),
    feature.column = 2
  )
  
  CreateSeuratObject(
    counts = mat,
    project = s,
    min.cells = 0,
    min.features = 0
  )
})

names(seurat_list) <- names(samples)


# =========================================================
# 4. QC for each sample (adaptive threshold capped at 40%)
# =========================================================

filtered_list <- lapply(names(seurat_list), function(s) {
  
  obj <- seurat_list[[s]]
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  meta <- obj@meta.data
  
  mt_upper <- min(
    median(meta$percent.mt) + 3 * mad(meta$percent.mt),
    40
  )
  
  feat_lower <- pmax(
    200,
    median(meta$nFeature_RNA) - 3 * mad(meta$nFeature_RNA)
  )
  
  feat_upper <- median(meta$nFeature_RNA) + 3 * mad(meta$nFeature_RNA)
  
  keep <- meta$percent.mt < mt_upper &
    meta$nFeature_RNA > feat_lower &
    meta$nFeature_RNA < feat_upper
  
  obj[, keep]
})

names(filtered_list) <- names(seurat_list)

lapply(filtered_list, ncol)


# =========================================================
# 5. Merge the 5 QC-filtered samples
# =========================================================

merged_qc <- merge(
  x = filtered_list[[1]],
  y = filtered_list[-1],
  add.cell.ids = names(filtered_list)
)

dim(merged_qc)

merged_qc <- JoinLayers(merged_qc)


# =========================================================
# 6. Detect doublets using scDblFinder
# =========================================================

sce <- as.SingleCellExperiment(merged_qc)

set.seed(100)

sce <- scDblFinder(
  sce,
  samples = "orig.ident"
)


# =========================================================
# 7. Add doublet results to Seurat object
# =========================================================

merged_qc$doublet_score <- colData(sce)$scDblFinder.score
merged_qc$doublet_class <- colData(sce)$scDblFinder.class

table(merged_qc$doublet_class)


# =========================================================
# 8. Remove doublets
# =========================================================

merged_singlets <- subset(
  merged_qc,
  subset = doublet_class == "singlet"
)

ncol(merged_singlets)


# =========================================================
# 9. Normalization
# =========================================================

merged_singlets <- NormalizeData(
  merged_singlets,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# =========================================================
# 10. Find Variable Features
# =========================================================

merged_singlets <- FindVariableFeatures(
  merged_singlets,
  selection.method = "vst",
  nfeatures = 2000
)

VariableFeaturePlot(merged_singlets)

head(VariableFeatures(merged_singlets), 10)

plot1 <- VariableFeaturePlot(merged_singlets)

top10 <- head(VariableFeatures(merged_singlets), 10)

LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE
)



# =========================================================
# 11. Scale Data
# =========================================================

merged_singlets <- ScaleData(merged_singlets)


# =========================================================
# 12. PCA
# =========================================================

merged_singlets <- RunPCA(
  merged_singlets,
  features = VariableFeatures(merged_singlets)
)

print(merged_singlets[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(
  merged_singlets,
  dims = 1:2,
  reduction = "pca"
)

DimPlot(
  merged_singlets,
  reduction = "pca"
)

DimHeatmap(
  merged_singlets,
  dims = 1,
  cells = 500,
  balanced = TRUE
)

ElbowPlot(
  merged_singlets,
  ndims = 50
)


# =========================================================
# 13. Harmony Integration
# =========================================================

merged_singlets <- RunHarmony(
  merged_singlets,
  group.by.vars = "orig.ident",
  dims.use = 1:30,
  theta = 4
)


# =========================================================
# 14. Clustering & UMAP BEFORE Harmony
# =========================================================

merged_singlets <- FindNeighbors(
  merged_singlets,
  reduction = "pca",
  dims = 1:30
)

merged_singlets <- FindClusters(
  merged_singlets,
  resolution = 0.5
)

merged_singlets <- RunUMAP(
  merged_singlets,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap_pca"
)

DimPlot(
  merged_singlets,
  reduction = "umap_pca",
  group.by = "orig.ident"
)


# =========================================================
# 15. Clustering & UMAP AFTER Harmony
# =========================================================

merged_singlets <- FindNeighbors(
  merged_singlets,
  reduction = "harmony",
  dims = 1:30
)

merged_singlets <- FindClusters(
  merged_singlets,
  resolution = 0.5
)

merged_singlets <- RunUMAP(
  merged_singlets,
  reduction = "harmony",
  dims = 1:30
)

DimPlot(
  merged_singlets,
  reduction = "umap",
  label = TRUE
)

DimPlot(
  merged_singlets,
  reduction = "umap",
  group.by = "orig.ident"
)

DimPlot(
  merged_singlets,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE
)
# =========================================================
# 16. t-SNE after Harmony integration
# =========================================================

merged_singlets <- RunTSNE(
  merged_singlets,
  reduction = "harmony",
  dims = 1:30,
  reduction.name = "tsne_harmony"
)

DimPlot(
  merged_singlets,
  reduction = "tsne_harmony",
  group.by = "seurat_clusters",
  label = TRUE
)

DimPlot(
  merged_singlets,
  reduction = "tsne_harmony",
  group.by = "orig.ident"
)


# =========================================================
# 17. Find marker genes for each cluster
# =========================================================

markers <- FindAllMarkers(
  merged_singlets,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

top_markers

print(top_markers, n = 80)

# =========================================================
# 18. Subclustering: Myeloid cluster
# =========================================================

myeloid_markers_check <- c("TREM2", "C1QB", "GPR34", "FOLR2", "C1QC")

avg_expr <- AverageExpression(
  merged_singlets,
  features = myeloid_markers_check,
  group.by = "seurat_clusters"
)$RNA

myeloid_score <- colMeans(avg_expr)

myeloid_cluster_id <- names(which.max(myeloid_score))
myeloid_cluster_id <- "0"
cat("Candidate myeloid cluster:", myeloid_cluster_id, "\n")
print(sort(myeloid_score, decreasing = TRUE))

myeloid_cells <- WhichCells(
  merged_singlets,
  expression = seurat_clusters == myeloid_cluster_id
)

myeloid <- subset(
  merged_singlets,
  cells = myeloid_cells
)

ncol(myeloid)

table(myeloid$orig.ident)


myeloid <- FindVariableFeatures(
  myeloid,
  selection.method = "vst",
  nfeatures = 2000
)

myeloid <- ScaleData(myeloid)

myeloid <- RunPCA(
  myeloid,
  features = VariableFeatures(myeloid)
)


ElbowPlot(
  myeloid,
  ndims = 30
)


n_dims <- 1:20


myeloid <- RunHarmony(
  myeloid,
  group.by.vars = "orig.ident",
  dims.use = n_dims,
  theta = 4
)


myeloid <- FindNeighbors(
  myeloid,
  reduction = "harmony",
  dims = n_dims
)

myeloid <- FindClusters(
  myeloid,
  resolution = 0.2
)

myeloid <- RunUMAP(
  myeloid,
  reduction = "harmony",
  dims = n_dims
)


DimPlot(
  myeloid,
  reduction = "umap",
  label = TRUE
)

DimPlot(
  myeloid,
  reduction = "umap",
  group.by = "orig.ident"
)


# =========================================================
# 19. Find marker genes for each myeloid subcluster
# =========================================================

myeloid_markers <- FindAllMarkers(
  myeloid,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top_myeloid_markers <- myeloid_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 5
  )

print(top_myeloid_markers, n = Inf)


FeaturePlot(
  myeloid,
  features = c(
    "C1QC",
    "C1QB",
    "TREM2",
    "FOLR2",
    "LYZ",
    "CD3D",
    "CD3E",
    "NKG7",
    "GNLY",
    "MKI67"
  ),
  reduction = "umap"
)

#=============================================================================
# 20. Automated Annotation of merged_singlets
#=============================================================================
BiocManager::install(c("SingleR", "celldex"))
library(SingleR)
library(celldex)

hpca_ref   <- celldex::HumanPrimaryCellAtlasData()
monaco_ref <- celldex::MonacoImmuneData()

SCE_merged_singlets <- as.SingleCellExperiment(merged_singlets)

singleR_merged_singlets <- SingleR(
  test   = SCE_merged_singlets,
  ref    = list(HPCA = hpca_ref, Monaco = monaco_ref),
  labels = list(hpca_ref$label.main, monaco_ref$label.main)
)

merged_singlets$singleR_annotation <- singleR_merged_singlets$labels

DimPlot(
  merged_singlets,
  reduction = "umap",
  group.by = "singleR_annotation",
  label = TRUE,
  repel = TRUE,
  label.size = 3.5
) +
  ggtitle("Merged Singlets - SingleR Automated Annotation")


#==============================================================================
# Differential expression annotation of merged-singlets
#==============================================================================
markers_merged_singlets <- FindAllMarkers(
  merged_singlets, only.pos = TRUE,
  min.pct = 0.25, logfc.threshold = 0.25
)

top_markers_merged <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 25)

top_markers_merged

write.csv(top_markers_merged, "top_markers_merged.csv", row.names = FALSE)


cluster_labels <- c(
  "0"  = "Myloid",
  "1"  = "Tumor_basal_like",
  "2"  = "Tumor_luminal_progenitor_like",
  "3"  = "Tumor_proliferating",
  "4"  = "Tumor_hormone_receptor_luminal",
  "5"  = "Tumor_secretory",
  "6"  = "T_NK_cell",
  "7"  = "Oligodendrocyte",
  "8"  = "Tumor_HER2_luminal",
  "9"  = "Ependymal",
  "10" = "Tumor_hypoxic",
  "11" = "Pericyte_fibroblast",
  "12" = "Endothelial",
  "13" = "Monocyte_neutrophil",
  "14" = "B_plasma_cell"
)

setdiff(levels(merged_singlets$seurat_clusters), names(cluster_labels))  # should return character(0)

merged_singlets$cell_type <- unname(cluster_labels[as.character(merged_singlets$seurat_clusters)])
Idents(merged_singlets) <- "cell_type"

DimPlot(merged_singlets, reduction = "umap", group.by = "cell_type",
        label = TRUE, repel = TRUE) + ggtitle("Cell type annotation")

#===============================================================================
#Differential expression annotation of Myeloid cells
#===============================================================================
top_myeloid_markers <- myeloid_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 25)

write.csv(top_myeloid_markers, "top25_myeloid_markers.csv", row.names = FALSE)


myeloid_labels <- c(
  "0" = "Macrophage",
  "1" = "Macrophage_MMP9_high_mito",
  "2" = "Microglia",
  "3" = "cDC2",
  "4" = "Proliferating_myeloid"
)

myeloid$myeloid_subtype <- unname(myeloid_labels[as.character(myeloid$seurat_clusters)])
Idents(myeloid) <- "myeloid_subtype"

DimPlot(myeloid, reduction = "umap", label = TRUE, repel = TRUE) +
  ggtitle("Myeloid subtypes and states")

#===============================================================================
# Determining of M1 and M2 states
#=========================================================================
m1_genes <- c("CD80","CD86","TNF","IL1B","IL6","IL12A","IL12B","CXCL9","CXCL10",
              "CXCL11","NOS2","STAT1","IRF5","FCGR1A","IDO1","CCL5","HLA-DRA")
m2_genes <- c("CD163","MRC1","MSR1","IL10","TGFB1","ARG1","CCL18","CCL22",
              "CD209","VEGFA","CLEC7A","TGM2","STAB1","MARCO")

m1_genes <- intersect(m1_genes, rownames(myeloid))
m2_genes <- intersect(m2_genes, rownames(myeloid))

myeloid <- AddModuleScore(myeloid, features = list(m1_genes), name = "M1_score")
myeloid <- AddModuleScore(myeloid, features = list(m2_genes), name = "M2_score")

# compare clusters
VlnPlot(myeloid, features = c("M1_score1", "M2_score1"),
        group.by = "myeloid_subtype", pt.size = 0)

myeloid@meta.data %>%
  group_by(myeloid_subtype) %>%
  summarise(M1 = mean(M1_score1), M2 = mean(M2_score1), diff = M1 - M2)
#===============================================================================
# Annotation of M1 and M2 states
#===============================================================================
myeloid$cell_type_mm <- dplyr::case_when(
  myeloid$myeloid_subtype == "Microglia" ~ "Microglia",
  myeloid$myeloid_subtype %in% c("Macrophage", "Macrophage_MMP9_high_mito") ~ "Macrophage",
  TRUE ~ as.character(myeloid$myeloid_subtype)          
)

myeloid$M1_minus_M2 <- myeloid$M1_score1 - myeloid$M2_score1
cut <- 0.1   

myeloid$state <- ifelse(myeloid$M1_minus_M2 >  cut, "M1-like",
                        ifelse(myeloid$M1_minus_M2 < -cut, "M2-like", "Intermediate"))


myeloid$annotation_state <- ifelse(
  myeloid$cell_type_mm %in% c("Macrophage", "Microglia"),
  paste(myeloid$cell_type_mm, myeloid$state, sep = "_"),
  myeloid$cell_type_mm
)

table(myeloid$annotation_state)

#Plot
DimPlot(myeloid, reduction = "umap", group.by = "annotation_state",
        label = TRUE, repel = TRUE) +
  ggtitle("Myeloid cell types and M1/M2-like states")
#=============================================================================
# Different visualizations of M1 and M2 states
#=============================================================================

# Plot 1
DimPlot(myeloid, reduction = "umap", group.by = "cell_type_mm",
        split.by = "state", label = FALSE) +
  ggtitle("Cell types split by M1/M2-like state")

# Plot 2 (Blue cells lean M1-like and red cells lean M2-like.)

FeaturePlot(myeloid, features = "M1_minus_M2", reduction = "umap") +
  scale_color_gradient2(low = "red", mid = "grey90", high = "blue", midpoint = 0) +
  ggtitle("M1 score minus M2 score")

# Plot 3

mm_cells <- myeloid@meta.data %>%
  filter(cell_type_mm %in% c("Macrophage", "Microglia"))

ggplot(mm_cells, aes(x = cell_type_mm, fill = state)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "State") +
  theme_classic()


round(prop.table(table(mm_cells$cell_type_mm, mm_cells$state), margin = 1) * 100, 1)

# Plot 4 

DimPlot(myeloid, reduction = "umap", group.by = "cell_type_mm",
        label = TRUE, repel = TRUE) + ggtitle("Myeloid cell types")












