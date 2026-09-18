# =========================================================
# 1. Load packages
# =========================================================

library(Seurat)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
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

data_dir <- "/Users/Technology - Laptoop/Downloads/GSE234832_RAW"


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
  resolution = 0.5
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
