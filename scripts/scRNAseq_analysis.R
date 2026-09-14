#===========================================
# SAMPLE 1 . BRBMET2
#===========================================
library(Seurat)

# Read BRBMET2 data
sc_data <- Read10X(
  data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/BRBMET2/"
)

# Create Seurat object
BRBMET2_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET2"
)

# Detect doublets
library(scDblFinder)
library(SingleCellExperiment)

# Convert Seurat object to SingleCellExperiment
sce_d1 <- as.SingleCellExperiment(BRBMET2_seurat)

# Run scDblFinder
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

# Add doublet results to Seurat object
BRBMET2_seurat$doublet_score <- colData(sce_d1)$scDblFinder.score
BRBMET2_seurat$doublet_class <- colData(sce_d1)$scDblFinder.class

# Check doublet classification
table(BRBMET2_seurat$doublet_class)

# Remove doublets
BRBMET2_seurat <- subset(
  BRBMET2_seurat,
  subset = doublet_class == "singlet"
)

# Check number of cells after removing doublets
ncol(BRBMET2_seurat)
#===========================================
# SAMPLE 2 . BRBMET3
#===========================================
library(Seurat)

# Read BRBMET3 data
sc_data <- Read10X(
  data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/BRBMET3/"
)

# Create Seurat object
BRBMET3_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET3"
)

# Detect doublets
library(scDblFinder)
library(SingleCellExperiment)

# Convert Seurat object to SingleCellExperiment
sce_d2 <- as.SingleCellExperiment(BRBMET3_seurat)

# Run scDblFinder
set.seed(100)
sce_d2 <- scDblFinder(sce_d2)

# Add doublet results to Seurat object
BRBMET3_seurat$doublet_score <- colData(sce_d2)$scDblFinder.score
BRBMET3_seurat$doublet_class <- colData(sce_d2)$scDblFinder.class

# Check doublet classification
table(BRBMET3_seurat$doublet_class)

# Remove doublets
BRBMET3_seurat <- subset(
  BRBMET3_seurat,
  subset = doublet_class == "singlet"
)

# Check number of cells after removing doublets
ncol(BRBMET3_seurat)
#===========================================
# SAMPLE 3 . BRBMET87
#===========================================
library(Seurat)

# Read BRBMET87 data
sc_data <- Read10X(
  data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/BRBMET87/"
)

# Create Seurat object
BRBMET87_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET87"
)

# Detect doublets
library(scDblFinder)
library(SingleCellExperiment)

# Convert Seurat object to SingleCellExperiment
sce_d3 <- as.SingleCellExperiment(BRBMET87_seurat)

# Run scDblFinder
set.seed(100)
sce_d3 <- scDblFinder(sce_d3)

# Add doublet results to Seurat object
BRBMET87_seurat$doublet_score <- colData(sce_d3)$scDblFinder.score
BRBMET87_seurat$doublet_class <- colData(sce_d3)$scDblFinder.class

# Check doublet classification
table(BRBMET87_seurat$doublet_class)

# Remove doublets
BRBMET87_seurat <- subset(
  BRBMET87_seurat,
  subset = doublet_class == "singlet"
)

# Check number of cells after removing doublets
ncol(BRBMET87_seurat)
#===========================================
# SAMPLE 4 . LUBMET1
#===========================================
library(Seurat)

# Read LUBMET1 data
sc_data <- Read10X(
  data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/LUBMET1/"
)

# Create Seurat object
LUBMET1_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "LUBMET1"
)

# Detect doublets
library(scDblFinder)
library(SingleCellExperiment)

# Convert Seurat object to SingleCellExperiment
sce_d5 <- as.SingleCellExperiment(LUBMET1_seurat)

# Run scDblFinder
set.seed(100)
sce_d5 <- scDblFinder(sce_d5)

# Add doublet results to Seurat object
LUBMET1_seurat$doublet_score <- colData(sce_d5)$scDblFinder.score
LUBMET1_seurat$doublet_class <- colData(sce_d5)$scDblFinder.class

# Check doublet classification
table(LUBMET1_seurat$doublet_class)

# Remove doublets
LUBMET1_seurat <- subset(
  LUBMET1_seurat,
  subset = doublet_class == "singlet"
)

# Check number of cells after removing doublets
ncol(LUBMET1_seurat)
#===========================================
# SAMPLE 5 . LUBMET7
#===========================================
library(Seurat)

# Read LUBMET7 data
sc_data <- Read10X(
  data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/LUBMET7/"
)

# Create Seurat object
LUBMET7_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "LUBMET7"
)

# Detect doublets
library(scDblFinder)
library(SingleCellExperiment)

# Convert Seurat object to SingleCellExperiment
sce_d4 <- as.SingleCellExperiment(LUBMET7_seurat)

# Run scDblFinder
set.seed(100)
sce_d4 <- scDblFinder(sce_d4)

# Add doublet results to Seurat object
LUBMET7_seurat$doublet_score <- colData(sce_d4)$scDblFinder.score
LUBMET7_seurat$doublet_class <- colData(sce_d4)$scDblFinder.class

# Check doublet classification
table(LUBMET7_seurat$doublet_class)

# Remove doublets
LUBMET7_seurat <- subset(
  LUBMET7_seurat,
  subset = doublet_class == "singlet"
)

# Check number of cells after removing doublets
ncol(LUBMET7_seurat)
#===========================================
# Merge the five samples
merged_seurat <- merge(
  BRBMET2_seurat,
  y = list(BRBMET3_seurat, BRBMET87_seurat, LUBMET7_seurat, LUBMET1_seurat),
  add.cell.ids = c("BRBMET2", "BRBMET3", "BRBMET87", "LUBMET7", "LUBMET1")
)

# Check the total number of cells
ncol(merged_seurat)

# Normalize the data
merged_seurat <- NormalizeData(
  merged_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Identify highly variable genes
merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

# Visualize variable features
VariableFeaturePlot(merged_seurat)

# Display the top 10 variable genes
head(VariableFeatures(merged_seurat), 10)

plot1 <- VariableFeaturePlot(merged_seurat)

top10 <- head(VariableFeatures(merged_seurat), 10)

LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE
)

# Scale the variable features
merged_seurat <- ScaleData(
  merged_seurat,
  features = VariableFeatures(merged_seurat)
)

# Run PCA on the scaled data
merged_seurat <- RunPCA(
  merged_seurat,
  features = VariableFeatures(merged_seurat)
)
print(merged_seurat[["pca"]], dims = 1:5, nfeatures = 5)
# Visualize loadings for the top PCs
VizDimLoadings(merged_seurat, dims = 1:2, reduction = "pca")

# PCA scatter plot (PC1 vs PC2)
DimPlot(merged_seurat, reduction = "pca")

# Heatmap to explore heterogeneity within a PC
DimHeatmap(merged_seurat, dims = 1, cells = 500, balanced = TRUE)

# Determine how many PCs to use downstream (elbow plot)
ElbowPlot(merged_seurat, ndims = 50)
# ---------------------------
#integration across samples using Harmony
#---------------------------------------------
library(harmony)
merged_seurat <- RunHarmony(
  merged_seurat,
  group.by.vars = "orig.ident",
  dims.use = 1:30,
  theta = 4  )

# --------------------------------------
# Step: Clustering & UMAP
# (Finding cell groups + 2D visualization based on top 30 PCs)
# -------------------------------------------------------------
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
DimPlot(merged_seurat, reduction = "umap", label = TRUE)

DimPlot(merged_seurat, reduction = "umap", group.by = "orig.ident")

merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])

#----------------------------------------------------------

# Find marker genes for each cluster
markers <- FindAllMarkers(
  merged_seurat,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Top markers per cluster
library(dplyr)
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

top_markers
----------------------------------------------------------------
#Sub-clustering the Fibroblasts population (cluster 7)
fibroblasts <- subset(merged_seurat, idents = "7")

fibroblasts <- FindVariableFeatures(fibroblasts)
fibroblasts <- ScaleData(fibroblasts)
fibroblasts <- RunPCA(fibroblasts)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:15)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.3)
fibroblasts <- RunUMAP(fibroblasts, dims = 1:15)

DimPlot(fibroblasts, label = TRUE)

fibro_markers <- FindAllMarkers(
  fibroblasts,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

fibro_top_markers <- fibro_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

fibro_top_markers
#---------------------------------------------------
