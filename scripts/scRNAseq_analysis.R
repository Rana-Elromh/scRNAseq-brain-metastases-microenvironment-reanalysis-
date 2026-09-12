#===========================================
library(Seurat)
sc_data <- Read10X(data.dir = "D:/oza/GSE234832_RAW/BRBMET2/")

BRBMET2_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET2"
)
BRBMET2_seurat[["percent.mt"]] <- PercentageFeatureSet(
  BRBMET2_seurat,
  pattern = "^MT-"
)

VlnPlot(
  BRBMET2_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(
  BRBMET2_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
FeatureScatter(
  BRBMET2_seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

BRBMET2_filtered <- subset(
  BRBMET2_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

VlnPlot(
  BRBMET2_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(BRBMET2_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

BRBMET2_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

BRBMET2_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

table(BRBMET2_filtered$doublet_class)

library(ggplot2)
ggplot(BRBMET2_filtered@meta.data,
       aes(x = nCount_RNA,
           y = doublet_score,
           color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score"
  ) +
  theme_classic()

VlnPlot(
  BRBMET2_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)
BRBMET2_filtered <- subset(BRBMET2_filtered, subset = doublet_class == "singlet")
#===========================================
library(Seurat)
sc_data <- Read10X(data.dir = "D:/oza/GSE234832_RAW/BRBMET3/")
BRBMET3_seurat <- CreateSeuratObject(counts = sc_data, project = "BRBMET3")
BRBMET3_seurat[["percent.mt"]] <- PercentageFeatureSet(
  BRBMET3_seurat,
  pattern = "^MT-"
)

VlnPlot(
  BRBMET3_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(
  BRBMET3_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
FeatureScatter(
  BRBMET3_seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

BRBMET3_filtered <- subset(
  BRBMET3_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

VlnPlot(
  BRBMET3_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(BRBMET3_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

BRBMET3_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

BRBMET3_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

table(BRBMET3_filtered$doublet_class)

library(ggplot2)

ggplot(BRBMET3_filtered@meta.data,
       aes(x = nCount_RNA,
           y = doublet_score,
           color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score"
  ) +
  theme_classic()

VlnPlot(
  BRBMET3_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)
BRBMET3_filtered <- subset(BRBMET3_filtered, subset = doublet_class == "singlet")
#===========================================
library(Seurat)
sc_data <- Read10X(data.dir = "D:/oza/GSE234832_RAW/BRBMET87/")

BRBMET87_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET87"
)
BRBMET87_seurat[["percent.mt"]] <- PercentageFeatureSet(
  BRBMET87_seurat,
  pattern = "^MT-"
)

VlnPlot(
  BRBMET87_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(
  BRBMET87_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
FeatureScatter(
  BRBMET87_seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

BRBMET87_filtered <- subset(
  BRBMET87_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

VlnPlot(
  BRBMET87_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(BRBMET87_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

BRBMET87_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

BRBMET87_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

table(BRBMET87_filtered$doublet_class)

library(ggplot2)
ggplot(BRBMET87_filtered@meta.data,
       aes(x = nCount_RNA,
           y = doublet_score,
           color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score"
  ) +
  theme_classic()

VlnPlot(
  BRBMET87_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

BRBMET87_filtered <- subset(BRBMET87_filtered, subset = doublet_class == "singlet")
#===========================================
library(Seurat)
sc_data <- Read10X(data.dir = "D:/oza/GSE234832_RAW/LUBMET1/")

LUBMET1_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "LUBMET1"
)

LUBMET1_seurat[["percent.mt"]] <- PercentageFeatureSet(
  LUBMET1_seurat,
  pattern = "^MT-"
)

VlnPlot(
  LUBMET1_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(
  LUBMET1_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

FeatureScatter(
  LUBMET1_seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

LUBMET1_filtered <- subset(
  LUBMET1_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

VlnPlot(
  LUBMET1_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)
library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(LUBMET1_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

LUBMET1_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

LUBMET1_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

table(LUBMET1_filtered$doublet_class)

library(ggplot2)

ggplot(LUBMET1_filtered@meta.data,
       aes(x = nCount_RNA,
           y = doublet_score,
           color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score"
  ) +
  theme_classic()

VlnPlot(
  LUBMET1_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

LUBMET1_filtered <- subset(LUBMET1_filtered, subset = doublet_class == "singlet")
#===========================================
library(Seurat)
sc_data <- Read10X(data.dir = "D:/oza/GSE234832_RAW/LUBMET7/")

LUBMET7_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "LUBMET7"
)
LUBMET7_seurat[["percent.mt"]] <- PercentageFeatureSet(
  LUBMET7_seurat,
  pattern = "^MT-"
)

VlnPlot(
  LUBMET7_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

FeatureScatter(
  LUBMET7_seurat,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)
FeatureScatter(
  LUBMET7_seurat,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)

LUBMET7_filtered <- subset(
  LUBMET7_seurat,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

VlnPlot(
  LUBMET7_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(LUBMET7_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

LUBMET7_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

LUBMET7_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

table(LUBMET7_filtered$doublet_class)

library(ggplot2)

ggplot(LUBMET7_filtered@meta.data,
       aes(x = nCount_RNA,
           y = doublet_score,
           color = doublet_class)) +
  geom_point(size = 1, alpha = 0.6) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score"
  ) +
  theme_classic()

VlnPlot(
  LUBMET7_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

LUBMET7_filtered <- subset(LUBMET7_filtered, subset = doublet_class == "singlet")
#===========================================
merged_seurat <- merge(
  x = BRBMET2_filtered,
  y = c(BRBMET3_filtered, BRBMET87_filtered, LUBMET7_filtered, LUBMET1_filtered),
  add.cell.ids = c("BRBMET2", "BRBMET3", "BRBMET87", "LUBMET7", "LUBMET1"),
  project = "GSE234832"
)
merged_seurat

merged_seurat <- NormalizeData(
  merged_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

merged_seurat <- ScaleData(merged_seurat)

merged_seurat <- FindVariableFeatures(
  merged_seurat,
  selection.method = "vst",
  nfeatures = 2000
)

VariableFeaturePlot(merged_seurat)
head(VariableFeatures(merged_seurat), 10)
plot1 <- VariableFeaturePlot(merged_seurat)

top10 <- head(VariableFeatures(merged_seurat), 10)
LabelPoints(
  plot = plot1,
  points = top10,
  repel = TRUE
)

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
# ============================
#integration across samples using Harmony
# ============================
library(harmony)
merged_seurat <- RunHarmony(
  merged_seurat,
  group.by.vars = "orig.ident",
  dims.use = 1:30,
  theta = 4  )

# ============================
# Step: Clustering & UMAP
# (Finding cell groups + 2D visualization based on top 30 PCs)
# ============================
merged_seurat <- FindNeighbors(merged_seurat, reduction = "harmony", dims = 1:30)
merged_seurat <- FindClusters(merged_seurat, resolution = 0.5)

merged_seurat <- RunUMAP(merged_seurat, reduction = "harmony", dims = 1:30)
DimPlot(merged_seurat, reduction = "umap", label = TRUE)

DimPlot(merged_seurat, reduction = "umap", group.by = "orig.ident")

merged_seurat[["RNA"]] <- JoinLayers(merged_seurat[["RNA"]])

#====================================================

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
#============================================
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
#==================================================
