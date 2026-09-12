#===========================================
# Load and Create Seurat Object

library(Seurat)
sc_data <- Read10X(data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/BRBMET2/")

BRBMET2_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET2"
)

# Calculate Mitochondrial Percentage

BRBMET2_seurat[["percent.mt"]] <- PercentageFeatureSet(
  BRBMET2_seurat,
  pattern = "^MT-"
)

# Visualize QC metrics

VlnPlot(
  BRBMET2_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Generate QC scatter plots

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

# Filter low-quality cells

BRBMET2_filtered <- subset(
  BRBMET2_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

# Visualize QC metrics after filtering

VlnPlot(
  BRBMET2_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Detect doublets

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(BRBMET2_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

# Add doublet scores and classifications

BRBMET2_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

BRBMET2_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

# Summarize doublet classifications

table(BRBMET2_filtered$doublet_class)

# Visualize doublet scores

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

# Visualize doublet score by classification

VlnPlot(
  BRBMET2_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

# Remove doublets

BRBMET2_filtered <- subset(BRBMET2_filtered, subset = doublet_class == "singlet")
#===========================================
# Load and Create Seurat Object

library(Seurat)
sc_data <- Read10X(data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/BRBMET3/")

BRBMET3_seurat <- CreateSeuratObject(counts = sc_data, project = "BRBMET3")

# Calculate Mitochondrial Percentage

BRBMET3_seurat[["percent.mt"]] <- PercentageFeatureSet(
  BRBMET3_seurat,
  pattern = "^MT-"
)

# Visualize QC metrics

VlnPlot(
  BRBMET3_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Generate QC scatter plots

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

# Filter low-quality cells

BRBMET3_filtered <- subset(
  BRBMET3_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

# Visualize QC metrics after filtering

VlnPlot(
  BRBMET3_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Detect doublets

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(BRBMET3_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

# Add doublet scores and classifications

BRBMET3_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

BRBMET3_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

# Summarize doublet classifications

table(BRBMET3_filtered$doublet_class)

# Visualize doublet scores

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

# Visualize doublet score by classification

VlnPlot(
  BRBMET3_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

# Remove doublets

BRBMET3_filtered <- subset(BRBMET3_filtered, subset = doublet_class == "singlet")
#===========================================
# Load and Create Seurat Object

library(Seurat)
sc_data <- Read10X(data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/BRBMET87/")

BRBMET87_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "BRBMET87"
)

# Calculate mitochondrial percentage

BRBMET87_seurat[["percent.mt"]] <- PercentageFeatureSet(
  BRBMET87_seurat,
  pattern = "^MT-"
)

# Visualize QC metrics

VlnPlot(
  BRBMET87_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Generate QC scatter plots

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

# Filter low-quality cells

BRBMET87_filtered <- subset(
  BRBMET87_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

# Visualize QC metrics after filtering

VlnPlot(
  BRBMET87_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Detect doublets

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(BRBMET87_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

# Add doublet scores and classifications

BRBMET87_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

BRBMET87_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

# Summarize doublet classifications

table(BRBMET87_filtered$doublet_class)

# Visualize doublet scores

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

# Visualize doublet score by classification

VlnPlot(
  BRBMET87_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

# Remove doublets

BRBMET87_filtered <- subset(BRBMET87_filtered, subset = doublet_class == "singlet")
#===========================================
# Load and Create Seurat Object

library(Seurat)
sc_data <- Read10X(data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/LUBMET1/")

LUBMET1_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "LUBMET1"
)

# Calculate mitochondrial percentage

LUBMET1_seurat[["percent.mt"]] <- PercentageFeatureSet(
  LUBMET1_seurat,
  pattern = "^MT-"
)

# Visualize QC metrics

VlnPlot(
  LUBMET1_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Generate QC scatter plots

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

# Filter low-quality cells

LUBMET1_filtered <- subset(
  LUBMET1_seurat,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

# Visualize QC metrics after filtering

VlnPlot(
  LUBMET1_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Detect doublets

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(LUBMET1_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

# Add doublet scores and classifications

LUBMET1_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

LUBMET1_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

# Summarize doublet classifications

table(LUBMET1_filtered$doublet_class)

# Visualize doublet scores

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

# Visualize doublet score by classification

VlnPlot(
  LUBMET1_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

# Remove doublets

LUBMET1_filtered <- subset(LUBMET1_filtered, subset = doublet_class == "singlet")
#===========================================
# Load and Create Seurat Object

library(Seurat)
sc_data <- Read10X(data.dir = "/Users/Technology - Laptoop/Downloads/GSE234832_RAW/LUBMET7/")

LUBMET7_seurat <- CreateSeuratObject(
  counts = sc_data,
  project = "LUBMET7"
)

# Calculate mitochondrial percentage

LUBMET7_seurat[["percent.mt"]] <- PercentageFeatureSet(
  LUBMET7_seurat,
  pattern = "^MT-"
)

# Visualize QC metrics

VlnPlot(
  LUBMET7_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Generate QC scatter plots

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

# Filter low-quality cells

LUBMET7_filtered <- subset(
  LUBMET7_seurat,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 40
)

# Visualize QC metrics after filtering

VlnPlot(
  LUBMET7_filtered,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Detect doublets

library(scDblFinder)
library(SingleCellExperiment)
sce_d1 <- as.SingleCellExperiment(LUBMET7_filtered)
set.seed(100)
sce_d1 <- scDblFinder(sce_d1)

# Add doublet scores and classifications

LUBMET7_filtered$doublet_score <-
  colData(sce_d1)$scDblFinder.score

LUBMET7_filtered$doublet_class <-
  colData(sce_d1)$scDblFinder.class

# Summarize doublet classifications

table(LUBMET7_filtered$doublet_class)

# Visualize doublet scores

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

# Visualize doublet score by classification

VlnPlot(
  LUBMET7_filtered,
  features = "doublet_score",
  group.by = "doublet_class"
)

# Remove doublets

LUBMET7_filtered <- subset(LUBMET7_filtered, subset = doublet_class == "singlet")
#===========================================
# Merge all samples

merged_seurat <- merge(
  x = BRBMET2_filtered,
  y = c(BRBMET3_filtered, BRBMET87_filtered, LUBMET7_filtered, LUBMET1_filtered),
  add.cell.ids = c("BRBMET2", "BRBMET3", "BRBMET87", "LUBMET7", "LUBMET1"),
  project = "GSE234832"
)
merged_seurat

# Normalize the data

merged_seurat <- NormalizeData(
  merged_seurat,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Scale the data

merged_seurat <- ScaleData(merged_seurat)

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