# Title: Uncovering Intercellular Communication Pathways Between Fibroblast and Myeloid Subclusters

## Project Description
Brain metastases represent a major clinical challenge characterized by a highly specialized tumor microenvironment (TME). To dissect the cellular architecture and intercellular cross-talk within this niche, this project presents a single-cell RNA sequencing (scRNA-seq) **re-analysis** of the human brain metastatic microenvironment. By integrating high-resolution single-cell transcriptomics with ligand-receptor network modeling, we systematically map the heterogeneous cell states and signaling crosstalk between specific tumor-associated fibroblast (TAF) cluster and M1 macrophage populations. Unraveling these specialized intercellular communication pathways provides critical insights into TME remodeling, offering prospective targets for therapeutic intervention and microenvironment-directed strategies.

### Background 


### Rationale


### Input Data
GEO? 
Single-cell RNA sequencing data from five brain metastasis samples were loaded and processed using Seurat. The dataset included three breast cancer brain metastasis samples and two lung cancer brain metastasis samples.
The samples were identified as:
•	BRBMET2
•	BRBMET3
•	BRBMET87
•	LUBMET7
•	LUBMET1
The expression matrices, feature annotations, and cell barcodes were loaded for each sample and converted into Seurat objects for downstream analysis.

### Methods

#### Data Loading
The five samples (BRBMET2, BRBMET3, BRBMET87, LUBMET7, and LUBMET1) were read from sparse matrix format (matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz) using the ReadMtx function in Seurat. Individual Seurat objects were constructed for each sample without initial threshold filtering (min.cells = 0, min.features = 0).

#### Quality Control and Doublet Removal
Cell quality control was performed independently for each sample using adaptive thresholding based on median absolute deviations (MAD). Mitochondrial gene expression percentages (percent.mt) were calculated using PercentageFeatureSet for genes matching the ^MT- prefix. Cells were retained if their mitochondrial transcript percentage was below an upper boundary defined as min({median} + 3*{MAD}, 40\%) and their detected gene counts (nFeature_RNA) fell within max(200, {median} - 3*{MAD}) and {median} + 3*{MAD}. High-quality cells across all five filtered samples were merged into a single Seurat object and joined across layers. To remove heterotypic doublets, the merged dataset was converted into a SingleCellExperiment object, and scDblFinder was executed with sample-level batch stratification (samples = "orig.ident"). Cells classified as doublets were removed, retaining only verified singlets for downstream processing.

#### Normalization, Identification of Variable Features, and Scaling
Gene expression counts of the singlet population were normalized using standard log-normalization via NormalizeData, scaling feature counts per cell by a factor of 10,000 and natural log-transforming the result (LogNormalize). Highly variable features were identified using the variance-stabilizing transformation method (selection.method = "vst") in FindVariableFeatures, selecting the top 2,000 variable genes. Linear scaling was subsequently applied across all genes using ScaleData to center expression measurements to a mean of zero and scale variance to one.

#### Principal Component Analysis (PCA)
Linear dimensional reduction was performed on the scaled expression matrix of the top 2,000 highly variable features using RunPCA. Dimensional loadings and primary sources of variation were evaluated across the top principal components. The proportion of variance explained by individual PCs was assessed using elbow plot visualization (ElbowPlot) over 50 dimensions to determine the optimal component cutoff for batch alignment and downstream embeddings.

#### Harmony Integration
To mitigate sample-specific batch effects across the five metastasis specimens, batch correction was executed on the top 30 principal components using the Harmony algorithm (RunHarmony) with sample identity (orig.ident) set as the batch variable and a diversity penalty parameter (theta = 4). Harmony embeddings (harmony) were generated for downstream dimensional reduction and graph-based clustering.

#### Clustering
Graph-based unsupervised clustering was performed before and after Harmony batch correction for comparison. Cell-cell similarity graphs were constructed using FindNeighbors based on either the first 30 uncorrected principal components or the first 30 Harmony-integrated dimensions. Community detection was subsequently executed using FindClusters at a modularity resolution of 0.5.

#### Clusters Visualization (UMAP and t-SNE)
For non-linear dimensional reduction and spatial visualization, Uniform Manifold Approximation and Projection (UMAP) embeddings were generated before (umap_pca) and after batch correction (umap) using the first 30 dimensions. In addition, t-Distributed Stochastic Neighbor Embedding (t-SNE) was performed on the 30 Harmony dimensions (RunTSNE, reduction.name = "tsne_harmony"). Cluster distributions and batch mixing across individual samples were evaluated on UMAP and t-SNE projections.

#### Identification of Clusters Marker Genes
Differentially expressed gene (DEG) markers defining each integrated cluster were identified using Seurat's FindAllMarkers function. Non-parametric Wilcoxon rank-sum tests were conducted evaluating positively expressed genes (only.pos = TRUE) detected in a minimum fraction of 25% of cells in either comparison group (min.pct = 0.25) with a natural log fold-change threshold of at least 0.25 (logfc.threshold = 0.25). The top 5 marker genes per cluster were ranked by average log2 fold-change (avg_log2FC).

#### Identification and Sub-clustering of Myeloid Cells
The myeloid cell population was identified by evaluating average cluster-level expression of canonical myeloid markers (TREM2, C1QB, GPR34, FOLR2, and C1QC) using AverageExpression. Cells belonging to the primary myeloid cluster (Cluster 0) were isolated and re-processed. Highly variable features (top 2,000 genes) were re-selected, scaled, and subjected to PCA. Harmony batch integration was re-applied across the top 20 principal components (dims.use = 1:20, theta = 4). Sub-clustering was performed at a resolution of 0.5 on the Harmony neighborhood graph, followed by UMAP projection (RunUMAP). Sub-cluster-specific marker genes were identified via FindAllMarkers, and canonical lineages/states were validated using feature expression plots (FeaturePlot) for myeloid lineage.

#### Automated Annotation of the defined clusters
Automated cell-type annotation of the single-cell transcriptomic dataset was performed using the SingleR package. The Seurat object was converted into a SingleCellExperiment object and queried against two reference datasets from the celldex package: the Human Primary Cell Atlas Data (HumanPrimaryCellAtlasData) and the Monaco Immune Data (MonacoImmuneData). Main cell-type labels (label.main) from both references were utilized simultaneously to assign predictions to individual cells. The resulting predicted labels were added back to the Seurat object metadata and visualized using Uniform Manifold Approximation and Projection (UMAP) dimensionality reduction plots generated via DimPlot.

#### Differential expression annotation of the defined clusters
To biologically annotate and refine the unsupervised clusters, marker genes for each cluster were identified using Seurat's FindAllMarkers function, restricting evaluation to upregulated genes (only.pos = TRUE) with a minimum expression fraction threshold (min.pct = 0.25) and a log2 fold-change threshold (logfc.threshold = 0.25). The top 25 marker genes per cluster based on average log2 fold-change (avg_log2FC) were extracted and saved. Clusters 0 through 14 were manually annotated into broad lineages and tumor states based on marker expression signatures.  Identified cell types were assigned to active cluster identities and visualized on UMAP plots.

#### A heatmap as QC to Check for NAs
Prior to expression visualization, quality control checks were conducted on the cell type metadata (cell_type) to verify data integrity and confirm the absence of missing values (NA). Cells with valid annotations were retained, and expression patterns of selected marker genes (genes_to_plot) across cell-type identities were visualized using Seurat's DoHeatmap. Gene expression values were scaled and displayed using a custom color gradient.

#### Differential expression annotation of the Myeloid sub-clusters
Sub-clustering analysis was focused on the myeloid lineage. Differential expression analysis was performed on the myeloid sub-clusters, and significant markers (padj < 0.05) were filtered. The top 25 markers per sub-cluster ranked by avg_log2FC were exported (top25_myeloid_markers.csv), followed by UMAP visualization.

#### Determining of M1 and M2 states
To characterize functional polarization states within the myeloid population, module scores for canonical M1 (pro-inflammatory) and M2 (anti-inflammatory/immunosuppressive) gene sets were calculated using Seurat's AddModuleScore function. Module scores were evaluated across myeloid sub-clusters via violin plots (VlnPlot), and mean scores as well as score differentials (M1 - M2) were calculated for each sub-cluster.

#### Annotation of M1 and M2 states
Myeloid sub-clusters were consolidated into broader lineages, merging cluster 0 and 1 into "Macrophage" while maintaining "Microglia" and other non-macrophage/microglia populations. A net polarization score was calculated for each cell as Delta S = {M1_score1} - {M2_score1}. 
Polarization states were categorized based on a threshold (= 0.1): cells with Delta S > 0.1 were assigned as "M1-like", cells with Delta S < -0.1 as "M2-like", and cells falling between (-0.1 and 0.1) as "Intermediate". For Macrophage and Microglia lineages, composite annotations combining cell identity and polarization state (e.g., "Macrophage_M1-like") were assigned and visualized using UMAP projections.

#### Visualizations of M1 and M2 states
Polarization dynamics and cell-state distributions were visualized through multiple complementary approaches: (1) UMAP projections of cell types split by polarization state (M1-like, M2-like, Intermediate); (2) UMAP feature plots displaying continuous Delta S values using a divergent color gradient (blue for M1-leaning, red for M2-leaning); (3) stacked bar charts (ggplot2) displaying the relative percentage distribution of M1-like, M2-like, and Intermediate states strictly within Macrophage and Microglia populations; and (4) sub-cluster UMAP feature maps paired with top marker heatmaps (top 5 markers per cluster, padj < 0.05).


### Results

#### Data Loading
Single-cell transcriptomic profiles were successfully loaded across all five brain metastasis specimens (BRBMET2, BRBMET3, BRBMET87, LUBMET7, and LUBMET1), generating independent initial Seurat objects for each dataset.

#### Quality Control and Doublet Removal
Sample-specific adaptive thresholding based on median absolute deviations effectively filtered out low-quality dying cells with high mitochondrial read content and unviable transcript counts. In silico doublet detection using scDblFinder identified and removed artificial cell multiplets from the merged dataset, yielding a clean singlet population for downstream integrative analyses.

#### Normalization, Identification of Variable Features, and Scaling
Log-normalization and variance-stabilizing transformation successfully identified 2,000 hypervariable genes across the singlet dataset. Linear scaling centered gene expression values appropriately across features, setting up the dataset for linear dimensional reduction.

#### Principal Component Analysis (PCA)
Principal component analysis effectively captured major biological axes of variation across the top dimensions. Inspection of the elbow plot confirmed that the top 30 principal components accounted for the majority of transcriptomic variance within the dataset.

#### Harmony Integration
Harmony integration effectively harmonized expression profiles across the five distinct metastasis samples. Inter-sample technical variation was successfully removed, placing equivalent biological cell states from different patients into shared coordinate spaces while preserving biological heterogeneity.

#### Clustering
Graph-based community detection at a resolution of 0.5 generated distinct cell clusters across the dataset. Comparative analysis demonstrated that clustering on Harmony-integrated dimensions yielded coherent cell groupings free from sample-specific batch artifacts compared to clustering on uncorrected PCA space.

#### Clusters Visualization (UMAP and t-SNE)
Both UMAP and t-SNE projections based on Harmony dimensions clearly resolved discrete cellular populations. Color-coding embeddings by sample identity (orig.ident) confirmed thorough inter-sample mixing across all shared clusters, whereas cluster-labeled projections demonstrated well-separated biological lineages.

#### Identification of Clusters Marker Genes
Differential expression testing with FindAllMarkers identified robust set-specific upregulated marker profiles for each cluster (padj < 0.05, log_2FC > 0.25). Extraction of top marker genes ranked by average log2 fold-change established clear cell-identity signatures defining each distinct parenchymal, immune, and stromal cluster.

#### Identification and Sub-clustering of Myeloid Cells
Expression profiling of canonical myeloid genes (TREM2, C1QB, GPR34, FOLR2, C1QC) unambiguously identified Cluster 0 as the primary myeloid lineage compartment. Re-clustering and sub-analysis of this isolated population following Harmony re-integration revealed distinct myeloid sub-clusters on UMAP projections. Differential marker analysis and feature plots confirmed distinct myeloid functional states—including C1QC/ C1QB/ TREM2/ FOLR2/ LYZ-positive macrophage and microglial subsets as well as MKI67-positive proliferating myeloid cells—while demonstrating complete absence of lymphoid lineage contaminants (CD3D, CD3E, NKG7, GNLY).


#### Automated Annotation of the defined clusters
Automated annotation using SingleR against the Human Primary Cell Atlas and Monaco Immune reference datasets successfully assigned broad reference labels across the dataset, providing an initial  overview of cell identities prior to marker-based cluster refinement.

#### Differential expression annotation of the defined clusters
Differential expression analysis via FindAllMarkers identified robust gene expression profiles distinguishing all 15 major clusters (clusters 0–14). 

#### Create a heatmap and Check for NAs first
Quality control evaluation confirmed complete metadata integrity across all annotated cells. Expression heatmap visualization (DoHeatmap) demonstrated distinct block-like marker gene expression patterns corresponding strictly to assigned cluster identities, confirming high specificity of marker panels and minimal cross-cluster contamination.

#### Differential expression annotation of the Myeloid sub-clusters
Focused sub-clustering of the myeloid lineage resolved five distinct sub-populations with distinct transcriptomic signatures. These were successfully categorized into classical Macrophages (sub-cluster 0), a stress/metabolically active MMP9-high mitochondrial Macrophage state (sub-cluster 1), Microglia (sub-cluster 2), conventional type 2 Dendritic Cells (cDC2, sub-cluster 3), and Proliferating Myeloid cells (sub-cluster 4).

#### Determining of M1 and M2 states
Gene module scoring revealed distinct transcriptomic continuum states of M1 pro-inflammatory and M2 anti-inflammatory activation across the myeloid subtypes. Comparison of mean score differences ({M1} - {M2}) via violin plots highlighted heterogeneous polarization preferences across sub-clusters, demonstrating that classical discrete M1/M2 dichotomies are better represented as continuous signature distributions in single-cell resolution.
#### Annotation of M1 and M2 states
Categorization based on net polarization score differentials (Delta S = {M1} - {M2}) with a 0.1 threshold successfully segregated macrophages and microglia into distinct functional subpopulations ("M1-like", "M2-like", and "Intermediate"). UMAP visualization of composite state annotations highlighted clear spatial clustering of M1-like versus M2-like polarization phenotypes within the tumor microenvironment.
#### Visualizations of M1 and M2 states
Multi-modal visualization confirmed spatial and proportional segregation of myeloid states:
- Split UMAP Projections: Demonstrated distinct spatial density differences between M1-like, M2-like, and Intermediate cell subsets across Macrophage and Microglia lineages.


- Polarization Gradient: Feature plots of continuous score differentials (Delta S) displayed a clear spectrum ranging from strongly M1-leaning (blue) to strongly M2-leaning (red) domains.


- Proportional Quantification: Proportion analysis demonstrated the relative breakdown of M1-like, M2-like, and Intermediate phenotypes within Macrophages versus Microglia, quantifying macrophage/microglial polarization dynamics within the microenvironment.


- Subtype Expression Heatmap: Marker gene visualization confirmed top lineage markers defining each sub-cluster state without expression overlap.


### Discussion

### Conclusion
