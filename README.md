# Title: Uncovering Intercellular Communication Pathways Between Fibroblast and Myeloid Subclusters

## Project Description
Brain metastases represent a major clinical challenge characterized by a highly specialized tumor microenvironment (TME). To dissect the cellular architecture and intercellular cross-talk within this niche, this project presents a single-cell RNA sequencing (scRNA-seq) **re-analysis** of the human brain metastatic microenvironment. By integrating high-resolution single-cell transcriptomics with ligand-receptor network modeling, we systematically map the heterogeneous cell states and signaling crosstalk between specific tumor-associated fibroblast (TAF) cluster and M1 macrophage populations. Unraveling these specialized intercellular communication pathways provides critical insights into TME remodeling, offering prospective targets for therapeutic intervention and microenvironment-directed strategies.

### Background 
Brain metastases affect approximately 200,000 patients annually in the US, occurring ten times more frequently than primary brain tumors and carrying high mortality and recurrence rates. While advances in targeted therapies and radiosurgery have improved patient management, effective biomarkers and therapeutic targets remain urgently needed. Evidence indicates that the tumor microenvironment (TME) significantly dictates response to radiotherapy and overall survival, with pan-brain metastasis markers, such as type I collagen genes (COL1A1/COL1A2), playing key roles across diverse primary cancer origins. Leveraging single-cell RNA sequencing (scRNA-seq) on five brain metastasis samples, this study constructs a high-resolution cell landscape, identifying type I collagen-high tumor-associated fibroblasts (TAFs) as central mediators of TME intercellular communication. The findings suggest that $M1$ activation of resident microglia and infiltrating macrophages establishes a pro-inflammatory, wound-healing response that metastatic tumor cells hijack—driving excess collagen expression to remodel the brain niche into a hospitable microenvironment and offering novel predictive biomarkers and therapeutic targets for metastatic brain cancer and glioblastoma.

### Rationale
We hypothesized that tumor-associated fibroblasts (TAFs) drives immunosuppressive and structural remodeling of the brain metastatic microenvironment by engaging in distinct ligand-receptor signaling crosstalk with M1 macrophage populations. We reasoned that high-resolution single-cell re-analysis and intercellular network modeling would resolve these specific communication axes, demonstrating that targeting this TAF–macrophage signaling circuit can disrupt pro-tumoral TME remodeling and reveal novel therapeutic vulnerabilities in human brain metastases.

### Input Data
Single-cell RNA-seq data is available in the NCBI Gene expression Omnibus database (GEO) with accession number GSE234832.
Data from five brain metastasis samples were loaded and processed using Seurat. The dataset included three breast cancer brain metastasis samples and two lung cancer brain metastasis samples.
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

#### Differential Expression Analysis and Volcano Plot Visualization
Differential expression analysis was performed using a two-sided Wilcoxon rank-sum test via the FindMarkers function in Seurat with a minimum cell fraction detection threshold (min.pct) of 0.1 and an absolute log2 fold-change threshold of 0. Differential expression results were categorized into significance tiers based on adjusted p-values calculated using the Benjamini-Hochberg false discovery rate (FDR) correction: highly significant (FDR < 0.01 and log_2FC > 1), significant (FDR < 0.05), or non-significant. Volcano plots were constructed using ggplot2 and ggrepel. To prevent visual distortion from extreme pvalues, -log10pvalues were capped at 300. For representative gene annotation, the top 12 upregulated and downregulated genes per contrast were selected based on adjusted pvalue and log2 fold-change filtering, excluding lowly expressed outliers (cell fraction < 0.1) and extreme fold-change fluctuations (log_2FC > 6).

#### Functional Enrichment Analysis (GO, KEGG, Reactome)
Directional functional enrichment analysis was conducted separately for significantly upregulated (log_2FC > 1, FDR < 0.05) and downregulated (log_2FC < -1, FDR < 0.05) genes. Gene symbols were mapped to Entrez Gene IDs using the org.Hs.eg.db R package via clusterProfiler::bitr. 
Over-representation analysis (ORA) was executed across Gene Ontology (GO) categories—Biological Process (BP), Molecular Function (MF), and Cellular Component (CC), as well as Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways using clusterProfiler, and Reactome pathways using ReactomePA. Multiple-testing adjustments were performed using the Benjamini-Hochberg method, with terms exhibiting an adjusted pvalue < 0.05 considered statistically enriched.

#### Pathway Enrichment Visualization
To facilitate intuitive comparison of biological themes across expression states, enriched functional terms were visualized using bidirectional diverging bar charts in ggplot2. The top 10 most significantly enriched terms (FDR < 0.05) for upregulated and downregulated gene sets were plotted along a continuous normalized axis. Enriched pathways associated with upregulated genes were represented along positive values, whereas downregulated pathways extended along negative values. To handle potential duplicate descriptions across categories, uniquely formatted labels were mapped back to their original term names while preserving factor ordering.

#### Differential Expression and Pathway Enrichment of Fibroblasts
Fibroblasts (pericyte/fibroblast cluster 11) were evaluated for differential gene expression against all remaining cell populations in the single-cell dataset using the Wilcoxon rank-sum test. Identified marker genes (FDR < 0.01, log_2FC > 1) were saved, and directional enrichment analysis was performed across REACTOME, KEGG, and GO ontologies (BP, MF, CC). Separate diverging bar plots were generated to visualize enriched biological functions characteristic of the fibroblast population.

#### Differential Expression of M1-like vs. M2-like Macrophages
Myeloid single-cell clusters were annotated and subsetted to isolate macrophages and microglia. Polarization module scores for M1-like and M2-like transcriptional signatures were computed per cell using AddModuleScore. Individual cells were categorized into polarization states based on the score differential (M1 - M2): M1-like (> 0.1), M2-like (< -0.1), or Intermediate. A direct differential expression contrast was then executed between M1-like macrophages (ident.1) and M2-like macrophages (ident.2) using the volcano plot.

#### Pathway Enrichment Analysis of Macrophage Polarization States
Genes differentially expressed between M1-like and M2-like macrophages (FDR < 0.05, log_2FC > 1) were subjected to ORA using clusterProfiler and ReactomePA. Biological enrichment patterns were evaluated across REACTOME pathways, KEGG pathways, and GO sub-ontologies. Enriched processes specific to M1-like (positive scores) versus M2-like (negative scores) states were displayed as diverging bar plots.

#### CellChat Model Construction and Ligand-Receptor Database Mapping
Single-cell transcriptomic profiles were integrated across datasets to focus on microenvironmental interaction dynamics. Myeloid sub-clusters annotated by functional polarization state replaced the myeloid cluster, followed by assay layer integration (JoinLayers) under the default RNA assay within Seurat. The dataset was subsequently filtered to isolate five key cell populations of interest: pericytes/fibroblasts (Pericyte_fibroblast), M1-like macrophages (Macrophage_M1-like), M2-like macrophages (Macrophage_M2-like), M1-like microglia (Microglia_M1-like), and M2-like microglia (Microglia_M2-like).

Inferred cell-cell communication networks were reconstructed using the CellChat R package. Expression matrices and cell metadata annotations were extracted to initialize a CellChat object grouped by DEG_annotation. Human ligand-receptor interactions were assigned from CellChatDB.human. To optimize computational performance without altering signal detection, the expression matrix was subsetted (subsetData) to retain only genes corresponding to established ligands and receptors present within the curated database.

#### Identification of Overexpressed Signal Pairs and Network Inference
Overexpressed signaling genes across the defined cell clusters were identified using identifyOverExpressedGenes paired with presto-accelerated differential expression. Overexpressed ligand-receptor interactions were subsequently determined (identifyOverExpressedInteractions) by filtering for pairs in which both signaling partners demonstrated significant group-wise overexpression. Communication probabilities between cell clusters were calculated (computeCommunProb) using law of mass action modeling based on average ligand and receptor expression levels. Interactions originating from rare cell groups containing fewer than 10 cells were removed using filterCommunication to mitigate low-confidence observational noise.

#### Pathway Communication Network Aggregation
Individual ligand-receptor communication probabilities were aggregated to evaluate signaling strengths at the pathway level (computeCommunProbPathway). Total network interactions were synthesized using aggregateNet to generate cell-group interaction matrices corresponding to raw signaling counts and weighted communication probabilities (strength). In addition, cell-type-specific signaling subsets were generated by filtering communication outputs specifically for pericyte/fibroblast source populations (sources.use = "Pericyte_fibroblast").

#### Network Visualization and Ligand-Receptor Interactions
Macro-level communication structure was visualized using circular network plots (netVisual_circle) displaying total interaction counts and weighted interaction strengths across all cell groups. Enriched signaling pathways were extracted, and individual pathway networks were rendered into high-resolution circular plots using automated layout functions (netVisual_aggregate). Communication data frames containing inferred probabilities and significance metrics (p < 0.01) were extracted via subsetCommunication. Custom processing pipelines generated faceted dot plots visualizing communication probabilities across target cell pairs, evaluated both globally and specifically for fibroblast-directed interactions, while systematically comparing networks with and without self-communication.

### Results

#### Data Loading
Single-cell transcriptomic profiles were successfully loaded across all five brain metastasis specimens (BRBMET2, BRBMET3, BRBMET87, LUBMET7, and LUBMET1), generating independent initial Seurat objects for each dataset.

#### Quality Control and Doublet Removal
Sample-specific adaptive thresholding based on median absolute deviations effectively filtered out low-quality dying cells with high mitochondrial read content and unviable transcript counts. In silico doublet detection using scDblFinder identified and removed artificial cell multiplets from the merged dataset, yielding a clean singlet population for downstream integrative analyses.

#### Normalization, Identification of Variable Features, and Scaling
After normalisation, identification of highly variable features across the single-cell dataset yielded 2,000 highly variable genes out of 33,538 total detected features, with ACY3 and COL3A1 exhibiting the highest standardized variance.

![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/02_Labeled%20VariableFeaturePlot.png)

Figure x. Variable feature selection and principal component analysis. Variance plot showing 2,000 highly variable genes (red) selected from 33,538 total detected genes, with top variable genes labeled (ACY3, COL3A1).
#### Principal Component Analysis (PCA)
Principal component analysis effectively captured major biological axes of variation across the top dimensions. Principal component analysis (PCA) based on these variable genes revealed distinct sample-level driver patterns along the major axes of variation. PC_1 clearly separated LUBMET1 samples along the negative axis from BRBMET87 samples aligned along positive PC_2, while BRBMET2, BRBMET3, and LUBMET7 populated intermediate distribution spaces. Top loading genes driving PC_1 highlighted opposing transcriptional profiles, with oligodendrocyte- and CNS-associated markers (PLP1, TF, MOG, APLP1, CARNS1, NKX6-2, STMN4, ENPP2) enriched on one spectrum, and epithelial/metastatic and inflammatory drivers (SPDEF, CREB3L4, IER3, LY6E, HMGB3, S100A14, FXYD3, S100A10) defining the opposing axis.

![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/05_PCA_Top_Genes_Heatmap.png)

Figure x. PCA scatter plot showing single cells projected onto PC_1 and PC_2, colored by sample origin.   (B) Heatmap displaying the top driving genes along PC_1, separating CNS/oligodendrocyte markers (PLP1, MOG) from metastatic/inflammatory markers (SPDEF, S100A14).   
Inspection of the elbow plot confirmed that the top 30 principal components accounted for the majority of transcriptomic variance within the dataset.
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/06_PCA_Elbow_Plot.png)

Figure x. Heatmap displaying the top driving genes along PC_1, separating CNS/oligodendrocyte markers (PLP1, MOG) from metastatic/inflammatory markers (SPDEF, S100A14).   

#### Harmony Integration, Clustering, and Visualization (UMAP and t-SNE)
Unsupervised clustering across the dataset identified 16 distinct transcriptomic clusters (Clusters 0–15). Prior to integration, uncorrected UMAP embedding exhibited pronounced batch-driven separation, with individual clusters composed almost exclusively of single samples (e.g., BRBMET2, BRBMET87, and LUBMET7 segregating into isolated islands). Application of Harmony batch correction successfully mitigated sample-specific variation, driving multi-sample alignment across central shared clusters while preserving true biological heterogeneity. Both post-Harmony UMAP and t-SNE projections demonstrated consistent spatial topology, effectively grouping homologous cell states across diverse patient origins.

![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/07_UMAP_Plot_by_Sample.png)

Figure x. Evaluation of batch correction and single-cell embedding. Uncorrected UMAP plot colored by sample identity, showing pre-integration batch effects.
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/08_UMAP_Clusters_Plot.png)

Figure x. Evaluation of batch correction and single-cell embedding. Post-Harmony UMAP plot colored by unsupervised Seurat clusters (0–15).
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/09_UMAP_Plot_by_Sample.png)

Figure x. Evaluation of batch correction and single-cell embedding. Post-Harmony UMAP plot colored by sample identity, demonstrating integration across shared cluster spaces.
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/10_tSNE_Clusters_Plot.png)

Figure x. Evaluation of batch correction and single-cell embedding. Post-Harmony t-SNE plot colored by Seurat clusters.
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/11_tSNE_Plot_by_Sample.png)

Figure x. Evaluation of batch correction and single-cell embedding. Post-Harmony t-SNE plot colored by sample identity.

#### Identification and Sub-clustering of Myeloid Cells
To determine the principal components (PCs) capturing the highest proportion of biological variance, elbow plot evaluation revealed a distinct inflection point around PC 8–10, beyond which standard deviation plateaued. Non-linear dimensionality reduction using UMAP on the top principal components resolved 7 distinct transcriptomic sub-clusters (Clusters 0–6). Post-integration overlay confirmed uniform distribution across all five patient samples (BRBMET2, BRBMET3, BRBMET87, LUBMET1, and LUBMET7), demonstrating effective batch correction without sample-specific subset bias. Expression profiling of canonical lineage markers revealed that the subset is heavily enriched for myeloid and macrophage populations. High expression of complement genes (C1QC, C1QB) alongside lipid metabolism and phagocytic markers (TREM2, FOLR2, LYZ) was localized across the primary contiguous clusters (Clusters 0, 1, and 3–5). Lymphocyte and cytotoxic markers (CD3D, CD3E, NKG7, GNLY) showed minimal or highly localized expression, while proliferation markers (MKI67) defined a distinct, minor cycling population (Cluster 6)

![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/12_PCA_Elbow_Plot.png)

Figure x. Dimensionality reduction and sub-clustering of myeloid populations. PCA elbow plot illustrating the standard deviation accounted for by each principal component, with an elbow point observed near PC 8–10.
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/13_UMAP_Clusters_Plot.png)

Figure x. Dimensionality reduction and sub-clustering of myeloid populations. Post-Harmony UMAP embedding displaying 7 distinct sub-clusters (Clusters 0–6).
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/14_UMAP_Plot_by_Sample.png)

Figure x. Dimensionality reduction and sub-clustering of myeloid populations. UMAP projection colored by sample identity (orig.ident), demonstrating homogenous sample integration and lack of batch effect across sub-clusters
![Aggregated Interactions Count](figures/01_Preprocessing_and_Clustering/15_UMAP_Marker_Genes_FeaturePlot.png)

Figure x. Expression landscape of canonical lineage and functional marker genes. Feature plots displaying normalized expression levels of key cell-type markers projected on UMAP space. Myeloid and macrophage-associated markers (C1QC, C1QB, TREM2, FOLR2, LYZ) show broad expression across main clusters. T-cell (CD3D, CD3E), NK cell/cytotoxic (NKG7, GNLY), and proliferative (MKI67) markers show sparse or cluster-restricted expression.

#### Differential expression annotation of the defined clusters
Differential expression analysis via FindAllMarkers identified robust gene expression profiles distinguishing all 15 major clusters (clusters 0–14). Cell identity assignment was confirmed via distinct lineage-specific marker expression on the expression heatmap, including FOLR2, TREM2, and C1QC in myeloid subsets, CD3D, TRBC1, and NKG7 in $T/\text{NK}$ cells, and COL3A1 in pericyte/fibroblasts.  

![Aggregated Interactions Count](figures/02_Annotation/Cell_type_annotation.png)

Figure x. Global cellular landscape and marker gene expression in brain metastases. UMAP visualization of integrated single cells colored by major cell type annotations, encompassing malignant epithelial subsets, stromal populations, neuro-glial lineage, and immune compartments.   

#### Differential expression annotation of the Myeloid sub-clusters
Sub-clustering restricted to the myeloid compartment further dissected four distinct sub-populations: Macrophages, Macrophage_MMP9_high_mito, Microglia, and type 2 conventional dendritic cells (cDC2).
![Aggregated Interactions Count](figures/02_Annotation/Myeloid_subtypes_and_states.png)

Figure x. Myeloid subtype characterization and M1/M2 polarization states.UMAP embedding of high-resolution sub-clustered myeloid populations identifying Macrophages, Macrophage_MMP9_high_mito, Microglia, and cDC2s.

#### Determining, Annotation, and Visualizations of M1 and M2 states
Module scoring for canonical pro-inflammatory (M1) versus anti-inflammatory/immunosuppressive (M2) transcriptional signatures revealed that microglia exhibit elevated baseline M1 score profiles compared to macrophages. 
Differential M1 minus M2 scoring projected across UMAP space highlighted continuous spatial gradients of polarization rather than discrete binary phenotypes. Quantitative proportion analysis confirmed that M1-like state activation predominated in both lineages, accounting for 85.4% of microglia and 50.7% of macrophages, whereas M2-like state polarization (24.5%) and intermediate states (24.9%) were substantially more prevalent within bone marrow-derived macrophages.
![Aggregated Interactions Count](figures/02_Annotation/M1_score_M2_score_plot.png)

Figure x. Myeloid subtype characterization and M1/M2 polarization states. Violin plots showing module score distributions for M1 (left) and M2 (right) transcriptomic signatures across myeloid sub-types.   
![Aggregated Interactions Count](figures/02_Annotation/M_%20score_minus_M2_score.png)

Figure x. Myeloid subtype characterization and M1/M2 polarization states. UMAP feature plot displaying the composite polarization score (M1 score minus M2 score), where positive values (blue) denote M1-like skewing and negative values (red) denote M2-like skewing.
![Aggregated Interactions Count](figures/02_Annotation/02_Annotation/Percentage_of_myeolid_cells.png)

Figure x. Myeloid subtype characterization and M1/M2 polarization states. Stacked bar plot illustrating the relative proportion of M1-like, M2-like, and intermediate polarization states among Macrophage and Microglia lineages.


#### Differential Expression Analysis and Volcano Plot Visualization
A customized volcano plot visualization framework effectively segregated statistically significant genes (FDR < 0.01 and log_2FC > 1) from minor or non-significant expression changes.

#### Functional Enrichment Analysis (GO, KEGG, Reactome)
Enrichment analysis of directionally separated differentially expressed genes revealed distinct functional programs governed by upregulated and downregulated gene cascades. Over-representation analysis across GO categories (BP, MF, CC), KEGG pathways, and Reactome pathways successfully identified significant biological terms (FDR < 0.05).

#### Pathway Enrichment Visualization
Bidirectional diverging bar plots provided a clear visual separation of upregulated and downregulated functional pathways. Terms associated with upregulated genes exhibited positive normalized enrichment scores, while downregulated pathways extended along negative values.

#### Differential Expression and Pathway Enrichment of Fibroblasts
Comparison of the fibroblast/pericyte population against all other cell types highlighted a localized transcriptomic signature enriched for extracellular matrix (ECM) remodeling and structural cell functions. Volcano plot analysis identified key upregulated fibroblast markers (FDR < 0.01, log_2FC > 1). Downstream functional profiling across REACTOME, KEGG, and GO sub-ontologies confirmed strong positive enrichment for biological processes related to collagen organization.

#### Differential Expression of M1-like vs. M2-like Macrophages
Transcriptional module scoring successfully stratified the macrophage cluster into distinct functional polarization states, including M1-like, M2-like, and intermediate sub-populations. Direct differential expression analysis contrasting M1-like against M2-like macrophages revealed clear marker separation. 

#### Pathway Enrichment Analysis of Macrophage Polarization States
Functional enrichment analysis of macrophage polarization states confirmed distinct metabolic and signaling profiles between M1-like and M2-like phenotypes. 


#### Overall Communication Network using CellChat Model
Aggregated cell-cell communication analysis across all signaling pathways revealed extensive crosstalk among microenvironmental cell populations. Overall interaction counts and weights were heavily dominated by bidirectional signaling between Macrophage_M1-like and Macrophage_M2-like subsets. Pericyte_fibroblasts additionally exhibited robust paracrine communication with both macrophage polarization states, whereas Microglia populations displayed comparatively fewer total interactions across the network
![Aggregated Interactions Count](figures/04_Cell_Cell_Communication/33.2.1%20Circle%20Plot-Overall%20Communication%20Network/Aggregated_Interactions_Count.png)

Figure x. Global cell-cell communication network across the brain metastasis microenvironment. Circular network plot displaying the total number of inferred interactions (counts) between pericytes/fibroblasts, polarized macrophages (M1-like and M2-like), and polarized microglia (M1-like and M2-like). Numbers on directed edges indicate total interaction counts between source and target populations.
![Aggregated Interactions Count](figures/04_Cell_Cell_Communication/33.2.1%20Circle%20Plot-Overall%20Communication%20Network/Aggregated_Interactions_Strength.png)

Figure x. Global cell-cell communication network across the brain metastasis microenvironment. Circular network plot displaying the total interaction strength (weighted communication probability) aggregated across all signaling pathways. Edge thickness corresponds to signaling strength

#### Specific Pathways Network Inference
Analysis of pathway-specific signaling networks demonstrated distinct cell-type-driven communication axes within the brain metastasis microenvironment. For the COLLAGEN pathway, Pericyte_fibroblasts functioned as the exclusive signal source, directing widespread paracrine interactions toward Macrophage_M1-like, Macrophage_M2-like, Microglia_M1-like, and Microglia_M2-like cells, alongside prominent autocrine signaling. 
Conversely, the Prostaglandin pathway operated primarily through dual signaling hubs: Macrophage_M1-like cells served as a primary sender directing signals to both M1-like and M2-like Microglia alongside autocrine feedback, while Pericyte_fibroblasts independently targeted Macrophage_M1-like, Microglia_M1-like, and Microglia_M2-like populations. Macrophage_M2-like cells showed no active participation in the Prostaglandin network
![Aggregated Interactions Count](figures/04_Cell_Cell_Communication/33.2.2%20Circular%20plots%20for%20Pathways%20Communication/Cell_cell_communications_Including_self_communications_COLLAGEN.png)

Figure x. Pathway-specific cell-cell communication networks. Circular network plot of the COLLAGEN signaling pathway, highlighting Pericyte_fibroblasts as the primary signal sender targeting all myeloid sub-populations alongside autocrine signaling.   (B) Circular network plot of the Prostaglandin signaling pathway, displaying dual signaling hubs centered on Macrophage_M1-like cells and Pericyte_fibroblasts directing paracrine signals toward M1-like and M2-like microglia. 
![Aggregated Interactions Count](figures/04_Cell_Cell_Communication/33.2.2%20Circular%20plots%20for%20Pathways%20Communication/Cell_cell_communications_Including_self_communications_Prostaglandin.png)

Figure x. Pathway-specific cell-cell communication networks. Circular network plot of the Prostaglandin signaling pathway, displaying dual signaling hubs centered on Macrophage_M1-like cells and Pericyte_fibroblasts directing paracrine signals toward M1-like and M2-like microglia. 

#### Ligand-Receptor Interactions
Detailed evaluation of ligand-receptor pairs (p < 0.01) demonstrated significant variation across target cell pairs. Evaluation of fibroblast-derived signaling revealed distinct receptor-binding preferences across recipient myeloid populations, particularly regarding Syndecan-4 (SDC4) interactions. Fibroblasts engaged M1-like macrophages (Pericyte_fibroblast -> Macrophage_M1-like) through collagen–SDC4 pairs (COL1A1–SDC4, COL6A1–SDC4, and COL9A3–SDC4) as well as MDK–SDC4 signaling. Crucially, these SDC4-mediated interactions were completely absent in M2-like macrophages (Pericyte_fibroblast -> Macrophage_M2-like), which interacted with fibroblast-derived collagens exclusively via CD44 (COL1A1–CD44, COL6A1–CD44, COL9A3–CD44). A similar polarization pattern was mirrored in microglia, where M1-like microglia retained low-level COL1A1–SDC4 communication while both M1- and M2-like microglia primarily engaged CD44.
![Aggregated Interactions Count](figures/04_Cell_Cell_Communication/33.2.4%20Ligand-Receptor%20Results%20and%20visualisations_fibroblasts_only/Including_Self_Ligand-receptor_interactions_grouped_pathways_fibro_only.png)

Figure x. Fibroblast-centric ligand-receptor communication landscape.Faceted dot plot displaying significant (p < 0.01) ligand-receptor interactions originating from Pericyte_fibroblasts across target recipient cell types. Color intensity represents communication probability. Collagen and Midkine signaling demonstrate selective SDC4 receptor engagement in M1-like macrophages compared to CD44-predominant binding in M2-like macrophages.


### Discussion

### Conclusion
