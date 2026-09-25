
library(ggplot2)
library(SingleR)
library(celldex)
set.seed(42)

#=============================================================================
# 20. Automated Annotation of merged_singlets
#=============================================================================
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
# 21. Differential expression annotation of merged-singlets
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
#=============================================================================
# 22. Create a heatmap and Check for NAs first
# =============================================================================
sum(is.na(merged_singlets$cell_type))
table(merged_singlets$cell_type, useNA = "always")
cells_keep <- colnames(merged_singlets)[!is.na(merged_singlets$cell_type)]


DoHeatmap(
  merged_singlets,
  features = genes_to_plot,
  cells = cells_keep,
  group.by = "cell_type",
  label = FALSE,       # removes text labels, keeps the colored strip
  size = 4,
  angle = 45,
  hjust = 0
) +
  scale_fill_gradientn(colors = c("navy", "white", "firebrick")) +
  theme(axis.text.y = element_text(size = 4, face = "bold"))

#===============================================================================
#23. Differential expression annotation of Myeloid cells
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
# 24. Determining of M1 and M2 states
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
# 25. Annotation of M1 and M2 states
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
# 26. Different visualizations of M1 and M2 states
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

mm_summary <- mm_cells %>%
  count(cell_type_mm, state) %>%
  group_by(cell_type_mm) %>%
  mutate(pct = n / sum(n) * 100,
         label = paste0(round(pct, 1), "%"))

ggplot(mm_summary, aes(x = cell_type_mm, y = n, fill = state)) +
  geom_bar(stat = "identity", position = "fill") +
  geom_text(aes(label = label),
            position = position_fill(vjust = 0.5),
            size = 3.5, color = "white", fontface = "bold") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Proportion of cells", fill = "State") +
  theme_classic()

round(prop.table(table(mm_cells$cell_type_mm, mm_cells$state), margin = 1) * 100, 1)


# Plot 4 

DimPlot(myeloid, reduction = "umap", group.by = "cell_type_mm",
        label = TRUE, repel = TRUE) + ggtitle("Myeloid cell types")


top_markers_heatmap <- markers %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%   # 5 instead of 10 -> ~75 rows total
  ungroup()








