#
# Title: "Cell-cell communication for Brain metastases microenvironment"
# Author: "Mostafa Hassanein"
# Date: "2026-09-27"
# Output: R.Script
#  *This script © 2026, by Mostafa Hassanein, is licensed under CC BY 4.0*
#  *To view a copy of this license, visit https://creativecommons.org/licenses/by/4.0/*

# Before run this script, you need to run all scripts of (01 - 03). Note about changing the directory when needed.

# ============================
# 33. Cell-Cell Communication#
# ============================
# 33.1.0 Install & Load Required Libraries#

# 33.1.1 Initialization & Data processing#
  # 33.1.2 Create CellChat object from the existing Seurat object# 
  # 33.1.3 Use the human CellChat database#
  # 33.1.4 Identify Overexpressed Genes & Interactions#
  # 33.1.5 Infer Cell-Cell Communication#
  # 33.1.6 Build the Communication Network#

# 33.2.0 Visualization & Results Extraction#
  # 33.2.1 Circle Plot — Overall Communication Network#
  # 33.2.2 Circular plots for Pathways Communication#
  # 33.2.3 Ligand-Receptor Results and visualisations# 
    # 33.2.3.a Ligand-Receptor Results and visualisations (Including self communication)#
    # 33.2.3.b Ligand-Receptor Results and visualisations (Excluding self communication)#
  # 33.2.4 Ligand-Receptor Results and visualisations_*fibroblasts_only* *For this part, Two adjustments in the functions script are needed*
    # 33.2.4.a Ligand-Receptor Results and visualisations (Including self communication)#
    # 33.2.4.b Ligand-Receptor Results and visualisations (Excluding self communication)#
###########################################
# 33.1.0 Install & Load Required Libraries#
########################################### 

# List of CRAN/Bioconductor packages
cran_pkgs <- c("Seurat", "dplyr", "patchwork", "ggplot2", "tidyseurat", "RColorBrewer", "tibble")

# Check, install missing *CRAN* packages, and load
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Check, install missing CellChat from *GitHub*, and load
if (!requireNamespace("CellChat", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github(c("jinworks/CellChat", "immunogenomics/presto"))
}
library(CellChat)

devtools::install_github('immunogenomics/presto')
library(presto)
##########################################
# 33.1.1 Initialization & Data Processing#
##########################################

# Remove 'Myloid' cells from the first Seurat object:
merged_singlets_filtered <- subset(merged_singlets, subset = cell_type != "Myloid")

# Rename a unified 'DEG_annotation' column in both objects:
names(merged_singlets_filtered@meta.data)[names(merged_singlets_filtered@meta.data) == "cell_type"] <- "DEG_annotation"
names(myeloid@meta.data)[names(myeloid@meta.data) == "annotation_state"] <- "DEG_annotation"

## Merge the two objects (merged_singlets_filtered and myeloid) into one:
full_object <- merge(
  x = myeloid, 
  y = merged_singlets_filtered)

## Join split layers in the RNA assay (essential for Seurat v5):
full_object <- JoinLayers(full_object, assay = "RNA")

## Ensure default assay is set to RNA:
DefaultAssay(full_object) <- "RNA"


## Filter out cells where DEG_annotation is NA:
full_object <- full_object %>% 
  filter(!is.na(DEG_annotation))

## Subset for selected subclustered cells:
fibro_macro_micro_obj <- subset(
  full_object,
  subset = DEG_annotation %in% c("Pericyte_fibroblast", "Macrophage_M1-like", "Macrophage_M2-like", "Microglia_M1-like", "Microglia_M2-like")
)

################################################################
# 33.1.2 Create CellChat object from the existing Seurat object# 
################################################################

### *Initializes the CellChat object using expression data and cell group metadata*
cellchat_obj <- createCellChat(
  object = fibro_macro_micro_obj,
  group.by = "DEG_annotation"  # This is based on DEG_annotation
)
#########################################
# 33.1.3 Use the human CellChat database#
#########################################

### *Assigns human ligand-receptor reference interactions to the object*
cellchat_obj@DB <- CellChatDB.human

## Subset the expression data to save computational cost
### *subsets the normalized expression matrix stored inside the CellChat object to retain only the genes that are part of known ligand-receptor pairs present in the CellChat database (CellChatDB.human).*
### *It discards thousands of other background genes that are not involved in cell-cell communication signaling.*
cellchat_obj <- subsetData(cellchat_obj)

#####################################################
# 33.1.4 Identify Overexpressed Genes & Interactions#
#####################################################

## Identify signaling genes:
###*Finds significantly overexpressed genes in each cell group*
cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)

## Identify ligand-receptor interactions:
### *Filters the database to retain only ligand-receptor pairs where both genes are overexpressed.*
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)

#######################################
# 33.1.5 Infer Cell-Cell Communication#
#######################################

## Compute communication probability (This step may take a few minutes)
###*Calculates the probabilistic communication strength between cell groups.*
cellchat_obj <- computeCommunProb(cellchat_obj)

## Remove interactions from very small cell groups (less than 10 cells)
###*Excludes interactions from cell clusters containing fewer than 10 cells to remove low-confidence noise.*
cellchat_obj <- filterCommunication(
  cellchat_obj,
  min.cells = 10
)

######################################### 
# 33.1.6 Build the Communication Network#
#########################################

## Compute communication probability on a signaling pathway level
### *Combines individual ligand-receptor probabilities to determine signaling strength at the pathway level.*
cellchat_obj <- computeCommunProbPathway(cellchat_obj)

## Aggregate the network (calculate total interactions)
### *Calculates total communication counts and weights across all cell groups to build the macro-network.*
cellchat_obj <- aggregateNet(cellchat_obj)


cellchat_obj_fibro <- subsetCommunication(cellchat_obj, sources.use = "Pericyte_fibroblast")
############################################ 
# 33.2.0 Visualization & Results Extraction#
############################################ 

###############################################################################
# 33.2.1 Circle Plot — Overall Communication Network#
###############################################################################
### Draw the circle network plot
#### Total NUMBER of interactions 
png("Aggregated_Interactions_Count.png", width = 7, height = 10, units = "in", res = 600)
netVisual_circle(
  cellchat_obj@net$count, 
  weight.scale = TRUE, 
  label.edge = TRUE,           # Shows exact counts on the edges
  title.name = "Total Number of Interactions (All Pathways)"
)
dev.off()
#### Total Interaction strenghth(weight)(Probability) 
png("Aggregated_Interactions_Strength.png", width = 7, height = 10, units = "in", res = 600)
netVisual_circle(
  cellchat_obj@net$weight, 
  weight.scale = TRUE, 
  label.edge = FALSE,          
  title.name = "Total Interaction Strength (All Pathways)"
)
dev.off()

###################################################
# 33.2.2 Circular plots for Pathways Communication#
###################################################
###*Extracts a character vector containing the names of all significantly enriched signaling pathways inferred by CellChat and stores them in a variable*
pathways <- cellchat_obj@netP$pathways

### Extract Ligand-Receptor Results
####*Extracts the inferred ligand-receptor interactions into a structured data frame for downstream analysis and export.*
lr_results <- subsetCommunication(cellchat_obj)

#### *function automates the visualization of cell-cell communication pathways by cross-referencing a target list of pathways against a CellChat dataset and arranging circular network plots into a multi-page, 3*3 grid layout. Designed for high-throughput reporting, it dynamically handles output in multi-page PDF, while incorporating robust tryCatch error handling to cleanly skip or label pathways that lack significant interactions without crashing the rendering process.*
#*** Make sure that you already have the "plot_cellchat_pathway_grid" function in your workspace by running "Cell_cell_communication_functions.R" script.

### Run Pathway viusalisation function over all pathways (Including self-communcation) and save the output in the filename
plot_cellchat_pathway_png(
  object = cellchat_obj,
  pathways = pathways,
  filename = "Cell_cell_communications_Including_self_communications"
)

#################################################### 
# 33.2.3 Ligand-Receptor visualisations# 
#################################################### 
## Function to export and visualise the results
#*** Make sure that you already have the "process_lr_interactions" function in your workspace by running "Cell_cell_communication_functions.R" script.
###################################################################################
# 33.2.3.a Ligand-Receptor Results and visualisations (Including self communication)#
###################################################################################
process_lr_interactions(
  cellchat_obj = cellchat_obj, 
  exclude_self = FALSE
)
###################################################################################
# 33.2.3.b Ligand-Receptor Results and visualisations (Excluding self communication)#
###################################################################################
process_lr_interactions(
  cellchat_obj = cellchat_obj, 
  exclude_self = TRUE
)

# View session info directly in the console
sessionInfo()

# Save the session info directly to a text file
capture.output(sessionInfo(), file = "session_cell_cell_communication_info.txt")
