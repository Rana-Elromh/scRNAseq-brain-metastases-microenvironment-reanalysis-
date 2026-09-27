
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(stringr)


# ---------------------------------------------------------
# 27. VOLCANO PLOT — cleaner labeling, capped y-axis, up/down labeled separately
# ---------------------------------------------------------

run_volcano_plot <- function(seurat_obj, group_col, ident.1, ident.2 = NULL,
                             title = "Volcano Plot", label_n_each = 12,
                             pval_cap = 300, min_pct_label = 0.1,
                             fc_label_cap = 6) {
  
  Idents(seurat_obj) <- group_col
  
  degs <- FindMarkers(
    seurat_obj, ident.1 = ident.1, ident.2 = ident.2,
    logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox"
  )
  degs$gene <- rownames(degs)
  
  degs$category <- "Not significant"
  degs$category[degs$p_val_adj < 0.05] <- "FDR<0.05"
  degs$category[degs$p_val_adj < 0.01 & abs(degs$avg_log2FC) > 1] <- "FDR<0.01 & |log2FC|>1"
  degs$category <- factor(
    degs$category,
    levels = c("FDR<0.01 & |log2FC|>1", "FDR<0.05", "Not significant")
  )
  
  degs$neglog10p <- -log10(degs$p_val)
  was_capped <- max(degs$neglog10p, na.rm = TRUE) > pval_cap
  degs$neglog10p_plot <- pmin(degs$neglog10p, pval_cap)
  
  # cap fold change just for the labeling step -- keeps label picks away from
  # unstable, lowly-expressed extreme-FC outliers (the ones cluttering your plot)
  sig <- degs %>%
    filter(category == "FDR<0.01 & |log2FC|>1") %>%
    filter(pct.1 > min_pct_label | pct.2 > min_pct_label) %>%   # expressed in a meaningful fraction of cells
    filter(abs(avg_log2FC) < fc_label_cap)                       # drop extreme unstable FC from label picks
  
  # rank by significance (not raw FC) for label selection -- this is what
  # the published plot does: labels go to genes that are both significant
  # AND consistently expressed, not to lucky low-count outliers
  top_up   <- sig %>% filter(avg_log2FC > 0) %>% arrange(p_val, desc(avg_log2FC)) %>% slice_head(n = label_n_each)
  top_down <- sig %>% filter(avg_log2FC < 0) %>% arrange(p_val, avg_log2FC)       %>% slice_head(n = label_n_each)
  top_labels <- bind_rows(top_up, top_down)
  
  y_lab <- "-log10(p-value)"
  if (was_capped) y_lab <- paste0(y_lab, "  [capped at ", pval_cap, "]")
  
  p <- ggplot(degs, aes(x = avg_log2FC, y = neglog10p_plot, color = category)) +
    geom_point(alpha = 0.7, size = 0.9) +
    scale_color_manual(values = c(
      "FDR<0.01 & |log2FC|>1" = "red3",
      "FDR<0.05" = "steelblue",
      "Not significant" = "grey75"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.4) +
    geom_text_repel(
      data = top_labels, aes(label = gene), size = 3, color = "black",
      max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2,
      segment.size = 0.25, segment.color = "grey40",
      min.segment.length = 0, force = 2, force_pull = 0.5,
      max.iter = 20000, seed = 42
    ) +
    theme_classic(base_size = 13) +
    labs(title = title, x = "log2(Fold Change)", y = y_lab, color = NULL) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "top",
      legend.title = element_blank()
    )
  
  list(degs = degs, plot = p)
}

result <- run_volcano_plot(
  seurat_obj = merged_singlets,      # or myeloid, whatever object has your fibroblast cluster
  group_col  = "seurat_clusters",    # or "cell_type" if already annotated
  ident.1    = "11",                 # your fibroblast cluster
)

print(result$plot)
ggsave(
  "fibroblast_volcano.png",
  plot = result$plot,
  width = 8, height = 7, dpi = 300
)

# ---------------------------------------------------------
# 28. ENRICHMENT — run on UP and DOWN genes separately, GO split into BP/MF/CC
# ---------------------------------------------------------
run_enrichment_updown <- function(degs, pval_cutoff = 0.05, logfc_cutoff = 1) {
  
  up_genes   <- degs$gene[degs$p_val_adj < pval_cutoff & degs$avg_log2FC >  logfc_cutoff]
  down_genes <- degs$gene[degs$p_val_adj < pval_cutoff & degs$avg_log2FC < -logfc_cutoff]
  
  entrez_up   <- bitr(up_genes,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  entrez_down <- bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  
  # GO, split by ontology (BP / MF / CC) — use ont = "ALL" and separate later
  go_up_all   <- enrichGO(entrez_up$ENTREZID,   OrgDb = org.Hs.eg.db, ont = "ALL",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  go_down_all <- enrichGO(entrez_down$ENTREZID, OrgDb = org.Hs.eg.db, ont = "ALL",
                          pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
  
  kegg_up   <- enrichKEGG(entrez_up$ENTREZID,   organism = "hsa", pvalueCutoff = 0.05)
  kegg_down <- enrichKEGG(entrez_down$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
  
  reactome_up   <- enrichPathway(entrez_up$ENTREZID,   organism = "human", pvalueCutoff = 0.05, readable = TRUE)
  reactome_down <- enrichPathway(entrez_down$ENTREZID, organism = "human", pvalueCutoff = 0.05, readable = TRUE)
  
  list(
    go_up = go_up_all, go_down = go_down_all,
    kegg_up = kegg_up, kegg_down = kegg_down,
    reactome_up = reactome_up, reactome_down = reactome_down
  )
}


# ---------------------------------------------------------
# 29. DIVERGING BAR PLOT — orange = upregulated pathways, blue = downregulated
#    (matches Fig 3c/d style from the paper)
# ---------------------------------------------------------
plot_enrichment_updown <- function(res_up, res_down, top_n = 10,
                                   title = "Pathway enrichment",
                                   ontology_filter = NULL) {
  
  prep <- function(res, direction, sign) {
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (!is.null(ontology_filter) && "ONTOLOGY" %in% colnames(df)) {
      df <- df %>% filter(ONTOLOGY == ontology_filter)
    }
    if (nrow(df) == 0) return(NULL)
    df %>%
      slice_min(order_by = p.adjust, n = top_n) %>%
      mutate(direction = direction, score = sign * (-log10(p.adjust) / 10))
  }
  
  df_up   <- prep(res_up,   "Upregulated",   1)
  df_down <- prep(res_down, "Downregulated", -1)
  df <- bind_rows(df_up, df_down)
  
  if (is.null(df) || nrow(df) == 0) {
    message("No significant terms for: ", title)
    return(NULL)
  }
  
  df <- df %>% mutate(Description = factor(Description, levels = Description[order(score)]))
  
  ggplot(df, aes(x = score, y = Description, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("Upregulated" = "#E67E22", "Downregulated" = "#2E86C1")) +
    theme_classic(base_size = 12) +
    labs(x = "-log10(adjusted P-value)/10", y = NULL, fill = NULL, title = title)
}


# =========================================================
# 30. USAGE — Fibroblast, with up+down together everywhere
# =========================================================

fibro_label <- "Pericyte_fibroblast"
Idents(merged_singlets) <- "cell_type"

fibro_res <- run_volcano_plot(
  merged_singlets, group_col = "cell_type",
  ident.1 = fibro_label, ident.2 = NULL,
  title = "Fibroblast (cluster 11) vs Rest"
)
fibro_res$plot
ggsave("fibroblast_volcano.png", plot = fibro_res$plot, width = 9, height = 7, dpi = 300)
write.csv(fibro_res$degs, "fibroblast_DEGs.csv", row.names = FALSE)

fibro_enrich <- run_enrichment_updown(fibro_res$degs)

# REACTOME — up (orange) + down (blue) on one diverging plot
p_reactome <- plot_enrichment_updown(fibro_enrich$reactome_up, fibro_enrich$reactome_down,
                                     title = "REACTOME pathway (Fibroblast)")
p_reactome
ggsave("fibroblast_REACTOME.png", plot = p_reactome, width = 8, height = 6, dpi = 300)

# KEGG — up + down
p_kegg <- plot_enrichment_updown(fibro_enrich$kegg_up, fibro_enrich$kegg_down,
                                 title = "KEGG pathway (Fibroblast)")
p_kegg
ggsave("fibroblast_KEGG.png", plot = p_kegg, width = 8, height = 6, dpi = 300)

# GO — three separate plots: BP, MF, CC (each with up + down)
fibro_enrich$go_up@result$Description   <- str_to_sentence(fibro_enrich$go_up@result$Description)
fibro_enrich$go_down@result$Description <- str_to_sentence(fibro_enrich$go_down@result$Description)

p_go_bp <- plot_enrichment_updown(fibro_enrich$go_up, fibro_enrich$go_down,
                                  title = "GO Biological Process (Fibroblast)",
                                  ontology_filter = "BP")
p_go_mf <- plot_enrichment_updown(fibro_enrich$go_up, fibro_enrich$go_down,
                                  title = "GO Molecular Function (Fibroblast)",
                                  ontology_filter = "MF")

p_go_cc <- plot_enrichment_updown(fibro_enrich$go_up, fibro_enrich$go_down,
                                  title = "GO Cellular Component (Fibroblast)",
                                  ontology_filter = "CC")
print(p_go_bp)
print(p_go_mf)
print(p_go_cc)

p_go_bp; p_go_mf; p_go_cc
ggsave("Fibroblast_GO_BP.png", plot = p_go_bp, width = 8, height = 6, dpi = 300)
ggsave("Fibroblast_GO_MF.png", plot = p_go_mf, width = 8, height = 6, dpi = 300)
ggsave("Fibroblast_GO_CC.png", plot = p_go_cc, width = 8, height = 6, dpi = 300)

# ======================================================================================
# 31.M1-macrophage CELL COUNT CHECK + FULL PIPELINE (volcano + GO/KEGG/REACTOME)
# ======================================================================================
# 1. Define/update the function (paste this in first)
run_volcano_plot <- function(seurat_obj, group_col, ident.1, ident.2 = NULL,
                             title = "Volcano Plot", label_n_each = 12,
                             pval_cap = 300, min_pct_label = 0.1,
                             fc_label_cap = 6) {
  
  Idents(seurat_obj) <- group_col
  
  degs <- FindMarkers(
    seurat_obj, ident.1 = ident.1, ident.2 = ident.2,
    logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox"
  )
  degs$gene <- rownames(degs)
  
  degs$category <- "Not significant"
  degs$category[degs$p_val_adj < 0.05] <- "FDR<0.05"
  degs$category[degs$p_val_adj < 0.01 & abs(degs$avg_log2FC) > 1] <- "FDR<0.01 & |log2FC|>1"
  degs$category <- factor(
    degs$category,
    levels = c("FDR<0.01 & |log2FC|>1", "FDR<0.05", "Not significant")
  )
  
  degs$neglog10p <- -log10(degs$p_val)
  was_capped <- max(degs$neglog10p, na.rm = TRUE) > pval_cap
  degs$neglog10p_plot <- pmin(degs$neglog10p, pval_cap)
  
  sig <- degs %>%
    filter(category == "FDR<0.01 & |log2FC|>1") %>%
    filter(pct.1 > min_pct_label | pct.2 > min_pct_label) %>%
    filter(abs(avg_log2FC) < fc_label_cap)
  
  top_up   <- sig %>% filter(avg_log2FC > 0) %>% arrange(p_val, desc(avg_log2FC)) %>% slice_head(n = label_n_each)
  
  top_down <- sig %>% filter(avg_log2FC < 0) %>% arrange(p_val, avg_log2FC)       %>% slice_head(n = label_n_each)
  top_labels <- bind_rows(top_up, top_down)
  
  y_lab <- "-log10(p-value)"
  if (was_capped) y_lab <- paste0(y_lab, "  [capped at ", pval_cap, "]")
  
  p <- ggplot(degs, aes(x = avg_log2FC, y = neglog10p_plot, color = category)) +
    geom_point(alpha = 0.7, size = 0.9) +
    scale_color_manual(values = c(
      "FDR<0.01 & |log2FC|>1" = "red3",
      "FDR<0.05" = "steelblue",
      "Not significant" = "grey75"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.4) +
    geom_text_repel(
      data = top_labels, aes(label = gene), size = 3, color = "black",
      max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2,
      segment.size = 0.25, segment.color = "grey40",
      min.segment.length = 0, force = 2, force_pull = 0.5,
      max.iter = 20000, seed = 42
    ) +
    theme_bw(base_size = 13) +
    labs(title = title, x = "log2(Fold Change)", y = y_lab, color = NULL) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
  
  
  list(degs = degs, plot = p)
}

# 1. myeloid_subtype
myeloid_labels <- c(
  "0" = "Macrophage",
  "1" = "Macrophage_MMP9_high_mito",
  "2" = "Microglia",
  "3" = "cDC2",
  "4" = "Proliferating_myeloid"
)

myeloid$myeloid_subtype <- unname(myeloid_labels[as.character(myeloid$seurat_clusters)])
Idents(myeloid) <- "myeloid_subtype"

table(myeloid$myeloid_subtype)


# 2. M1_score1 / M2_score1
m1_genes <- c("CD80","CD86","TNF","IL1B","IL6","IL12A","IL12B","CXCL9","CXCL10",
              "CXCL11","NOS2","STAT1","IRF5","FCGR1A","IDO1","CCL5","HLA-DRA")
m2_genes <- c("CD163","MRC1","MSR1","IL10","TGFB1","ARG1","CCL18","CCL22",
              "CD209","VEGFA","CLEC7A","TGM2","STAB1","MARCO")

m1_genes <- intersect(m1_genes, rownames(myeloid))
m2_genes <- intersect(m2_genes, rownames(myeloid))

myeloid <- AddModuleScore(myeloid, features = list(m1_genes), name = "M1_score")
myeloid <- AddModuleScore(myeloid, features = list(m2_genes), name = "M2_score")


# 3. cell_type_mm, state, annotation_state
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

# 3. Call it on your macrophage comparison
macrophage_res <- run_volcano_plot(
  myeloid,
  group_col = "annotation_state",
  ident.1   = "Macrophage_M1-like",
  ident.2   = "Macrophage_M2-like",
  title     = "M1-like vs M2-like Macrophages"
)

# 4. View it
macrophage_res$plot

# 5. Save it
ggsave("macrophage_M1_vs_M2_volcano_framed.png",
       plot = macrophage_res$plot, width = 9, height = 7, dpi = 300)
  

macrophage_res <- run_volcano_plot(
  myeloid,
  group_col = "annotation_state",
  ident.1   = "Macrophage_M1-like",
  ident.2   = "Macrophage_M2-like",
  title     = "M1-like vs M2-like Macrophages"
)

macrophage_res$plot

ggsave("macrophage_M1_vs_M2_volcano_framed.png",
       plot = macrophage_res$plot, width = 9, height = 7, dpi = 300)
write.csv(macrophage_res$degs, "Macrophage_M1_vs_M2_DEGs.csv", row.names = FALSE)
#===============================================================
# 32. ENRICHMENT
# ===============================================================

enrich <- run_enrichment_updown(macrophage_res$degs)

# REACTOME
p_reactome <- plot_enrichment_updown(
  enrich$reactome_up, enrich$reactome_down,
  title = "REACTOME pathway (M1-like vs M2-like Macrophages)"
)
p_reactome
ggsave("macrophage_REACTOME.png", plot = p_reactome, width = 8, height = 6, dpi = 300)

# KEGG
p_kegg <- plot_enrichment_updown(
  enrich$kegg_up, enrich$kegg_down,
  title = "KEGG pathway (M1-like vs M2-like Macrophages)"
)
p_kegg
ggsave("macrophage_KEGG.png", plot = p_kegg, width = 8, height = 6, dpi = 300)

# GO -- BP / MF / CC

plot_enrichment_updown <- function(res_up, res_down, top_n = 10,
                                   title = "Pathway enrichment",
                                   ontology_filter = NULL) {
  
  prep <- function(res, direction, sign) {
    if (is.null(res)) return(NULL)
    df <- as.data.frame(res)
    if (!is.null(ontology_filter) && "ONTOLOGY" %in% colnames(df)) {
      df <- df %>% filter(ONTOLOGY == ontology_filter)
    }
    if (nrow(df) == 0) return(NULL)
    df %>%
      slice_min(order_by = p.adjust, n = top_n) %>%
      mutate(direction = direction, score = sign * (-log10(p.adjust) / 10))
  }
  
  df_up   <- prep(res_up,   "Upregulated",   1)
  df_down <- prep(res_down, "Downregulated", -1)
  df <- bind_rows(df_up, df_down)
  
  if (is.null(df) || nrow(df) == 0) {
    message("No significant terms for: ", title)
    return(NULL)
  }
  
# make labels unique in case the same term appears in both up and down,
# or is duplicated within one direction -- factor() can't have duplicate levels
  df <- df %>%
    mutate(row_id = row_number()) %>%
    mutate(Description_unique = make.unique(Description))
  
  df <- df %>% mutate(Description_unique = factor(Description_unique, levels = Description_unique[order(score)]))
  
  ggplot(df, aes(x = score, y = Description_unique, fill = direction)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    scale_fill_manual(values = c("Upregulated" = "#E67E22", "Downregulated" = "#2E86C1")) +
    scale_y_discrete(labels = df$Description[order(match(df$Description_unique, levels(df$Description_unique)))]) +
    theme_classic(base_size = 12) +
    labs(x = "-log10(adjusted P-value)/10", y = NULL, fill = NULL, title = title)
}
p_go_bp <- plot_enrichment_updown(enrich$go_up, enrich$go_down,
                                  title = "GO Biological Process (M1-like vs M2-like Macrophages)", ontology_filter = "BP")
p_go_mf <- plot_enrichment_updown(enrich$go_up, enrich$go_down,
                                  title = "GO Molecular Function (M1-like vs M2-like Macrophages)", ontology_filter = "MF")
p_go_cc <- plot_enrichment_updown(enrich$go_up, enrich$go_down,
                                  title = "GO Cellular Component (M1-like vs M2-like Macrophages)", ontology_filter = "CC")

# view each one
print(p_go_bp)
print(p_go_mf)
print(p_go_cc)

ggsave("macrophage_GO_BP.png", plot = p_go_bp, width = 8, height = 6, dpi = 300)
ggsave("macrophage_GO_MF.png", plot = p_go_mf, width = 8, height = 6, dpi = 300)
ggsave("macrophage_GO_CC.png", plot = p_go_cc, width = 8, height = 6, dpi = 300)