# Function: Circular plots for Pathways Communication (33.2.2)
plot_cellchat_pathway_png <- function(object, 
                                      pathways, 
                                      filename = "CellChat_Pathway",
                                      width = 10, 
                                      height = 10,
                                      res = 300) {
  # Filter pathways to only those present in the CellChat object
  valid_pathways <- intersect(pathways, object@netP$pathways)
  if (length(valid_pathways) == 0) {
    stop("None of the specified pathways were found in the CellChat object.")
  }
  
  # Iterate through each pathway and save to individual PNG files
  for (i in seq_along(valid_pathways)) {
    p <- valid_pathways[i]
    
    # Generate unique filename incorporating pathway name
    png_filename <- paste0(filename, "_", p, ".png")
    
    png(filename = png_filename, 
        width = width, 
        height = height, 
        units = "in", 
        res = res)
    
    # Set standard margin parameters
    par(mfrow = c(1, 1), mar = c(2, 2, 3, 2) + 0.1, xpd = TRUE)
    
    # Render pathway visualization (adjust layout argument if needed)
    netVisual_aggregate(object, signaling = p, layout = "circle")
    
    dev.off()
  }
  
  cat("Finished! Saved", length(valid_pathways), "pathway PNG file(s).\n")
}

################################################################################

# Function for Ligand-Receptor Results and visualisations (33.2.3)
process_lr_interactions <- function(cellchat_obj, exclude_self = FALSE, file_prefix = "") {
  # Extract Ligand-Receptor Results
  lr_results <- subsetCommunication(cellchat_obj)
  # Filter out self-communication if requested
  if (exclude_self) {
    # Remove rows where source and target are the same cell type
    lr_results <- lr_results %>% filter(source != target)
    prefix_str <- paste0(file_prefix, "Excluding_Self_")
    message("Processing: EXCLUDING self-communication...")
  } else {
    prefix_str <- paste0(file_prefix, "Including_Self_")
    message("Processing: INCLUDING self-communication...")
  }
  #  Export lr_results as a csv file
  csv_filename <- paste0(prefix_str, "Cell_cell_communications_Ligand_Receptor_pathways_data.csv")
  write.csv(lr_results, file = csv_filename, row.names = FALSE)
  #  Visualisations of LR interactions *without* pathways
  plot_data <- lr_results %>%
    mutate(Cell_pair = paste(source, "->", target)) %>%
    filter(pval < 0.01) 
  # %>%    #Uncomment these 2 lines for fibroblast only visualisations
  #   filter(source == "Pericyte_fibroblast") 
  # Safety check in case filtering removes all rows
  if (nrow(plot_data) > 0) {
    # png(paste0(prefix_str, "Ligand-receptor_interactions_fibro_only.png"), width = 4, height = 8, units = "in", res = 1200) #Uncomment this line for fibroblast only visualisations
    png(paste0(prefix_str, "Ligand-receptor_interactions.png"), width = 8, height = 16, units = "in", res = 1200)
    p1 <- ggplot(plot_data, aes(x = Cell_pair, y = interaction_name_2)) +
      geom_point(aes(color = prob), size = 3.5) +
      scale_color_gradientn(
        colors = rev(brewer.pal(11, "Spectral")),
        name = "Commun. Prob.",
        breaks = range(plot_data$prob),
        labels = c("min", "max"),
        guide = guide_colorbar(ticks = FALSE, barheight = 5)
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(face = "italic", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.key = element_blank()
      )
    print(p1) # Explicitly print inside function
    dev.off()
  } else {
    message("No significant interactions found for Plot 1. Skipping.")
  }
  # Visualisations of LR interactions *with* pathways
  plot_data_with_pathways <- lr_results %>%
    mutate(Cell_pair = paste(source, "->", target)) %>%
    filter(pval < 0.01) %>%
  # %>% #Uncomment these 2 lines for fibroblast only visualisations
  #   filter(source == "Pericyte_fibroblast") %>% 
    mutate(pathway_label = paste(pathway_name, "Pathway through", annotation)) %>%
    arrange(pathway_label, interaction_name_2) %>%
    mutate(
      interaction_name_2 = factor(interaction_name_2, levels = rev(unique(interaction_name_2))),
      pathway_label = factor(pathway_label, levels = unique(pathway_label))
    )
  if (nrow(plot_data_with_pathways) > 0) {
    # png(paste0(prefix_str, "Ligand-receptor_interactions_grouped_pathways_fibro_only.png"), width = 7, height = 9, units = "in", res = 1200) #Uncomment this line for fibroblast only visualisations
    png(paste0(prefix_str, "Ligand-receptor_interactions_grouped_pathways.png"), width = 12, height = 18, units = "in", res = 1200)
    p2 <- ggplot(plot_data_with_pathways, aes(x = Cell_pair, y = interaction_name_2)) +
      geom_point(aes(color = prob), size = 3.5) +
      facet_grid(pathway_label ~ ., scales = "free_y", space = "free_y") +
      scale_color_gradientn(
        colors = rev(brewer.pal(11, "Spectral")),
        name = "Commun. Prob.",
        breaks = range(plot_data_with_pathways$prob),
        labels = c("min", "max"),
        guide = guide_colorbar(
          ticks = FALSE, 
          barwidth = 0.3, 
          barheight = 3, 
          title.position = "top",
          title.hjust = 0.5
        )
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.text.y = element_text(face = "italic", color = "black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(angle = 0, face = "plain", color = "black", size = 8, hjust = 0),
        strip.background = element_rect(fill = "grey95", color = "black", linewidth = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        panel.grid.major = element_line(color = "grey92"),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        legend.key = element_blank(),
        legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 8)
      )
    print(p2) # Explicitly print inside function
    dev.off()
  } else {
    message("No significant interactions found for Plot 2. Skipping.")
  }
  message("Done! Files saved successfully.\n")
}
