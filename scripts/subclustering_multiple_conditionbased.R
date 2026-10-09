datainput_subclustering_multiple_conditionbased <- function(index_subclustering_multiple_conditionbased_input, index_subclustering_multiple_sample_normalization_method, index_m_subclustering_conditionbased1, index_m_subclustering_conditionbased2, index_m_subclustering_conditionbased3, index_m_subclustering_conditionbased4, index_m_subclustering_conditionbased5, index_m_subclustering_conditionbased6, index_m_subclustering_conditionbased7, index_m_subclustering_conditionbased8, index_m_subclustering_conditionbased9, index_m_subclustering_conditionbased10, index_m_subclustering_conditionbased_violin_dots = TRUE, index_m_subclustering_conditionbased_cell_column = NULL, index_m_subclustering_conditionbased_cells = NULL) {
  subclustering_multiple_sample_clustering_sample_based <- index_subclustering_multiple_conditionbased_input
  violin_pt_size <- ifelse(isTRUE(as.logical(index_m_subclustering_conditionbased_violin_dots)), 1, 0)
  # Optional: restrict the analysis to the selected cluster(s) / cell type(s).
  # All cells are used when every group (or no valid group) is selected, as before.
  if (!is.null(index_m_subclustering_conditionbased_cell_column) && isTRUE(index_m_subclustering_conditionbased_cell_column %in% colnames(subclustering_multiple_sample_clustering_sample_based@meta.data))) {
    cell_groups <- as.character(subclustering_multiple_sample_clustering_sample_based@meta.data[[index_m_subclustering_conditionbased_cell_column]])
    all_groups <- unique(cell_groups[!is.na(cell_groups)])
    selected_groups <- intersect(as.character(index_m_subclustering_conditionbased_cells), all_groups)
    if (length(selected_groups) > 0 && length(selected_groups) < length(all_groups)) {
      subclustering_multiple_sample_clustering_sample_based <- subset(subclustering_multiple_sample_clustering_sample_based, cells = colnames(subclustering_multiple_sample_clustering_sample_based)[cell_groups %in% selected_groups])
      cells_per_condition <- table(factor(as.character(subclustering_multiple_sample_clustering_sample_based@meta.data$condition), levels = c(index_m_subclustering_conditionbased1, index_m_subclustering_conditionbased2)))
      if (any(cells_per_condition < 3)) {
        stop(paste0("The selected cluster(s) / cell type(s) have fewer than 3 cells in one of the two conditions (",
                    paste(names(cells_per_condition), as.integer(cells_per_condition), sep = ": ", collapse = ", "),
                    "). Please select more clusters / cell types."))
      }
    }
  }
  
  Idents(subclustering_multiple_sample_clustering_sample_based) <- "condition"
  subclustering_multiple_sample_clustering_sample_based <-subset(subclustering_multiple_sample_clustering_sample_based, idents = c(index_m_subclustering_conditionbased1, index_m_subclustering_conditionbased2))
  if (index_subclustering_multiple_sample_normalization_method == "LogNormalize"){
  subclustering_multiple_condition_differential_expressed <- FindMarkers(subclustering_multiple_sample_clustering_sample_based, ident.1 = index_m_subclustering_conditionbased1, ident.2 = index_m_subclustering_conditionbased2, min.pct = index_m_subclustering_conditionbased3, logfc.threshold = index_m_subclustering_conditionbased4, test.use = index_m_subclustering_conditionbased5, only.pos = index_m_subclustering_conditionbased6)
 }
  else if(index_subclustering_multiple_sample_normalization_method == "SCTransform") {
    subclustering_multiple_sample_clustering_sample_based <- PrepSCTFindMarkers(subclustering_multiple_sample_clustering_sample_based, assay = "SCT", verbose = TRUE)
  subclustering_multiple_condition_differential_expressed <- FindMarkers(subclustering_multiple_sample_clustering_sample_based, ident.1 = index_m_subclustering_conditionbased1, ident.2 = index_m_subclustering_conditionbased2, min.pct = index_m_subclustering_conditionbased3, logfc.threshold = index_m_subclustering_conditionbased4, test.use = index_m_subclustering_conditionbased5, only.pos = index_m_subclustering_conditionbased6, assay = "SCT", recorrect_umi = FALSE) 
  }
  subclustering_multiple_condition_differential_expressed$gene <- rownames(subclustering_multiple_condition_differential_expressed)
  #top25
  subclustering_multiple_condition_differential_expressed$direction = ifelse(subclustering_multiple_condition_differential_expressed$avg_log2FC > index_m_subclustering_conditionbased4, "UP", "DOWN")
  #top25_sub  <- subclustering_multiple_condition_differential_expressed %>% group_by(direction) %>%  slice_max(n = index_m_subclustering_conditionbased9, order_by = avg_log2FC)
  
  if ((index_m_subclustering_conditionbased9 == 1 | index_m_subclustering_conditionbased9 == 2 | index_m_subclustering_conditionbased9 == 3 | index_m_subclustering_conditionbased9 == 4 | index_m_subclustering_conditionbased9 == 5 | index_m_subclustering_conditionbased9 == 6 | index_m_subclustering_conditionbased9 == 7 | index_m_subclustering_conditionbased9 == 8 | index_m_subclustering_conditionbased9 == 9 | index_m_subclustering_conditionbased9 == 10 | index_m_subclustering_conditionbased9 == 11 | index_m_subclustering_conditionbased9 == 12 | index_m_subclustering_conditionbased9 == 13 | index_m_subclustering_conditionbased9 == 14 | index_m_subclustering_conditionbased9 == 15 | index_m_subclustering_conditionbased9 == 16 | index_m_subclustering_conditionbased9 == 17 | index_m_subclustering_conditionbased9 == 18 | index_m_subclustering_conditionbased9 == 19 | index_m_subclustering_conditionbased9 == 20 | index_m_subclustering_conditionbased9 == 21 | index_m_subclustering_conditionbased9 == 22 | index_m_subclustering_conditionbased9 == 23 | index_m_subclustering_conditionbased9 == 24	| index_m_subclustering_conditionbased9 == 25)){
  top25_sub  <-subclustering_multiple_condition_differential_expressed %>% group_by(direction) %>%  slice_max(n = as.numeric(index_m_subclustering_conditionbased9), order_by = avg_log2FC)
  }
  else if (index_m_subclustering_conditionbased9 == "gene_name_list") {
  gene_name <- unlist(strsplit(index_m_subclustering_conditionbased10, ","))
  top25_sub <- subclustering_multiple_condition_differential_expressed %>% dplyr::filter(gene %in% gene_name)
  }
  
  if (index_m_subclustering_conditionbased7 == "Dot Plot"){
  plots44 <- DotPlot(subclustering_multiple_sample_clustering_sample_based, features = rev(as.character(unique(top25_sub$gene))), group.by = index_m_subclustering_conditionbased8) + coord_flip() + theme(legend.position = 'right')
  }
  else if (index_m_subclustering_conditionbased7 == "VlnPlot"){
  plots44 <- VlnPlot(subclustering_multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, group.by = index_m_subclustering_conditionbased8, pt.size = violin_pt_size)  + theme(legend.position = 'right')
    }
  else if (index_m_subclustering_conditionbased7 == "RidgePlot"){
  plots44 <- RidgePlot(subclustering_multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, group.by = index_m_subclustering_conditionbased8)
  }
  else if (index_m_subclustering_conditionbased7 == "FeaturePlot"){
  Idents(subclustering_multiple_sample_clustering_sample_based) <- index_m_subclustering_conditionbased8
  plots44 <- FeaturePlot(subclustering_multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, split.by = index_m_subclustering_conditionbased8)+ theme(legend.position = 'right')
  }
  else if (index_m_subclustering_conditionbased7 == "VolcanoPlot"){
  plots44 <- EnhancedVolcano(subclustering_multiple_condition_differential_expressed, rownames(subclustering_multiple_condition_differential_expressed), x ="avg_log2FC", y ="p_val_adj")
  }
  return(list(plot1 = plots44, data1 = subclustering_multiple_condition_differential_expressed, data2 = subclustering_multiple_sample_clustering_sample_based)) 
}
  
  
  
