datainput_multiple_conditionbased <- function(index_multiple_conditionbased_input, index_multiple_sample_normalization_method, index_m_conditionbased1, index_m_conditionbased2, index_m_conditionbased3, index_m_conditionbased4, index_m_conditionbased5, index_m_conditionbased6, index_m_conditionbased7, index_m_conditionbased8, index_m_conditionbased9, index_m_conditionbased10) {
  multiple_sample_clustering_sample_based <- index_multiple_conditionbased_input
  
  Idents(multiple_sample_clustering_sample_based) <- "condition"
  multiple_sample_clustering_sample_based <-subset(multiple_sample_clustering_sample_based, idents = c(index_m_conditionbased1, index_m_conditionbased2))
  if (index_multiple_sample_normalization_method == "LogNormalize"){
  multiple_condition_differential_expressed <- FindMarkers(multiple_sample_clustering_sample_based, ident.1 = index_m_conditionbased1, ident.2 = index_m_conditionbased2, min.pct = index_m_conditionbased3, logfc.threshold = index_m_conditionbased4, test.use = index_m_conditionbased5, only.pos = index_m_conditionbased6)
 }
  else if(index_multiple_sample_normalization_method == "SCTransform") {
    multiple_sample_clustering_sample_based <- PrepSCTFindMarkers(multiple_sample_clustering_sample_based, assay = "SCT", verbose = TRUE)
  multiple_condition_differential_expressed <- FindMarkers(multiple_sample_clustering_sample_based, ident.1 = index_m_conditionbased1, ident.2 = index_m_conditionbased2, min.pct = index_m_conditionbased3, logfc.threshold = index_m_conditionbased4, test.use = index_m_conditionbased5, only.pos = index_m_conditionbased6, assay = "SCT", recorrect_umi = FALSE) 
  }
  multiple_condition_differential_expressed$gene <- rownames(multiple_condition_differential_expressed)
  #top25
  multiple_condition_differential_expressed$direction = ifelse(multiple_condition_differential_expressed$avg_log2FC > index_m_conditionbased4, "UP", "DOWN")
  #top25_sub  <- multiple_condition_differential_expressed %>% group_by(direction) %>%  slice_max(n = index_m_conditionbased9, order_by = avg_log2FC)
  
  if ((index_m_conditionbased9 == 1 | index_m_conditionbased9 == 2 | index_m_conditionbased9 == 3 | index_m_conditionbased9 == 4 | index_m_conditionbased9 == 5 | index_m_conditionbased9 == 6 | index_m_conditionbased9 == 7 | index_m_conditionbased9 == 8 | index_m_conditionbased9 == 9 | index_m_conditionbased9 == 10 | index_m_conditionbased9 == 11 | index_m_conditionbased9 == 12 | index_m_conditionbased9 == 13 | index_m_conditionbased9 == 14 | index_m_conditionbased9 == 15 | index_m_conditionbased9 == 16 | index_m_conditionbased9 == 17 | index_m_conditionbased9 == 18 | index_m_conditionbased9 == 19 | index_m_conditionbased9 == 20 | index_m_conditionbased9 == 21 | index_m_conditionbased9 == 22 | index_m_conditionbased9 == 23 | index_m_conditionbased9 == 24	| index_m_conditionbased9 == 25)){
  top25_sub  <- multiple_condition_differential_expressed %>% group_by(direction) %>%  slice_max(n = as.numeric(index_m_conditionbased9), order_by = avg_log2FC)
  }
  else if (index_m_conditionbased9 == "gene_name_list") {
  gene_name <- unlist(strsplit(index_m_conditionbased10, ","))
  top25_sub <- multiple_condition_differential_expressed %>% dplyr::filter(gene %in% gene_name)
  }
  
  if (index_m_conditionbased7 == "Dot Plot"){
  plots44 <- DotPlot(multiple_sample_clustering_sample_based, features = rev(as.character(unique(top25_sub$gene))), group.by = index_m_conditionbased8) + coord_flip() + theme(legend.position = 'right')
  }
  else if (index_m_conditionbased7 == "VlnPlot"){
  plots44 <- VlnPlot(multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, group.by = index_m_conditionbased8)  + theme(legend.position = 'right')
    }
  else if (index_m_conditionbased7 == "RidgePlot"){
  plots44 <- RidgePlot(multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, group.by = index_m_conditionbased8)
  }
  else if (index_m_conditionbased7 == "FeaturePlot"){
  Idents(multiple_sample_clustering_sample_based) <- index_m_conditionbased8
  plots44 <- FeaturePlot(multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, split.by = index_m_conditionbased8)+ theme(legend.position = 'right')
  }
  else if (index_m_conditionbased7 == "VolcanoPlot"){
  plots44 <- EnhancedVolcano(multiple_condition_differential_expressed, rownames(multiple_condition_differential_expressed), x ="avg_log2FC", y ="p_val_adj")
  }
  return(list(plot1 = plots44, data1 = multiple_condition_differential_expressed, data2 = multiple_sample_clustering_sample_based)) 
}
  
  
  