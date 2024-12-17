datainput_subclustering_multiple_conditionbased <- function(index_subclustering_multiple_conditionbased_input, index_subclustering_multiple_sample_normalization_method, index_m_subclustering_conditionbased1, index_m_subclustering_conditionbased2, index_m_subclustering_conditionbased3, index_m_subclustering_conditionbased4, index_m_subclustering_conditionbased5, index_m_subclustering_conditionbased6, index_m_subclustering_conditionbased7, index_m_subclustering_conditionbased8, index_m_subclustering_conditionbased9) {
  subclustering_multiple_sample_clustering_sample_based <- index_subclustering_multiple_conditionbased_input
  
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
  top25_sub  <- subclustering_multiple_condition_differential_expressed %>% group_by(direction) %>%  slice_max(n = index_m_subclustering_conditionbased9, order_by = avg_log2FC)
  
  if (index_m_subclustering_conditionbased7 == "Dot Plot"){
  plots44 <- DotPlot(subclustering_multiple_sample_clustering_sample_based, features = rev(as.character(unique(top25_sub$gene))), group.by = index_m_subclustering_conditionbased8) + coord_flip() + theme(legend.position = 'right')
  }
  else if (index_m_subclustering_conditionbased7 == "VlnPlot"){
  plots44 <- VlnPlot(subclustering_multiple_sample_clustering_sample_based, features = as.character(unique(top25_sub$gene)), ncol = 5, group.by = index_m_subclustering_conditionbased8)  + theme(legend.position = 'right')
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
  
  
  