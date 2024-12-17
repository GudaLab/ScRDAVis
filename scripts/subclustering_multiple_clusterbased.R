datainput_subclustering_multiple_clusterbased <- function(index_subclustering_multiple_clusterbased_input, index_subclustering_multiple_clusterbased_features, index_m_subclustering_celltype_method, index_m_subclustering_clusterbased1, index_m_subclustering_clusterbased2, index_m_subclustering_clusterbased3, index_m_subclustering_clusterbased4, index_m_subclustering_clusterbased5) {
  subclustering_multiple_sample_clustering <- index_subclustering_multiple_clusterbased_input
  if (index_m_subclustering_clusterbased1 == "gene_name_list" & index_m_subclustering_clusterbased5 == "NULL"){
    top_markers_plot <- unlist(strsplit(index_m_subclustering_clusterbased2, ","))
    
    if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_clusterbased4) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_clusterbased4)
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_clusterbased4
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5)
    }
	
	else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_celltype_method) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_celltype_method)
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_celltype_method
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5)
    }
  }
 
  else if ((index_m_subclustering_clusterbased1 == 1 | index_m_subclustering_clusterbased1 == 2 | index_m_subclustering_clusterbased1 == 3 | index_m_subclustering_clusterbased1 == 4 | index_m_subclustering_clusterbased1 == 5 | index_m_subclustering_clusterbased1 == 6 | index_m_subclustering_clusterbased1 == 7 | index_m_subclustering_clusterbased1 == 8 | index_m_subclustering_clusterbased1 == 9 | index_m_subclustering_clusterbased1 == 10) &  index_m_subclustering_clusterbased5 == "NULL"){
    index_m_subclustering_clusterbased1 <- as.numeric(index_m_subclustering_clusterbased1)
    subclustering_multiple_sample_clustering_markers <- index_subclustering_multiple_clusterbased_features
    top_markers <- subclustering_multiple_sample_clustering_markers %>% group_by(cluster) %>%  slice_max(n = index_m_subclustering_clusterbased1, order_by = avg_log2FC)
    top_markers_plot <- unique(top_markers$gene)
    if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_clusterbased4) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_clusterbased4)
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_clusterbased4
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5)
    }
	
	else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_celltype_method) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_celltype_method)
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_celltype_method
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5)
    }
	
  }
  else if (index_m_subclustering_clusterbased1 == "gene_name_list" &  (index_m_subclustering_clusterbased5 == "condition"| index_m_subclustering_clusterbased5 == "orig.ident")){
    top_markers_plot <- unlist(strsplit(index_m_subclustering_clusterbased2, ","))
    if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_clusterbased4, split.by = index_m_subclustering_clusterbased5, cols=c(rep("blue",150), "white")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_clusterbased4, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_clusterbased4
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
	
	else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_celltype_method, split.by = index_m_subclustering_clusterbased5, cols=c(rep("blue",150), "white")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_celltype_method, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_celltype_method
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
  }
  else if ((index_m_subclustering_clusterbased1 == 1 | index_m_subclustering_clusterbased1 == 2 | index_m_subclustering_clusterbased1 == 3 | index_m_subclustering_clusterbased1 == 4 | index_m_subclustering_clusterbased1 == 5 | index_m_subclustering_clusterbased1 == 6 | index_m_subclustering_clusterbased1 == 7 | index_m_subclustering_clusterbased1 == 8 | index_m_subclustering_clusterbased1 == 9 | index_m_subclustering_clusterbased1 == 10) &  (index_m_subclustering_clusterbased5 == "condition"| index_m_subclustering_clusterbased5 == "orig.ident")){
    index_m_subclustering_clusterbased1 <- as.numeric(index_m_subclustering_clusterbased1)
    subclustering_multiple_sample_clustering_markers <- index_subclustering_multiple_clusterbased_features
    top_markers <- subclustering_multiple_sample_clustering_markers %>% group_by(cluster) %>%  slice_max(n = index_m_subclustering_clusterbased1, order_by = avg_log2FC)
    top_markers_plot <- unique(top_markers$gene)
    if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_clusterbased4, split.by = index_m_subclustering_clusterbased5, cols=c(rep("blue",150), "white")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_clusterbased4, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
    else if(index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_clusterbased4
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
	
	if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "Dot Plot") {
      plots43 <- DotPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, group.by = index_m_subclustering_celltype_method, split.by = index_m_subclustering_clusterbased5, cols=c(rep("blue",150), "white")) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + coord_flip()
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "VlnPlot") {
      plots43 <- VlnPlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_celltype_method, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
    else if(index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased3 == "FeaturePlot") {
      Idents(subclustering_multiple_sample_clustering) <- index_m_subclustering_celltype_method
      plots43 <- FeaturePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, split.by = index_m_subclustering_clusterbased5)+ theme(legend.position = 'right')
    }
  }
  
  
  if (index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased1 == "gene_name_list" & index_m_subclustering_clusterbased3 == "RidgePlot"){
    top_markers_plot <- unlist(strsplit(index_m_subclustering_clusterbased2, ","))
    plots43 <- RidgePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_clusterbased4)
  }
   else if (index_m_subclustering_clusterbased4 == "seurat_clusters" & (index_m_subclustering_clusterbased1 == 1 | index_m_subclustering_clusterbased1 == 2 | index_m_subclustering_clusterbased1 == 3 | index_m_subclustering_clusterbased1 == 4 | index_m_subclustering_clusterbased1 == 5 | index_m_subclustering_clusterbased1 == 6 | index_m_subclustering_clusterbased1 == 7 | index_m_subclustering_clusterbased1 == 8 | index_m_subclustering_clusterbased1 == 9 | index_m_subclustering_clusterbased1 == 10) &  index_m_subclustering_clusterbased3 == "RidgePlot"){
    index_m_subclustering_clusterbased1 <- as.numeric(index_m_subclustering_clusterbased1)
    subclustering_multiple_sample_clustering_markers <- index_subclustering_multiple_clusterbased_features
    top_markers <- subclustering_multiple_sample_clustering_markers %>% group_by(cluster) %>%  slice_max(n = index_m_subclustering_clusterbased1, order_by = avg_log2FC)
    top_markers_plot <- unique(top_markers$gene)
    plots43 <- RidgePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_clusterbased4)
  } 
  
  else if (index_m_subclustering_clusterbased4 == "predicted" & index_m_subclustering_clusterbased1 == "gene_name_list" & index_m_subclustering_clusterbased3 == "RidgePlot"){
    top_markers_plot <- unlist(strsplit(index_m_subclustering_clusterbased2, ","))
    plots43 <- RidgePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_celltype_method)
  }
  else if (index_m_subclustering_clusterbased4 == "predicted" & (index_m_subclustering_clusterbased1 == 1 | index_m_subclustering_clusterbased1 == 2 | index_m_subclustering_clusterbased1 == 3 | index_m_subclustering_clusterbased1 == 4 | index_m_subclustering_clusterbased1 == 5 | index_m_subclustering_clusterbased1 == 6 | index_m_subclustering_clusterbased1 == 7 | index_m_subclustering_clusterbased1 == 8 | index_m_subclustering_clusterbased1 == 9 | index_m_subclustering_clusterbased1 == 10) &  index_m_subclustering_clusterbased3 == "RidgePlot"){
    index_m_subclustering_clusterbased1 <- as.numeric(index_m_subclustering_clusterbased1)
    subclustering_multiple_sample_clustering_markers <- index_subclustering_multiple_clusterbased_features
    top_markers <- subclustering_multiple_sample_clustering_markers %>% group_by(cluster) %>%  slice_max(n = index_m_subclustering_clusterbased1, order_by = avg_log2FC)
    top_markers_plot <- unique(top_markers$gene)
    plots43 <- RidgePlot(subclustering_multiple_sample_clustering, features = top_markers_plot, ncol = 5, group.by = index_m_subclustering_celltype_method)
  } 
  
  
  if (index_m_subclustering_clusterbased4 == "seurat_clusters"){
    if(index_m_subclustering_clusterbased5 == 'NULL'){
      subclustering_multiple_sample_cell_count_proportion <- PrctCellExpringGene(subclustering_multiple_sample_clustering, genes = top_markers_plot)
    }
    else if(index_m_subclustering_clusterbased5 == 'condition'){
      subclustering_multiple_sample_cell_count_proportion <- PrctCellExpringGene(subclustering_multiple_sample_clustering, genes = top_markers_plot, group.by = 'condition')
    }
    else if(index_m_subclustering_clusterbased5 == 'orig.ident'){
    subclustering_multiple_sample_cell_count_proportion <- PrctCellExpringGene(subclustering_multiple_sample_clustering, genes = top_markers_plot, group.by = 'orig.ident')
    }
  }
  else if (index_m_subclustering_clusterbased4 == "seurat_clusters" & index_m_subclustering_clusterbased3 == "RidgePlot"){
    subclustering_multiple_sample_cell_count_proportion <- PrctCellExpringGene(subclustering_multiple_sample_clustering, genes = top_markers_plot)
  }
  else{
    subclustering_multiple_sample_cell_count_proportion <-  NA 
  }
  return(list(plot1 = plots43, data1 = subclustering_multiple_sample_clustering, data2 = subclustering_multiple_sample_cell_count_proportion))
}