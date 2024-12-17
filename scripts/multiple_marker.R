datainput_multiple_marker <- function(index_multiple_marker_input, index_m_marker1, index_m_marker2, index_m_marker3, index_m_marker4, index_m_marker5, index_m_marker6, index_m_marker7, index_m_marker8, index_m_marker9, index_m_marker10, index_multiple_sample_normalization_method){
  multiple_sample_clustering <- index_multiple_marker_input
  index_m_marker5 <-as.logical(index_m_marker5)
  if(index_m_marker1 == 1){
    #All markers
    # multiple_sample_clustering_markers <- FindAllMarkers(multiple_sample_clustering, min.pct = 0.25, logfc.threshold = 0.25, test.use = 'wilcox', only.pos = TRUE)
    if (index_multiple_sample_normalization_method == "LogNormalize"){
    multiple_sample_clustering_markers <- FindAllMarkers(multiple_sample_clustering, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5)
    }
    else if (index_multiple_sample_normalization_method == "SCTransform"){
    multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
    multiple_sample_clustering_markers <- FindAllMarkers(multiple_sample_clustering, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, assay = "SCT", recorrect_umi = FALSE)
    }
    
    top_markers <- multiple_sample_clustering_markers %>% group_by(cluster) %>%  slice_max(n = 5, order_by = avg_log2FC)
    top_markers_plot <- unique(top_markers$gene)
    
    #cluster_heatmap
    plots39 <- DoHeatmap(multiple_sample_clustering, features = top_markers$gene, angle = 90, group.by ="seurat_clusters")
    
    return(list(data1 = multiple_sample_clustering_markers, data2 = multiple_sample_clustering, plot1 = plots39))
  }
  
  else if(index_m_marker1 == 2){
    #Markers in one cluster
    if (index_multiple_sample_normalization_method == "LogNormalize"){
    multiple_sample_clustering_markers_one_cluster <- FindMarkers(multiple_sample_clustering, ident.1 = index_m_marker6, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5)
    multiple_sample_clustering_markers_one_cluster$gene <- rownames(multiple_sample_clustering_markers_one_cluster)
	}
    else if (index_multiple_sample_normalization_method == "SCTransform"){
    multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
    multiple_sample_clustering_markers_one_cluster <- FindMarkers(multiple_sample_clustering, ident.1 = index_m_marker6, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, assay = "SCT", recorrect_umi = FALSE) 
    multiple_sample_clustering_markers_one_cluster$gene <- rownames(multiple_sample_clustering_markers_one_cluster)
	} 
    return(list(data1 = multiple_sample_clustering_markers_one_cluster, data2 = multiple_sample_clustering))
    
  }
  else if(index_m_marker1 == 3){
    #Markers in two cluster
    if (index_multiple_sample_normalization_method == "LogNormalize"){
    multiple_sample_clustering_markers_two_cluster <- FindMarkers(multiple_sample_clustering, ident.1 = index_m_marker6,  ident.2 = index_m_marker7, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5)
    multiple_sample_clustering_markers_two_cluster$gene <- rownames(multiple_sample_clustering_markers_two_cluster)
	}
    else if (index_multiple_sample_normalization_method == "SCTransform"){
    multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
    multiple_sample_clustering_markers_two_cluster <- FindMarkers(multiple_sample_clustering, ident.1 = index_m_marker6,  ident.2 = index_m_marker7, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, assay = "SCT", recorrect_umi = FALSE)
    multiple_sample_clustering_markers_two_cluster$gene <- rownames(multiple_sample_clustering_markers_two_cluster)
	} 
    return(list(data1 = multiple_sample_clustering_markers_two_cluster, data2 = multiple_sample_clustering))
  }
  else if(index_m_marker1 == 4){
    #Conserved Markers in one cluster
    if (index_multiple_sample_normalization_method == "LogNormalize"){
    multiple_sample_clustering_markers_one_cluster <- FindConservedMarkers(multiple_sample_clustering, ident.1 = index_m_marker8, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, grouping.var = index_m_marker10)
    multiple_sample_clustering_markers_one_cluster$gene <- rownames(multiple_sample_clustering_markers_one_cluster)
	}
    else if (index_multiple_sample_normalization_method == "SCTransform"){
    multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
    multiple_sample_clustering_markers_one_cluster <- FindConservedMarkers(multiple_sample_clustering, ident.1 = index_m_marker8, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, grouping.var = index_m_marker10, assay = "SCT", recorrect_umi = FALSE) 
    multiple_sample_clustering_markers_one_cluster$gene <- rownames(multiple_sample_clustering_markers_one_cluster)
	} 
    return(list(data1 = multiple_sample_clustering_markers_one_cluster, data2 = multiple_sample_clustering))
    
  }
  else if(index_m_marker1 == 5){
    #Conserved Markers in two cluster
    if (index_multiple_sample_normalization_method == "LogNormalize"){
    multiple_sample_clustering_markers_two_cluster <- FindConservedMarkers(multiple_sample_clustering, ident.1 = index_m_marker8,  ident.2 = index_m_marker9, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, grouping.var = index_m_marker10)
    multiple_sample_clustering_markers_two_cluster$gene <- rownames(multiple_sample_clustering_markers_two_cluster)
	}
    else if (index_multiple_sample_normalization_method == "SCTransform"){
    multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
    multiple_sample_clustering_markers_two_cluster <- FindConservedMarkers(multiple_sample_clustering, ident.1 = index_m_marker8,  ident.2 = index_m_marker9, min.pct = index_m_marker2, logfc.threshold = index_m_marker3, test.use = index_m_marker4, only.pos = index_m_marker5, grouping.var = index_m_marker10, assay = "SCT", recorrect_umi = FALSE)
    multiple_sample_clustering_markers_two_cluster$gene <- rownames(multiple_sample_clustering_markers_two_cluster)
	} 
    return(list(data1 = multiple_sample_clustering_markers_two_cluster, data2 = multiple_sample_clustering))
  }
  
}
