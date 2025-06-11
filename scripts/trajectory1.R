datainput_single_multiple_sample_trajectory1<- function(index_multiple_sample_input, index_subclustering_multiple_sample_input, index_multiple_sample_input2, index_subclustering_multiple_sample_input2, index_s_trajectory1, index_s_trajectory2, index_s_trajectory3, index_s_trajectory4, index_s_trajectory5, index_s_trajectory6, index_s_trajectory7, index_s_trajectory8){
  index_s_trajectory3 <- as.logical(index_s_trajectory3)
  index_s_trajectory4 <- as.logical(index_s_trajectory4)
  index_s_trajectory5 <- as.logical(index_s_trajectory5)
  index_s_trajectory6 <- as.logical(index_s_trajectory6)
  index_s_trajectory7 <- as.logical(index_s_trajectory7)
  index_s_trajectory8 <- as.logical(index_s_trajectory8)
  
  if (index_s_trajectory1 == "multiple_sample" & index_s_trajectory2 == "seurat_clusters"){
  single_multiple_sample_clustering <- index_multiple_sample_input 
  }
  else if (index_s_trajectory1 == "multiple_sample_subclustering" & index_s_trajectory2 == "seurat_clusters"){
  single_multiple_sample_clustering <- index_subclustering_multiple_sample_input  
  }
  else if (index_s_trajectory1 == "multiple_sample" & index_s_trajectory2 == "predicted"){
  single_multiple_sample_clustering <- index_multiple_sample_input
  Idents(single_multiple_sample_clustering) <- index_multiple_sample_input2
  }
  else if (index_s_trajectory1 == "multiple_sample_subclustering" & index_s_trajectory2 == "predicted"){
  single_multiple_sample_clustering <- index_subclustering_multiple_sample_input
  Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_input2
  }
  
  #Trajectory Analysis with Monocle3
  #Converting seuratobject to celldataset object for Monocle3
  cds <- as.cell_data_set(single_multiple_sample_clustering)
  fData(cds)$gene_short_name <- rownames(fData(cds))
  
  # data<- GetAssayData(single_multiple_sample_clustering,assay='RNA',layer = "data")
  # cell_metadata<-data.frame(cell_data=colnames(data))
  # gene_annotation<- data.frame(gene_short_name=rownames(data))
  # rownames(gene_annotation)<- rownames(data)
  # rownames(cell_metadata)<- colnames(data)
  # cds<- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)
  # fData(cds)$gene_short_name <- rownames(fData(cds))
  
  #Retrieve clustering information from Surat object
  # Assign partitions
  recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
  names(recreate.partitions) <- cds@colData@rownames
  recreate.partitions <- as.factor(recreate.partitions)
  recreate.partitions
  cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
  
  #Assign cluster information
  list.cluster <- single_multiple_sample_clustering@active.ident
  cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
    
  #Assign UMAP coordinates
  cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- single_multiple_sample_clustering@reductions$umap@cell.embeddings
  
  #Learn Trajectory
  cds <- learn_graph(cds, use_partition = index_s_trajectory3, close_loop = index_s_trajectory4)
  plots100 <- plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = index_s_trajectory5, label_branch_points = index_s_trajectory6, label_roots = index_s_trajectory7, label_leaves = index_s_trajectory8, group_label_size = 10, graph_label_size = 4)
  
  return(list(data1 = single_multiple_sample_clustering, data2 = cds, plot1 = plots100))
}