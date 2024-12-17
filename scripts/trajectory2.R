datainput_single_multiple_sample_trajectory2 <- function(index_trajectory2_input1, index_trajectory2_input2,  index_trajectory2_multiple_sample_input2, index_trajectory2_subclustering_multiple_sample_input2,  index_s_trajectory1, index_s_trajectory2, index_s_trajectory10, index_s_trajectory11, index_s_trajectory12, index_s_trajectory13, index_s_trajectory14){
  index_s_trajectory11 <- as.logical(index_s_trajectory11)
  index_s_trajectory12 <- as.logical(index_s_trajectory12)
  index_s_trajectory13 <- as.logical(index_s_trajectory13)
  index_s_trajectory14 <- as.logical(index_s_trajectory14)
  single_multiple_sample_clustering <- index_trajectory2_input1
  cds <- index_trajectory2_input2
  
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[, monocle3::clusters(cds) == index_s_trajectory10]))
  plots101 <- plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = index_s_trajectory11, label_branch_points = index_s_trajectory12, label_roots = index_s_trajectory13, label_leaves = index_s_trajectory14, graph_label_size = 5, group_label_size = 5)
  #Cells ordered by Monocle3 Pseudotime
  cds$monocle3_pseudotime <- pseudotime(cds)
  table_pseudotime <- as.data.frame(colData(cds))
  
  if (index_s_trajectory1 == "multiple_sample" & index_s_trajectory2 == "seurat_clusters"){
    plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
    plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()
    plots101+plots103
  }
  else if (index_s_trajectory1 == "multiple_sample_subclustering" & index_s_trajectory2 == "seurat_clusters"){
    plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, seurat_clusters, fill = seurat_clusters)) + geom_boxplot()
    plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(seurat_clusters, monocle3_pseudotime), fill = seurat_clusters)) + geom_boxplot()
    plots101+plots103
  }
  
  else if (index_s_trajectory1 == "multiple_sample" & index_s_trajectory2 == "predicted"){
    if (index_trajectory2_multiple_sample_input2 == "sctype_classification"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, sctype_classification, fill = sctype_classification)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(sctype_classification, monocle3_pseudotime), fill = sctype_classification)) + geom_boxplot()
      plots101+plots103
    }
    else if (index_trajectory2_multiple_sample_input2 == "singleR_labels"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, singleR_labels, fill = singleR_labels)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(singleR_labels, monocle3_pseudotime), fill = singleR_labels)) + geom_boxplot()
      plots101+plots103
    }
    else if (index_trajectory2_multiple_sample_input2 == "GPTCelltype"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, GPTCelltype, fill = GPTCelltype)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(GPTCelltype, monocle3_pseudotime), fill = GPTCelltype)) + geom_boxplot()
      plots101+plots103
    }
    else if (index_trajectory2_multiple_sample_input2 == "cell_type"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, cell_type, fill = cell_type)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(cell_type, monocle3_pseudotime), fill = cell_type)) + geom_boxplot()
      plots101+plots103
    }
  }
  else if (index_s_trajectory1 == "multiple_sample_subclustering" & index_s_trajectory2 == "predicted"){
    if (index_trajectory2_subclustering_multiple_sample_input2 == "sctype_classification"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, sctype_classification, fill = sctype_classification)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(sctype_classification, monocle3_pseudotime), fill = sctype_classification)) + geom_boxplot()
      plots101+plots103
    }
    else if (index_trajectory2_subclustering_multiple_sample_input2 == "singleR_labels"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, singleR_labels, fill = singleR_labels)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(singleR_labels, monocle3_pseudotime), fill = singleR_labels)) + geom_boxplot()
      plots101+plots103
    }
    else if (index_trajectory2_subclustering_multiple_sample_input2 == "GPTCelltype"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, GPTCelltype, fill = GPTCelltype)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(GPTCelltype, monocle3_pseudotime), fill = GPTCelltype)) + geom_boxplot()
      plots101+plots103
    }
    else if (index_trajectory2_subclustering_multiple_sample_input2 == "cell_type"){
      plots102 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, cell_type, fill = cell_type)) + geom_boxplot()
      plots103 <- ggplot(table_pseudotime, aes(monocle3_pseudotime, reorder(cell_type, monocle3_pseudotime), fill = cell_type)) + geom_boxplot()
      plots101+plots103
    }
  }
  
  return(list(data1 = single_multiple_sample_clustering, data2 = cds, plot1 = plots101, plot2 = plots102+plots103))
}
