datainput_subclustering_multiple_sample <- function(index_subclustering_multiple_sample_file, index_subclustering_multiple_sample_celltype, index_m_subclustering1, index_m_subclustering2, index_m_subclustering3){
  multiple_sample_clustering  <- index_subclustering_multiple_sample_file
  
  if (index_m_subclustering1 == "seurat_clusters"){
  Idents(multiple_sample_clustering) <- "seurat_clusters"
  subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering2)
  }
  else if (index_m_subclustering1 == "predicted"){
	if (index_subclustering_multiple_sample_celltype == "sctype_classification"){
    Idents(multiple_sample_clustering) <- "sctype_classification"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
   else if (index_subclustering_multiple_sample_celltype == "singleR_labels"){
    Idents(multiple_sample_clustering) <- "singleR_labels"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
   else if (index_subclustering_multiple_sample_celltype == "GPTCelltype"){
    Idents(multiple_sample_clustering) <- "GPTCelltype"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
   else if (index_subclustering_multiple_sample_celltype == "cell_type"){
    Idents(multiple_sample_clustering) <- "cell_type"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
  }
    
  table1 <- table(subclustering_multiple_sample$orig.ident) %>% as.data.frame 
  colnames(table1) <- c("Sample names", "Cell counts")
  
  subclustering_multiple_sample[["percent.mt"]] <- PercentageFeatureSet(subclustering_multiple_sample, pattern = "^MT-")


  plots1 <- VlnPlot(subclustering_multiple_sample, features = "nFeature_RNA", ncol = 1)
  plots2 <- VlnPlot(subclustering_multiple_sample, features = "nCount_RNA", ncol = 1)
  plots3 <- VlnPlot(subclustering_multiple_sample, features = "percent.mt", ncol = 1)
  
  return(list(plot = plots1+plots2+plots3, data1 = table1, data2 = subclustering_multiple_sample))
}