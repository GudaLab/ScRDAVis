datainput_multiple_doublet2 <- function(index_multiple_doublet1_input, index_m_doublet2, index_m_doublet3, index_m_clustering6, index_multiple_sample_normalization_method){
  
  multiple_sample_clustering <- index_multiple_doublet1_input 
  
  if (index_m_doublet2 == "remove_doublet"){
    #remove doublet
    multiple_sample_clustering <- subset(multiple_sample_clustering, subset = doublet_finder == "Singlet")
    table7 <- table(multiple_sample_clustering$orig.ident) %>% as.data.frame
    colnames(table7) <- c("Sample names", "No of Singlets used for further analysis")
  }
  
  else if (index_m_doublet2 == "keep_doublet"){
    multiple_sample_clustering <- index_multiple_doublet1_input
    table7 <- table(multiple_sample_clustering$orig.ident) %>% as.data.frame
    colnames(table7) <- c("Sample names", "No of Singlets / Doublet used for further analysis")
  }
  
  plots26 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', raster=FALSE, label = index_m_doublet3)
  plots27 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', split.by= "condition", raster=FALSE, label = index_m_doublet3, ncol = 6)
  plots28 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', split.by= "orig.ident", raster=FALSE, label = index_m_doublet3, ncol = 6)
  plots29 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', split.by= "seurat_clusters", raster=FALSE, label = index_m_doublet3, ncol = 5)
  plots30 <-DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_doublet3, group.by = "seurat_clusters")
  plots31 <-DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_doublet3, group.by = "condition")
  plots32 <-DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_doublet3, group.by = "orig.ident")
  plots33 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'condition', split.by= "seurat_clusters", raster=FALSE, label = index_m_doublet3, ncol = 5)
  
  multiple_sample_clustering_cell_couts_in_custer1 <- table(multiple_sample_clustering@meta.data$seurat_clusters) %>% as.data.table
  colnames(multiple_sample_clustering_cell_couts_in_custer1) <- c("Clusters", "Counts")
  
  plots34 <- ggplot(multiple_sample_clustering@meta.data, aes(seurat_clusters, fill = seurat_clusters)) +
    geom_bar(stat="count", position = position_dodge())+
    geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, position = position_dodge(0.9), size=3.5)+
    theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
    theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  #cell_couts_in_custer_for_each_condition
  multiple_sample_clustering_cell_couts_in_condition1 <- multiple_sample_clustering@meta.data %>% as.data.table
  multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition1 <- data.frame(t(multiple_sample_clustering_cell_couts_in_condition1[, .N, by = c("condition", "seurat_clusters")] %>% dcast(., condition ~ seurat_clusters, value.var = "N"))) %>%  rownames_to_column(var = "condition") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
  
  plots35 <- ggplot(multiple_sample_clustering@meta.data, aes(seurat_clusters, fill = condition)) +
    geom_bar(stat="count")+
    geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=3.5)+
    theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
    theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Condition"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  #cell_couts_in_custer_for_each_samples
  multiple_sample_clustering_cell_couts_in_samples1 <- multiple_sample_clustering@meta.data %>% as.data.table
  multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples1 <- data.frame(t(multiple_sample_clustering_cell_couts_in_condition1[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N"))) %>%  rownames_to_column(var = "Clusters") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
  colnames(multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples1)[1] <- "Clusters"
  
  plots36 <- ggplot(multiple_sample_clustering@meta.data, aes(seurat_clusters, fill = orig.ident)) +
    geom_bar(stat="count")+
    geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=3.5)+
    theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
    theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Clusters"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  plots37 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'seurat_clusters', split.by= "condition", raster=FALSE, label = index_m_doublet3, ncol = 6)
  plots38 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'seurat_clusters', split.by= "orig.ident", raster=FALSE, label = index_m_doublet3, ncol = 6)
  
   if (index_multiple_sample_normalization_method == "SCTransform"){
     multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
   }
  
  return(list(plot1 = plots26+plots30, plot2 = plots27+plots31, plot3 = plots28+plots32, plot4 = plots29+plots33, plot5=plots34+plots35+plots36, data1= table7, data2 = multiple_sample_clustering_cell_couts_in_custer1, data3 = multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition1, data4 = multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples1, data5 = multiple_sample_clustering, plot6=plots37+plots38))
}