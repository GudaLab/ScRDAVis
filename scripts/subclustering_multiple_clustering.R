datainput_subclustering_multiple_clustering <- function(index_subclustering_multiple_clustering_input, index_subclustering_multiple_sample_normalization_method, index_m_subclustering_clustering1, index_m_subclustering_clustering2, index_m_subclustering_clustering3, index_m_subclustering_clustering4, index_m_subclustering_clustering5, index_m_subclustering_clustering6, index_m_subclustering_clustering7, index_m_subclustering_clustering8, index_m_subclustering_clustering9, index_m_subclustering_clustering10, index_m_subclustering_clustering11, index_m_subclustering_clustering12, index_m_subclustering_clustering13){
  index_m_subclustering_clustering5 <- as.numeric(index_m_subclustering_clustering5)
  label_size <- ifelse(as.logical(index_m_subclustering_clustering10) | as.logical(index_m_subclustering_clustering12), 3.5, 0)
  
  if (index_m_subclustering_clustering13 == "None"){
    subclustering_multiple_sample_clustering<- FindNeighbors(index_subclustering_multiple_clustering_input, dims = 1:index_m_subclustering_clustering1 , k.param = index_m_subclustering_clustering2, n.trees = index_m_subclustering_clustering3)
    subclustering_multiple_sample_clustering<- FindClusters(subclustering_multiple_sample_clustering, resolution = index_m_subclustering_clustering4, algorithm = index_m_subclustering_clustering5)
    if (index_m_subclustering_clustering6 == "umap")
    { 
      subclustering_multiple_sample_clustering<- RunUMAP(subclustering_multiple_sample_clustering, dims = 1:index_m_subclustering_clustering7, n.neighbors = index_m_subclustering_clustering8, min.dist = index_m_subclustering_clustering9)
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_m_subclustering_clustering10, raster=FALSE, group.by = "seurat_clusters")
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_m_subclustering_clustering10, raster=FALSE, group.by = "condition")
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_m_subclustering_clustering10, raster=FALSE, group.by = "orig.ident")
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_m_subclustering_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_m_subclustering_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    else if (index_m_subclustering_clustering6 == "tsne")
    {
      subclustering_multiple_sample_clustering<- RunTSNE(subclustering_multiple_sample_clustering, dims = 1:index_m_subclustering_clustering11)
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters")
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "condition")
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "orig.ident")
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    #cell_couts_in_custer
    subclustering_multiple_sample_clustering_cell_couts_in_custer <- table(subclustering_multiple_sample_clustering@meta.data$seurat_clusters) %>% as.data.table
    colnames(subclustering_multiple_sample_clustering_cell_couts_in_custer) <- c("Clusters", "Counts")
    
    plots19 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(seurat_clusters, fill = seurat_clusters)) +
      geom_bar(stat="count", position = position_dodge())+
      geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, position = position_dodge(0.9), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_condition
    subclustering_multiple_sample_clustering_cell_couts_in_condition <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("condition", "seurat_clusters")] %>% dcast(., condition ~ seurat_clusters, value.var = "N"))) %>%  rownames_to_column(var = "condition") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    
    plots20 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(seurat_clusters, fill = condition)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Condition"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_samples
    subclustering_multiple_sample_clustering_cell_couts_in_samples <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_samples[, .N, by = c("orig.ident", "seurat_clusters")] %>% dcast(., orig.ident ~ seurat_clusters, value.var = "N"))) %>%  rownames_to_column(var = "Clusters") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    colnames(subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples)[1] <- "Clusters"
    
    plots21 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(seurat_clusters, fill = orig.ident)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Clusters"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    

    #subclustering_multiple_sample_clustering <- JoinLayers(subclustering_multiple_sample_clustering)
    
    cluster_type = "seurat_clusters"
     }
  
  else if (index_subclustering_multiple_clustering13 == "HarmonyIntegration")
  {
    if (index_subclustering_multiple_sample_normalization_method == "LogNormalize"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = HarmonyIntegration, normalization.method = "LogNormalize", orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
    }
    else if (index_subclustering_multiple_sample_normalization_method == "SCTransform"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = HarmonyIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
    }
    #subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "harmony", verbose = FALSE)
    subclustering_multiple_sample_clustering<- FindNeighbors(subclustering_multiple_sample_clustering, reduction = "harmony", dims = 1:index_subclustering_multiple_clustering1 , k.param = index_subclustering_multiple_clustering2, n.trees = index_subclustering_multiple_clustering3)
    subclustering_multiple_sample_clustering<- FindClusters(subclustering_multiple_sample_clustering, resolution = index_subclustering_multiple_clustering4, algorithm = index_subclustering_multiple_clustering5, cluster.name = "harmony_clusters")
    if (index_subclustering_multiple_clustering6 == "umap")
    { 
      subclustering_multiple_sample_clustering<- RunUMAP(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering7, n.neighbors = index_subclustering_multiple_clustering8, min.dist = index_subclustering_multiple_clustering9, reduction = "harmony", reduction.name = "umap")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'harmony_clusters', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'condition', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'orig.ident', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    else if (index_subclustering_multiple_clustering6 == "tsne")
    {
      subclustering_multiple_sample_clustering<- RunTSNE(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering11, reduction = "harmony", reduction.name = "tsne")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'harmony_clusters', raster=FALSE, label = index_m_subclustering_clustering12)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'condition', raster=FALSE, label = index_m_subclustering_clustering12)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'orig.ident', raster=FALSE, label = index_m_subclustering_clustering12)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    #cell_couts_in_custer
    subclustering_multiple_sample_clustering_cell_couts_in_custer <- table(subclustering_multiple_sample_clustering@meta.data$harmony_clusters) %>% as.data.table
    colnames(subclustering_multiple_sample_clustering_cell_couts_in_custer) <- c("Clusters", "Counts")
    plots19 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(harmony_clusters, fill = harmony_clusters)) +
      geom_bar(stat="count", position = position_dodge())+
      geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, position = position_dodge(0.9), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_condition
    subclustering_multiple_sample_clustering_cell_couts_in_condition <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("condition", "harmony_clusters")] %>% dcast(., condition ~ harmony_clusters, value.var = "N"))) %>%  rownames_to_column(var = "condition") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    
    plots20 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(harmony_clusters, fill = condition)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Condition"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_samples
    subclustering_multiple_sample_clustering_cell_couts_in_samples <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("orig.ident", "harmony_clusters")] %>% dcast(., orig.ident ~ harmony_clusters, value.var = "N"))) %>%  rownames_to_column(var = "Clusters") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    colnames(subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples)[1] <- "Clusters"
    
    plots21 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(harmony_clusters, fill = orig.ident)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Clusters"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    
      #subclustering_multiple_sample_clustering <- JoinLayers(subclustering_multiple_sample_clustering)
    
    cluster_type = "harmony_clusters"
  }
  
  else if (index_subclustering_multiple_clustering13 == "CCAIntegration")
  {
    if (index_subclustering_multiple_sample_normalization_method == "LogNormalize"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = CCAIntegration, normalization.method = "LogNormalize", orig.reduction = "pca", new.reduction = "cca", verbose = FALSE)
    }
    else if (index_subclustering_multiple_sample_normalization_method == "SCTransform"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = CCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "cca", verbose = FALSE)
    }
    subclustering_multiple_sample_clustering<- FindNeighbors(subclustering_multiple_sample_clustering, reduction = "cca", dims = 1:index_subclustering_multiple_clustering1 , k.param = index_subclustering_multiple_clustering2, n.trees = index_subclustering_multiple_clustering3)
    subclustering_multiple_sample_clustering<- FindClusters(subclustering_multiple_sample_clustering, resolution = index_subclustering_multiple_clustering4, algorithm = index_subclustering_multiple_clustering5, cluster.name = "cca_clusters")
    if (index_subclustering_multiple_clustering6 == "umap")
    { 
      subclustering_multiple_sample_clustering<- RunUMAP(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering7, n.neighbors = index_subclustering_multiple_clustering8, min.dist = index_subclustering_multiple_clustering9, reduction = "cca", reduction.name = "umap")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'cca_clusters', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'condition', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'orig.ident', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    else if (index_subclustering_multiple_clustering6 == "tsne")
    {
      subclustering_multiple_sample_clustering<- RunTSNE(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering11, reduction = "cca", reduction.name = "tsne")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'cca_clusters', raster=FALSE, label = index_m_subclustering_clustering12)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'condition', raster=FALSE, label = index_m_subclustering_clustering12)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'orig.ident', raster=FALSE, label = index_m_subclustering_clustering12)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    #cell_couts_in_custer
    subclustering_multiple_sample_clustering_cell_couts_in_custer <- table(subclustering_multiple_sample_clustering@meta.data$cca_clusters) %>% as.data.table
    colnames(subclustering_multiple_sample_clustering_cell_couts_in_custer) <- c("Clusters", "Counts")
    plots19 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(cca_clusters, fill = cca_clusters)) +
      geom_bar(stat="count", position = position_dodge())+
      geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, position = position_dodge(0.9), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_condition
    subclustering_multiple_sample_clustering_cell_couts_in_condition <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("condition", "cca_clusters")] %>% dcast(., condition ~ cca_clusters, value.var = "N"))) %>%  rownames_to_column(var = "condition") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    
    plots20 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(cca_clusters, fill = condition)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Condition"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_samples
    subclustering_multiple_sample_clustering_cell_couts_in_samples <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("orig.ident", "cca_clusters")] %>% dcast(., orig.ident ~ cca_clusters, value.var = "N"))) %>%  rownames_to_column(var = "Clusters") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    colnames(subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples)[1] <- "Clusters"
    
    plots21 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(cca_clusters, fill = orig.ident)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Clusters"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #subclustering_multiple_sample_clustering <- JoinLayers(subclustering_multiple_sample_clustering)    
    cluster_type = "cca_clusters"
  }
  
  else if (index_subclustering_multiple_clustering13 == "RPCAIntegration")
  {
    if (index_subclustering_multiple_sample_normalization_method == "LogNormalize"){
    subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = RPCAIntegration, normalization.method = "LogNormalize", orig.reduction = "pca", new.reduction = "rpca", verbose = FALSE)
    }
    else if (index_subclustering_multiple_sample_normalization_method == "SCTransform"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = RPCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "rpca", verbose = FALSE)
    }
    subclustering_multiple_sample_clustering<- FindNeighbors(subclustering_multiple_sample_clustering, reduction = "rpca", dims = 1:index_subclustering_multiple_clustering1 , k.param = index_subclustering_multiple_clustering2, n.trees = index_subclustering_multiple_clustering3)
    subclustering_multiple_sample_clustering<- FindClusters(subclustering_multiple_sample_clustering, resolution = index_subclustering_multiple_clustering4, algorithm = index_subclustering_multiple_clustering5, cluster.name = "rpca_clusters")
    if (index_subclustering_multiple_clustering6 == "umap")
    { 
      subclustering_multiple_sample_clustering<- RunUMAP(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering7, n.neighbors = index_subclustering_multiple_clustering8, min.dist = index_subclustering_multiple_clustering9, reduction = "rpca", reduction.name = "umap")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'rpca_clusters', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'condition', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'orig.ident', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    else if (index_subclustering_multiple_clustering6 == "tsne")
    {
      subclustering_multiple_sample_clustering<- RunTSNE(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering11, reduction = "rpca", reduction.name = "tsne")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'rpca_clusters', raster=FALSE, label = index_m_subclustering_clustering12)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'condition', raster=FALSE, label = index_m_subclustering_clustering12)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'orig.ident', raster=FALSE, label = index_m_subclustering_clustering12)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    #cell_couts_in_custer
    subclustering_multiple_sample_clustering_cell_couts_in_custer <- table(subclustering_multiple_sample_clustering@meta.data$rpca_clusters) %>% as.data.table
    colnames(subclustering_multiple_sample_clustering_cell_couts_in_custer) <- c("Clusters", "Counts")
    plots19 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(rpca_clusters, fill = rpca_clusters)) +
      geom_bar(stat="count", position = position_dodge())+
      geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, position = position_dodge(0.9), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_condition
    subclustering_multiple_sample_clustering_cell_couts_in_condition <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("condition", "rpca_clusters")] %>% dcast(., condition ~ rpca_clusters, value.var = "N"))) %>%  rownames_to_column(var = "condition") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    
    plots20 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(rpca_clusters, fill = condition)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Condition"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_samples
    subclustering_multiple_sample_clustering_cell_couts_in_samples <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("orig.ident", "rpca_clusters")] %>% dcast(., orig.ident ~ rpca_clusters, value.var = "N"))) %>%  rownames_to_column(var = "Clusters") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    colnames(subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples)[1] <- "Clusters"
    
    plots21 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(rpca_clusters, fill = orig.ident)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Clusters"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #subclustering_multiple_sample_clustering <- JoinLayers(subclustering_multiple_sample_clustering)    
    cluster_type = "rpca_clusters"
  }
  
  else if (index_subclustering_multiple_clustering13 == "JointPCAIntegration")
  {
    if (index_subclustering_multiple_sample_normalization_method == "LogNormalize"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = JointPCAIntegration, normalization.method = "LogNormalize", orig.reduction = "pca", new.reduction = "jointpca", verbose = FALSE)
    }
    else if (index_subclustering_multiple_sample_normalization_method == "SCTransform"){
      subclustering_multiple_sample_clustering<- IntegrateLayers(object = index_subclustering_multiple_clustering_input, method = JointPCAIntegration, normalization.method = "SCT", orig.reduction = "pca", new.reduction = "jointpca", verbose = FALSE)
    }
    subclustering_multiple_sample_clustering<- FindNeighbors(subclustering_multiple_sample_clustering, reduction = "jointpca", dims = 1:index_subclustering_multiple_clustering1 , k.param = index_subclustering_multiple_clustering2, n.trees = index_subclustering_multiple_clustering3)
    subclustering_multiple_sample_clustering<- FindClusters(subclustering_multiple_sample_clustering, resolution = index_subclustering_multiple_clustering4, algorithm = index_subclustering_multiple_clustering5, cluster.name = "jointpca_clusters")
    if (index_subclustering_multiple_clustering6 == "umap")
    { 
      subclustering_multiple_sample_clustering<- RunUMAP(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering7, n.neighbors = index_subclustering_multiple_clustering8, min.dist = index_subclustering_multiple_clustering9, reduction = "jointpca", reduction.name = "umap")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'jointpca_clusters', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'condition', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'umap', group.by = 'orig.ident', raster=FALSE, label = index_subclustering_multiple_clustering10)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "umap", label = index_subclustering_multiple_clustering10, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    else if (index_subclustering_multiple_clustering6 == "tsne")
    {
      subclustering_multiple_sample_clustering<- RunTSNE(subclustering_multiple_sample_clustering, dims = 1:index_subclustering_multiple_clustering11, reduction = "jointpca", reduction.name = "tsne")
      plots16 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'jointpca_clusters', raster=FALSE, label = index_m_subclustering_clustering12)
      plots17 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'condition', raster=FALSE, label = index_m_subclustering_clustering12)
      plots18 <-DimPlot(subclustering_multiple_sample_clustering, reduction = 'tsne', group.by = 'orig.ident', raster=FALSE, label = index_m_subclustering_clustering12)
      plots22 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "condition", ncol = 6)
      plots23 <-DimPlot(subclustering_multiple_sample_clustering, reduction = "tsne", label = index_m_subclustering_clustering12, raster=FALSE, group.by = "seurat_clusters", split.by= "orig.ident", ncol = 6)
    }
    #cell_couts_in_custer
    subclustering_multiple_sample_clustering_cell_couts_in_custer <- table(subclustering_multiple_sample_clustering@meta.data$jointpca_clusters) %>% as.data.table
    colnames(subclustering_multiple_sample_clustering_cell_couts_in_custer) <- c("Clusters", "Counts")
    plots19 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(jointpca_clusters, fill = jointpca_clusters)) +
      geom_bar(stat="count", position = position_dodge())+
      geom_text(stat='count', aes(label=after_stat(count)), vjust=-0.5, position = position_dodge(0.9), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_condition
    subclustering_multiple_sample_clustering_cell_couts_in_condition <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("condition", "jointpca_clusters")] %>% dcast(., condition ~ jointpca_clusters, value.var = "N"))) %>%  rownames_to_column(var = "condition") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    
    plots20 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(jointpca_clusters, fill = condition)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Condition"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #cell_couts_in_custer_for_each_samples
    subclustering_multiple_sample_clustering_cell_couts_in_samples <- subclustering_multiple_sample_clustering@meta.data %>% as.data.table
    subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples <- data.frame(t(subclustering_multiple_sample_clustering_cell_couts_in_condition[, .N, by = c("orig.ident", "jointpca_clusters")] %>% dcast(., orig.ident ~ jointpca_clusters, value.var = "N"))) %>%  rownames_to_column(var = "Clusters") %>% `colnames<-`(.[1, ]) %>%  .[-1, ]
    colnames(subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples)[1] <- "Clusters"
    
    plots21 <- ggplot(subclustering_multiple_sample_clustering@meta.data, aes(jointpca_clusters, fill = orig.ident)) +
      geom_bar(stat="count")+
      geom_text(stat='count', aes(label=after_stat(count)), position = position_stack(vjust = 0.5), size=label_size)+
      theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
      theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Clusters"))+
      theme(axis.text.x = element_text(angle = 90, vjust = 1))
    
    #subclustering_multiple_sample_clustering <- JoinLayers(subclustering_multiple_sample_clustering)
        cluster_type = "jointpca_clusters"
  }
  if (index_subclustering_multiple_sample_normalization_method == "LogNormalize"){
    subclustering_multiple_sample_clustering <- JoinLayers(subclustering_multiple_sample_clustering)
  }
  else if (index_subclustering_multiple_sample_normalization_method == "SCTransform"){
  subclustering_multiple_sample_clustering <- PrepSCTFindMarkers(subclustering_multiple_sample_clustering, assay = "SCT", verbose = TRUE)
  }
  return(list(plot1 = plots16, plot2 = plots19, plot3 = plots17, plot4 = plots20, plot5 = plots18, plot6 = plots21, data1 = subclustering_multiple_sample_clustering_cell_couts_in_custer, data2 = subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_condition, data3 = subclustering_multiple_sample_clustering_total_cell_couts_in_custer_for_each_samples, data4 = subclustering_multiple_sample_clustering, data5 = unique(subclustering_multiple_sample_clustering@meta.data$seurat_clusters), data6 = max(as.numeric(subclustering_multiple_sample_clustering@meta.data$seurat_clusters)), data7 = unique(subclustering_multiple_sample_clustering@meta.data$condition), text_summary=cluster_type, plot7=plots22+plots23))
}