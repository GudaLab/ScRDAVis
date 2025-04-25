datainput_multiple_doublet <- function(index_multiple_doublet_input, index_m_doublet1, index_m_doublet3, index_m_clustering6, index_multiple_sample_normalization_method, index_multiple_sample_pca_dim){
  multiple_sample_clustering <- index_multiple_doublet_input
  # #doubletfinder
   multiple_sample_clustering.split <- SplitObject(multiple_sample_clustering, split.by = "orig.ident") 
  # loop through samples to find doublets
   
   if (index_multiple_sample_normalization_method == "LogNormalize"){
    for (i in 1:length(multiple_sample_clustering.split)) {
    # print the sample we are on
    print(paste0("orig.ident",i))

    multiple_sample_clustering.sample <- multiple_sample_clustering.split[[i]]

    #min.pc = 30
    
    sweep.list <- paramSweep(multiple_sample_clustering.sample, PCs = 1:index_multiple_sample_pca_dim, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)


    ## Homotypic doublet proportion estimate
    annotations <- multiple_sample_clustering.sample@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations)
    nExp.poi <- round(index_m_doublet1 * nrow(multiple_sample_clustering.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))

    # run DoubletFinder
    
    multiple_sample_clustering.sample <- doubletFinder(seu = multiple_sample_clustering.sample,
                                                       PCs = 1:index_multiple_sample_pca_dim,
                                                       pK = 0.09,
                                                       nExp = nExp.poi.adj, sct = FALSE)
    metadata <- multiple_sample_clustering.sample@meta.data
    colnames(metadata)[9] <- "doublet_finder"
    multiple_sample_clustering.sample@meta.data <- metadata
    
    # subset and save
    multiple_sample_clustering.doublet <- subset(multiple_sample_clustering.sample)
    multiple_sample_clustering.split[[i]] <- multiple_sample_clustering.doublet
    remove(multiple_sample_clustering.doublet)
  }
}

   else if (index_multiple_sample_normalization_method == "SCTransform"){
     for (i in 1:length(multiple_sample_clustering.split)) {
       # print the sample we are on
       print(paste0("orig.ident",i))
       
       multiple_sample_clustering.sample <- multiple_sample_clustering.split[[i]]
       
       #min.pc = 30
       
       sweep.list <- paramSweep(multiple_sample_clustering.sample, PCs = 1:index_multiple_sample_pca_dim, sct = TRUE)
       sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
       bcmvn <- find.pK(sweep.stats)
       
       
       ## Homotypic doublet proportion estimate
       annotations <- multiple_sample_clustering.sample@meta.data$seurat_clusters
       homotypic.prop <- modelHomotypic(annotations)
       nExp.poi <- round(index_m_doublet1 * nrow(multiple_sample_clustering.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
       nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
       
       # run DoubletFinder
       
       multiple_sample_clustering.sample <- doubletFinder(seu = multiple_sample_clustering.sample,
                                                          PCs = 1:index_multiple_sample_pca_dim,
                                                          pK = 0.09,
                                                          nExp = nExp.poi.adj, sct = TRUE)
       metadata <- multiple_sample_clustering.sample@meta.data
       colnames(metadata)[11] <- "doublet_finder"
       multiple_sample_clustering.sample@meta.data <- metadata
       
       # subset and save
       multiple_sample_clustering.doublet <- subset(multiple_sample_clustering.sample)
       #multiple_sample_clustering.doublet <- PrepSCTFindMarkers(multiple_sample_clustering.doublet, assay = "SCT", verbose = TRUE)
       multiple_sample_clustering.split[[i]] <- multiple_sample_clustering.doublet
       remove(multiple_sample_clustering.doublet)
       multiple_sample_clustering.split[[i]]@assays$RNA <- NULL
       
     }
   }
   
  # convert multiple_sample_clustering.split
  multiple_sample_clustering.doublet <- merge(x = multiple_sample_clustering.split[[1]], y = multiple_sample_clustering.split[-1], add.cell.ids = names(multiple_sample_clustering.split), project="merged")

  #adding doublet to clustering file
  metadata <- multiple_sample_clustering.doublet@meta.data$doublet_finder
  #metadata <- multiple_sample_clustering.sample@meta.data
  multiple_sample_clustering@meta.data$doublet_finder <- metadata


  multiple_sample_clustering.singlets <- subset(multiple_sample_clustering, doublet_finder == "Singlet")
  multiple_sample_clustering.doublets <- subset(multiple_sample_clustering, doublet_finder == "Doublet")

  #counts
  table5 <- table(multiple_sample_clustering.singlets$orig.ident) %>% as.data.frame
  colnames(table5) <- c("Sample names", "Singlets")
  table6 <- table(multiple_sample_clustering.doublets$orig.ident)  %>% as.data.frame
  colnames(table6) <- c("Sample names", "Doublets")

  table(multiple_sample_clustering$doublet_finder)

  Doublet_count <- inner_join(table5, table6)

  plots22 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', raster=FALSE, label = index_m_doublet3)
  plots23 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', split.by= "condition", raster=FALSE, label = index_m_doublet3, ncol = 6)
  plots24 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', split.by= "orig.ident", raster=FALSE, label = index_m_doublet3, ncol = 6)
  plots25 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, group.by = 'doublet_finder', split.by= "seurat_clusters", raster=FALSE, label = index_m_doublet3, ncol = 5)

   if (index_multiple_sample_normalization_method == "SCTransform"){
   multiple_sample_clustering <- PrepSCTFindMarkers(multiple_sample_clustering, assay = "SCT", verbose = TRUE)
   }

  return(list(plot1 = plots22, plot2 = plots23, plot3 = plots24, plot4 = plots25, data1 = Doublet_count, data2 = multiple_sample_clustering))
}