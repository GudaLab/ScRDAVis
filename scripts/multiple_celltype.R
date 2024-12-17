datainput_multiple_celltype <- function(index_multiple_celltype_input, index_cell_markers, index_m_celltype1, index_m_celltype2, index_m_celltype3, index_m_celltype4, index_m_celltype5, index_m_celltype6, index_m_celltype7, index_m_celltype8, index_m_clustering6, index_multiple_sample_normalization_method){
  multiple_sample_clustering <- index_multiple_celltype_input
  multiple_sample_clustering_markers <- index_cell_markers
  if(index_m_celltype1 == 1){
    # load gene set preparation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # load cell type annotation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    # DB file
    db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
    tissue <- index_m_celltype2 # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus
    # prepare gene sets
    gs_list <- gene_sets_prepare(db_, tissue)
    
    if (index_multiple_sample_normalization_method == "LogNormalize"){
      # check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
      seurat_package_v5 <- isFALSE('counts' %in% names(attributes(multiple_sample_clustering[["RNA"]])));
      # extract scaled scRNA-seq matrix
      scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(multiple_sample_clustering[["RNA"]]$scale.data) else as.matrix(multiple_sample_clustering[["RNA"]]@scale.data)
    }
    else if (index_multiple_sample_normalization_method == "SCTransform"){
      # check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
      seurat_package_v5 <- isFALSE('counts' %in% names(attributes(multiple_sample_clustering[["SCT"]])));
      # extract scaled scRNA-seq matrix
      scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(multiple_sample_clustering[["SCT"]]$scale.data) else as.matrix(multiple_sample_clustering[["SCT"]]@scale.data)
    }
    # run ScType
    es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
    # merge by cluster
    cL_resutls <- do.call("rbind", lapply(unique(multiple_sample_clustering@meta.data$seurat_clusters), function(cl){
      es.max.cl = sort(rowSums(es.max[ ,rownames(multiple_sample_clustering@meta.data[multiple_sample_clustering@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(multiple_sample_clustering@meta.data$seurat_clusters==cl)), 10)
    }))
    sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    # set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
    #overlay the identified cell types on UMAP plot:
    multiple_sample_clustering@meta.data$sctype_classification = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,];
      multiple_sample_clustering@meta.data$sctype_classification[multiple_sample_clustering@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
    }
    
    plots40 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_celltype8, repel = TRUE, group.by = c('sctype_classification', 'seurat_clusters'))
    return(list(data1 = multiple_sample_clustering, data2 = unique(multiple_sample_clustering@meta.data$seurat_clusters), data3 = unique(multiple_sample_clustering@meta.data$sctype_classification), text_summary = 'sctype_classification', plot1 = plots40, data2 = sctype_scores))
  }
  
  else if(index_m_celltype1 == 2){
    ref <- fetchReference(index_m_celltype3, "2024-02-26")
    if (index_multiple_sample_normalization_method == "LogNormalize"){
      multiple_sample_clustering.counts<- GetAssayData(multiple_sample_clustering, layer="counts", assay = "RNA")
      #perform annotation
      Pred <- SingleR(as.SingleCellExperiment(multiple_sample_clustering, assay = "RNA"), ref = ref, labels = ref$label.main)
      
    }
    else if (index_multiple_sample_normalization_method == "SCTransform"){
      multiple_sample_clustering.counts<- GetAssayData(multiple_sample_clustering, layer="counts", assay = "SCT")
      #perform annotation
      Pred <- SingleR(as.SingleCellExperiment(multiple_sample_clustering, assay = "SCT"), ref = ref, labels = ref$label.main)
    }
    Pred <- SingleR(test = multiple_sample_clustering.counts, ref = ref, labels = ref$label.main,  assay.type.test=1, de.method = index_m_celltype4)
    table(Pred$labels)
    Pred2 <- table(Label=Pred$labels, Lost=is.na(Pred$pruned.labels)) %>% as.data.table
    Pred1 <- Pred %>% as.data.table
    Pred1$cells <- rownames(Pred)
    multiple_sample_clustering$singleR_labels <- Pred$labels[match(rownames(multiple_sample_clustering@meta.data), rownames(Pred))]
    
    plots40 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_celltype8, repel = TRUE, group.by=c('singleR_labels', 'seurat_clusters'))
    plots41 <- plotScoreHeatmap(Pred)
    plots42 <- plotDeltaDistribution(Pred)
    return(list(data1 = multiple_sample_clustering, data2 = unique(multiple_sample_clustering@meta.data$seurat_clusters), data3 = unique(multiple_sample_clustering@meta.data$singleR_labels), text_summary = 'singleR_labels', plot1 = plots40, data2 = Pred1, plot2 = plots41, plot3 = plots42))
  }
  
  else if(index_m_celltype1 == 3){
    #GPTCelltype
    # IMPORTANT! Assign your OpenAI API key. See Vignette for details
    
    res <- gptcelltype(multiple_sample_clustering_markers, model = index_m_celltype5, topgenenumber = index_m_celltype6)
    
    # Assign cell type annotation back to Seurat object
    multiple_sample_clustering@meta.data$GPTCelltype <- as.factor(res[as.character(Idents(multiple_sample_clustering))])
    
    # Visualize cell type annotation on UMAP
    plots40 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_celltype8, repel = TRUE,  group.by = c('GPTCelltype', 'seurat_clusters'))
    return(list(data1 = multiple_sample_clustering, data2 = unique(multiple_sample_clustering@meta.data$seurat_clusters), data3 = unique(multiple_sample_clustering@meta.data$GPTCelltype), text_summary = 'GPTCelltype', plot1 = plots40))
  }
  else if(index_m_celltype1 == 4){
    cell_type <- index_m_celltype7
    names(cell_type) <- levels(multiple_sample_clustering)
    multiple_sample_clustering@meta.data$cell_type <- as.factor(cell_type[as.character(Idents(multiple_sample_clustering))])
    
    plots40 <- DimPlot(multiple_sample_clustering, reduction = index_m_clustering6, label = index_m_celltype8, repel = TRUE,  group.by = c('cell_type','seurat_clusters'), order = as.numeric(sort(unique('seurat_clusters'))))
    return(list(data1 = multiple_sample_clustering, data2 = unique(multiple_sample_clustering@meta.data$seurat_clusters), data3 = unique(multiple_sample_clustering@meta.data$cell_type), text_summary = 'cell_type', plot1 = plots40))
    }
  
}
