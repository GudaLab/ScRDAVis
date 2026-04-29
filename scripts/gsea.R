datainput_single_multiple_sample_gsea<- function(index_multiple_sample_gsea_input, index_subclustering_multiple_sample_gsea_input, index_multiple_sample_gsea_input2, index_subclustering_multiple_sample_gsea_input2, index_multiple_sample_gsea_input3, index_subclustering_multiple_sample_gsea_input3, index_s_gsea1, index_s_gsea2, index_s_gsea3, index_s_gsea4, index_s_gsea5, index_s_gsea6, index_s_gsea7, index_s_gsea8, index_s_gsea9, index_s_gsea10, index_s_gsea11, index_s_gsea12, progress_callback = NULL){
  notify_progress <- function(value = NULL, detail = NULL) {
    if (is.function(progress_callback)) {
      progress_callback(value = value, detail = detail)
    }
  }
  
  notify_progress(0.1, "Preparing the ranked gene list.")
  
    if (index_s_gsea1 == "multiple_sample" & index_s_gsea2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_multiple_sample_gsea_input
    single_sample_clustering_markers <- index_multiple_sample_gsea_input3 
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
    gene_ranking <- deg_genes$avg_log2FC
    names(gene_ranking) <- deg_genes$gene
  }
  else if (index_s_gsea1 == "multiple_sample_subclustering" & index_s_gsea2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_subclustering_multiple_sample_gsea_input 
    single_sample_clustering_markers <- index_subclustering_multiple_sample_gsea_input3 
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
    gene_ranking <- deg_genes$avg_log2FC
    names(gene_ranking) <- deg_genes$gene
  }
    
  else if (index_s_gsea1 == "multiple_sample" & index_s_gsea2 == "predicted"){
    if (index_multiple_sample_gsea_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_gsea3  & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
    else if (index_multiple_sample_gsea_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
    else if (index_multiple_sample_gsea_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
    else if (index_multiple_sample_gsea_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
  }
  
  else if (index_s_gsea1 == "multiple_sample_subclustering" & index_s_gsea2 == "predicted"){
    if (index_subclustering_multiple_sample_gsea_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_gsea_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_gsea_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_gsea_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_gsea_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_gsea_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_gsea_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_gsea3 & p_val_adj < index_s_gsea4))
      gene_ranking <- deg_genes$avg_log2FC
      names(gene_ranking) <- deg_genes$gene
    }
  }
  notify_progress(0.25, "Cleaning ranked genes for GSEA.")
  gene_ranking <- sort(gene_ranking, decreasing = TRUE)
  gene_ranking <- gene_ranking[!duplicated(names(gene_ranking))]
  
  if (length(gene_ranking) == 0) {
    stop("No ranked genes were available for GSEA. Please adjust the selected cluster(s) or p_val_adj cutoff.")
  }
  
  notify_progress(0.4, "Loading the selected MSigDB collection.")
  msigdb_db_species <- if (identical(index_s_gsea5, "Mus musculus")) "MM" else "HS"
  msigdb_gene_sets <- msigdbr(
    db_species = msigdb_db_species,
    species = index_s_gsea5,
    collection = index_s_gsea6
  )
  gene_sets <- split(msigdb_gene_sets$gene_symbol, msigdb_gene_sets$gs_name)
  
  if (length(gene_sets) == 0) {
    stop("No MSigDB gene sets were returned for the selected organism and collection.")
  }
  
  notify_progress(0.65, "Running fgsea on the selected pathways.")
  fgsea_results <- fgsea(pathways = gene_sets, stats = gene_ranking, scoreType = index_s_gsea7, minSize = index_s_gsea8, maxSize = index_s_gsea9, nPermSimple = index_s_gsea10)
  #fgsea_results  <- subset(fgsea_results, (fgsea_results$padj <=0.05))
  fgsea_results <- as.data.frame(fgsea_results)
  
  if (nrow(fgsea_results) == 0) {
    stop("fgsea returned no results. Please try a different collection or adjust the gene set size settings.")
  }
  
  # View significant pathways
  pathway_limit <- max(1L, floor(as.integer(index_s_gsea12) / 2))
  topPathwaysup <- fgsea_results %>% dplyr::filter(ES > 0) %>% dplyr::arrange(padj) %>% head(n = pathway_limit)
  topPathwaysdown <- fgsea_results %>% dplyr::filter(ES < 0) %>% dplyr::arrange(padj) %>% head(n = pathway_limit)
  topPathways <- rbind(topPathwaysup, topPathwaysdown)
  topPathways1Up <- topPathwaysup$pathway
  topPathways1Down <- topPathwaysdown$pathway
  topPathways1 <- unique(c(topPathways1Up, rev(topPathways1Down)))
  topPathways1 <- topPathways1[!is.na(topPathways1)]
  
  if (length(topPathways1) == 0) {
    stop("No enriched pathways were available to plot for the selected inputs.")
  }
  
  plot_count <- min(length(topPathways1), max(1L, as.integer(index_s_gsea12)))
  
  # Visualize GSEA results
  notify_progress(0.85, "Generating the selected GSEA plots.")
   if(index_s_gsea11 == "GSEA_plot"){
    plot_list <- list()
    for (i in seq_len(plot_count)) {
      plot<-plotEnrichment(gene_sets[[topPathways1[i]]], gene_ranking)+ labs(title=topPathways1[i])
      plot_list[[i]] <- plot
    }
    plots501<-grid.arrange(grobs = plot_list[seq_len(plot_count)])
    
  }
  else if (index_s_gsea11 == "barplot"){
    topPathways$Significant <- ifelse(topPathways$ES > 0, "Positive", "Negative")
    plots501 <- ggplot(topPathways, aes(reorder(pathway, ES), ES, fill = Significant)) +
      geom_col() + coord_flip() +
      labs(title = "Top Pathways Enriched", x = "Pathway", y = "Normalized Enrichment Score")+theme_bw()
    }
  else if (index_s_gsea11 == "plotGseaTable"){
    plots501 <-plotGseaTable(gene_sets[topPathways1], gene_ranking, fgsea_results, gseaParam=0.5)
  }
    fgsea_results$leadingEdge <- sapply(fgsea_results$leadingEdge, function(x) paste(unlist(x), collapse = ", "))
  notify_progress(1, "GSEA results are ready.")
	
  return(list(plot1 = plots501, data1 = fgsea_results))
}
