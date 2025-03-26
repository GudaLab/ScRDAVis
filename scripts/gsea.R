datainput_single_multiple_sample_gsea<- function(index_multiple_sample_gsea_input, index_subclustering_multiple_sample_gsea_input, index_multiple_sample_gsea_input2, index_subclustering_multiple_sample_gsea_input2, index_multiple_sample_gsea_input3, index_subclustering_multiple_sample_gsea_input3, index_s_gsea1, index_s_gsea2, index_s_gsea3, index_s_gsea4, index_s_gsea5, index_s_gsea6, index_s_gsea7, index_s_gsea8, index_s_gsea9, index_s_gsea10, index_s_gsea11, index_s_gsea12){
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
  
  gene_ranking <- sort(gene_ranking, decreasing = TRUE)
  
  # Get the gene sets from MSigDB (e.g., Hallmark gene sets)
  msigdb_gene_sets <- msigdbr(species = index_s_gsea5, collection = index_s_gsea6)
  gene_sets <- split(msigdb_gene_sets$gene_symbol, msigdb_gene_sets$gs_name)
  
  fgsea_results <- fgsea(pathways = gene_sets, stats = gene_ranking, scoreType = index_s_gsea7, minSize = index_s_gsea8, maxSize = index_s_gsea9, nPermSimple = index_s_gsea10)
  #fgsea_results  <- subset(fgsea_results, (fgsea_results$padj <=0.05))
  
  # View significant pathways
  topPathwaysup <- fgsea_results[ES > 0] %>% arrange(padj) %>% head(n = index_s_gsea12/2)
  topPathwaysdown <- fgsea_results[ES < 0] %>% arrange(padj) %>% head(n = index_s_gsea12/2)
  topPathways <- rbind(topPathwaysup, topPathwaysdown)
  topPathways1Up <- fgsea_results[ES > 0][head(order(padj), n=index_s_gsea12/2), pathway]
  topPathways1Down <- fgsea_results[ES < 0][head(order(padj), n=index_s_gsea12/2), pathway]
  topPathways1 <- c(topPathways1Up, rev(topPathways1Down))
  
  # Visualize GSEA results
   if(index_s_gsea11 == "GSEA_plot"){
    plot_list <- list()
    for (i in 1:index_s_gsea12) {
      plot<-plotEnrichment(gene_sets[[topPathways1[i]]], gene_ranking)+ labs(title=topPathways1[i])
      plot_list[[i]] <- plot
    }
    plots501<-grid.arrange(grobs = plot_list[1:index_s_gsea12])
    
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
  
  return(list(plot1 = plots501, data1 = fgsea_results))
}
