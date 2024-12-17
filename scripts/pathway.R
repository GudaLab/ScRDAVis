datainput_single_multiple_sample_pathway<- function(index_multiple_sample_pathway_input, index_subclustering_multiple_sample_pathway_input, index_multiple_sample_pathway_input2, index_subclustering_multiple_sample_pathway_input2, index_multiple_sample_pathway_input3, index_subclustering_multiple_sample_pathway_input3, index_s_pathway1, index_s_pathway2, index_s_pathway3, index_s_pathway4, index_s_pathway5, index_s_pathway6, index_s_pathway7, index_s_pathway8, index_s_pathway9, index_s_pathway10, index_s_pathway11, index_s_pathway12, index_s_pathway13, index_s_pathway14){
  if (index_s_pathway1 == "gene_name_list"){
    deg_genes <- unlist(strsplit(index_s_pathway14, ","))
  }
  else if (index_s_pathway1 == "multiple_sample" & index_s_pathway2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_multiple_sample_pathway_input
    single_sample_clustering_markers <- index_multiple_sample_pathway_input3 
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
    deg_genes <- deg_genes$gene
  }
  else if (index_s_pathway1 == "multiple_sample_subclustering" & index_s_pathway2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input 
    single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3 
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
    deg_genes <- deg_genes$gene
  }
  else if (index_s_pathway1 == "multiple_sample" & index_s_pathway2 == "predicted"){
    if (index_multiple_sample_pathway_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_pathway3  & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (index_multiple_sample_pathway_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (index_multiple_sample_pathway_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (index_multiple_sample_pathway_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
  }
  
  else if (index_s_pathway1 == "multiple_sample_subclustering" & index_s_pathway2 == "predicted"){
    if (index_subclustering_multiple_sample_pathway_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_pathway_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_pathway_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_pathway_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
  }
  
  entrez_ids <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = index_s_pathway5)
  entrez_gene_list <- entrez_ids$ENTREZID
  
  if (index_s_pathway6 == "KEGG"){
    if(index_s_pathway5 == "org.Hs.eg.db"){
    pathway_enrichment <- enrichKEGG(gene = entrez_gene_list, organism = 'hsa', pAdjustMethod = index_s_pathway7, pvalueCutoff = index_s_pathway8, qvalueCutoff = index_s_pathway9, minGSSize = index_s_pathway10, maxGSSize = index_s_pathway11)
    pathway_enrichment@result$geneID <- sapply(strsplit(pathway_enrichment@result$geneID, "/"), function(genes) {
    gene_symbols <- entrez_ids$SYMBOL[match(genes, entrez_ids$ENTREZID)]
    paste(gene_symbols, collapse = "/")
    }) 
	}
    else if(index_s_pathway5 == "org.Mm.eg.db"){
    pathway_enrichment <- enrichKEGG(gene = entrez_gene_list,organism = 'mmu', pAdjustMethod = index_s_pathway7, pvalueCutoff = index_s_pathway8, qvalueCutoff = index_s_pathway9, minGSSize = index_s_pathway10, maxGSSize = index_s_pathway11)
    pathway_enrichment@result$geneID <- sapply(strsplit(pathway_enrichment@result$geneID, "/"), function(genes) {
    gene_symbols <- entrez_ids$SYMBOL[match(genes, entrez_ids$ENTREZID)]
    paste(gene_symbols, collapse = "/")
    }) 
	}
    else if(index_s_pathway5 == "org.Rn.eg.db"){
    pathway_enrichment <- enrichKEGG(gene = entrez_gene_list, organism = 'rno', pAdjustMethod = index_s_pathway7, pvalueCutoff = index_s_pathway8, qvalueCutoff = index_s_pathway9, minGSSize = index_s_pathway10, maxGSSize = index_s_pathway11)
    pathway_enrichment@result$geneID <- sapply(strsplit(pathway_enrichment@result$geneID, "/"), function(genes) {
    gene_symbols <- entrez_ids$SYMBOL[match(genes, entrez_ids$ENTREZID)]
    paste(gene_symbols, collapse = "/")
    }) 
	}
    pathway_results <- as.data.frame(pathway_enrichment@result)
    #pathway_results$geneID_with_symbols <- sapply(pathway_results$geneID, function(x) {
    #  gene_ids <- unlist(strsplit(x, "/"))
    #  symbols <- entrez_ids$SYMBOL[match(gene_ids, entrez_ids$ENTREZID)]
    #  paste(symbols, collapse = "/")
    #})
    
  pathway_results <- pathway_results %>% dplyr::select(Description, geneID, everything())
    
  }
  else if (index_s_pathway6 == "Reactome"){
    if(index_s_pathway5 == "org.Hs.eg.db"){
      pathway_enrichment <- enrichPathway(gene = entrez_gene_list, organism = 'human', pAdjustMethod = index_s_pathway7, pvalueCutoff = index_s_pathway8, qvalueCutoff = index_s_pathway9, minGSSize = index_s_pathway10, maxGSSize = index_s_pathway11, readable = TRUE)
    }
    else if(index_s_pathway5 == "org.Mm.eg.db"){
      pathway_enrichment <- enrichPathway(gene = entrez_gene_list,organism = 'mouse', pAdjustMethod = index_s_pathway7, pvalueCutoff = index_s_pathway8, qvalueCutoff = index_s_pathway9, minGSSize = index_s_pathway10, maxGSSize = index_s_pathway11, readable = TRUE)
    }
    else if(index_s_pathway5 == "org.Rn.eg.db"){
      pathway_enrichment <- enrichPathway(gene = entrez_gene_list, organism = 'rat', pAdjustMethod = index_s_pathway7, pvalueCutoff = index_s_pathway8, qvalueCutoff = index_s_pathway9, minGSSize = index_s_pathway10, maxGSSize = index_s_pathway11, readable = TRUE)
    } 
    pathway_results <- pathway_enrichment@result
  }
    
   
  if(index_s_pathway12 == "dotplot"){
    plots301 <- dotplot(pathway_enrichment, showCategory = index_s_pathway13)
  }
  else if (index_s_pathway12 == "barplot"){
    plots301 <- barplot(pathway_enrichment, showCategory = index_s_pathway13)
  }
  else if (index_s_pathway12 == "goplot"){
    plots301 <- goplot(pathway_enrichment, showCategory = index_s_pathway13)
  }
  else if (index_s_pathway12 == "cnetplot"){
    plots301 <- cnetplot(pathway_enrichment, showCategory = index_s_pathway13)
  }
  else if (index_s_pathway12 == "upsetplot"){
    plots301 <- upsetplot(pathway_enrichment, n = index_s_pathway13)
  }
  
  return(list(plot1 = plots301, data1 = pathway_results))
}
