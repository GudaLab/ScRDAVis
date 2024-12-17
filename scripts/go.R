datainput_single_multiple_sample_go<- function(index_multiple_sample_go_input, index_subclustering_multiple_sample_go_input, index_multiple_sample_go_input2, index_subclustering_multiple_sample_go_input2, index_multiple_sample_go_input3, index_subclustering_multiple_sample_go_input3, index_s_go1, index_s_go2, index_s_go3, index_s_go4, index_s_go5, index_s_go6, index_s_go7, index_s_go8, index_s_go9, index_s_go10, index_s_go11, index_s_go12, index_s_go13, index_s_go14){
 
  if (index_s_go1 == "gene_name_list"){
    deg_genes <- unlist(strsplit(index_s_go14, ","))
  }
  
  else if (index_s_go1 == "multiple_sample" & index_s_go2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_multiple_sample_go_input
    single_sample_clustering_markers <- index_multiple_sample_go_input3 
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_go3 & p_val_adj < index_s_go4))
    deg_genes <- deg_genes$gene
  }
  else if (index_s_go1 == "multiple_sample_subclustering" & index_s_go2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_subclustering_multiple_sample_go_input 
    single_sample_clustering_markers <- index_subclustering_multiple_sample_go_input3 
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_go3 & p_val_adj < index_s_go4))
    deg_genes <- deg_genes$gene
  }
  
  
  else if (index_s_go1 == "multiple_sample" & index_s_go2 == "predicted"){
    if (index_multiple_sample_go_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_go_input2
      single_sample_clustering_markers <- index_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_go3  & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
    else if (index_multiple_sample_go_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_go_input2
      single_sample_clustering_markers <- index_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
    else if (index_multiple_sample_go_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_go_input2
      single_sample_clustering_markers <- index_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
    else if (index_multiple_sample_go_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_go_input2
      single_sample_clustering_markers <- index_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
  }
  
  else if (index_s_go1 == "multiple_sample_subclustering" & index_s_go2 == "predicted"){
    if (index_subclustering_multiple_sample_go_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_go_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_go_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_go_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_go_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_go_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
    else if (index_subclustering_multiple_sample_go_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_go_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_go_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_go_input3 
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_go3 & p_val_adj < index_s_go4))
      deg_genes <- deg_genes$gene
    }
  }
  
  entrez_ids <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = index_s_go5)

  entrez_gene_list <- entrez_ids$ENTREZID

  go_enrichment <- enrichGO(gene = entrez_gene_list,
                            OrgDb = index_s_go5,
                            ont = index_s_go6,
                            pAdjustMethod = index_s_go7,
                            pvalueCutoff = index_s_go8,
                            qvalueCutoff = index_s_go9,
                            minGSSize = index_s_go10,
                            maxGSSize = index_s_go11,
                            readable = TRUE)

  if(index_s_go12 == "dotplot"){
  plots201 <- dotplot(go_enrichment, showCategory = index_s_go13)
  }
  else if (index_s_go12 == "barplot"){
  plots201 <- barplot(go_enrichment, showCategory = index_s_go13)
  }
  else if (index_s_go12 == "goplot"){
    plots201 <- goplot(go_enrichment, showCategory = index_s_go13)
  }
  else if (index_s_go12 == "cnetplot"){
    plots201 <- cnetplot(go_enrichment, showCategory = index_s_go13)
  }
  else if (index_s_go12 == "upsetplot"){
    plots201 <- upsetplot(go_enrichment, n = index_s_go13)
  }

  go_enrichment_result <- go_enrichment@result

  return(list(plot1 = plots201, data1 = go_enrichment_result))
}
