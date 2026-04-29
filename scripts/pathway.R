datainput_single_multiple_sample_pathway<- function(index_multiple_sample_pathway_input, index_subclustering_multiple_sample_pathway_input, index_multiple_sample_pathway_input2, index_subclustering_multiple_sample_pathway_input2, index_multiple_sample_pathway_input3, index_subclustering_multiple_sample_pathway_input3, index_s_pathway1, index_s_pathway2, index_s_pathway3, index_s_pathway4, index_s_pathway5, index_s_pathway6, index_s_pathway7, index_s_pathway8, index_s_pathway9, index_s_pathway10, index_s_pathway11, index_s_pathway12, index_s_pathway13, index_s_pathway14){
  make_placeholder_plot <- function(message_text) {
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = message_text, size = 5) +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, 1)
  }

  empty_pathway_result <- function(message_text) {
    list(
      plot1 = make_placeholder_plot(message_text),
      data1 = data.frame(Message = message_text, stringsAsFactors = FALSE)
    )
  }

  reformat_kegg_gene_ids <- function(pathway_result, entrez_map) {
    if (nrow(pathway_result) == 0 || !"geneID" %in% colnames(pathway_result)) {
      return(pathway_result)
    }

    gene_id_values <- as.character(pathway_result$geneID)
    pathway_result$geneID <- vapply(
      gene_id_values,
      FUN.VALUE = character(1),
      function(gene_string) {
        if (is.na(gene_string) || !nzchar(gene_string)) {
          return("")
        }
        genes <- unlist(strsplit(gene_string, "/", fixed = TRUE), use.names = FALSE)
        genes <- genes[nzchar(genes)]
        gene_symbols <- entrez_map$SYMBOL[match(genes, entrez_map$ENTREZID)]
        gene_symbols <- gene_symbols[!is.na(gene_symbols) & nzchar(gene_symbols)]
        if (length(gene_symbols) == 0) {
          return(gene_string)
        }
        paste(unique(gene_symbols), collapse = "/")
      }
    )
    pathway_result
  }

  generate_pathway_plot <- function(pathway_enrichment, plot_type, show_category, fallback_message) {
    plot_object <- tryCatch(
      {
        if (plot_type == "dotplot") {
          dotplot(pathway_enrichment, showCategory = show_category)
        } else if (plot_type == "barplot") {
          barplot(pathway_enrichment, showCategory = show_category)
        } else if (plot_type == "goplot") {
          goplot(pathway_enrichment, showCategory = show_category)
        } else if (plot_type == "cnetplot") {
          cnetplot(pathway_enrichment, showCategory = show_category)
        } else if (plot_type == "upsetplot") {
          upsetplot(pathway_enrichment, n = show_category)
        } else {
          make_placeholder_plot("Selected pathway plot type is not supported.")
        }
      },
      error = function(e) {
        make_placeholder_plot(fallback_message)
      }
    )

    plot_object
  }

  deg_genes <- character(0)

  if (identical(index_s_pathway1, "gene_name_list")){
    gene_text <- if (length(index_s_pathway14) == 0 || is.null(index_s_pathway14)) "" else as.character(index_s_pathway14)[1]
    deg_genes <- unlist(strsplit(gene_text, ",", fixed = TRUE), use.names = FALSE)
  }
  else if (identical(index_s_pathway1, "multiple_sample") && identical(index_s_pathway2, "seurat_clusters")){
    single_sample_clustering_markers <- index_multiple_sample_pathway_input3
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
    deg_genes <- deg_genes$gene
  }
  else if (identical(index_s_pathway1, "multiple_sample_subclustering") && identical(index_s_pathway2, "seurat_clusters")){
    single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3
    deg_genes <- subset(single_sample_clustering_markers, (cluster ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
    deg_genes <- deg_genes$gene
  }
  else if (identical(index_s_pathway1, "multiple_sample") && identical(index_s_pathway2, "predicted")){
    if (identical(index_multiple_sample_pathway_input2, "sctype_classification")){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_pathway3  & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (identical(index_multiple_sample_pathway_input2, "singleR_labels")){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (identical(index_multiple_sample_pathway_input2, "GPTCelltype")){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (identical(index_multiple_sample_pathway_input2, "cell_type")){
      single_multiple_sample_clustering <- index_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
  }
  else if (identical(index_s_pathway1, "multiple_sample_subclustering") && identical(index_s_pathway2, "predicted")){
    if (identical(index_subclustering_multiple_sample_pathway_input2, "sctype_classification")){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, sctype_classification) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (sctype_classification ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (identical(index_subclustering_multiple_sample_pathway_input2, "singleR_labels")){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, singleR_labels) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (singleR_labels ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (identical(index_subclustering_multiple_sample_pathway_input2, "GPTCelltype")){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, GPTCelltype) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (GPTCelltype ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
    else if (identical(index_subclustering_multiple_sample_pathway_input2, "cell_type")){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_pathway_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_pathway_input2
      single_sample_clustering_markers <- index_subclustering_multiple_sample_pathway_input3
      export_df <- single_multiple_sample_clustering@meta.data %>% dplyr::select(seurat_clusters, cell_type) %>% distinct()
      single_sample_clustering_markers <- merge(single_sample_clustering_markers, export_df, by.x = "cluster", by.y = "seurat_clusters")
      deg_genes <- subset(single_sample_clustering_markers, (cell_type ==  index_s_pathway3 & p_val_adj < index_s_pathway4))
      deg_genes <- deg_genes$gene
    }
  }

  if (is.null(deg_genes)) {
    deg_genes <- character(0)
  }
  deg_genes <- as.character(deg_genes)
  deg_genes <- trimws(deg_genes)
  deg_genes <- unique(deg_genes[nzchar(deg_genes) & !is.na(deg_genes)])

  if (length(deg_genes) == 0) {
    return(empty_pathway_result("No genes were available for pathway enrichment."))
  }

  entrez_ids <- tryCatch(
    bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = index_s_pathway5),
    error = function(e) NULL
  )

  if (is.null(entrez_ids) || nrow(entrez_ids) == 0) {
    return(empty_pathway_result("No input genes could be mapped to Entrez IDs."))
  }

  entrez_ids <- entrez_ids[!is.na(entrez_ids$ENTREZID) & nzchar(entrez_ids$ENTREZID), , drop = FALSE]
  entrez_ids <- entrez_ids[!duplicated(entrez_ids$ENTREZID), , drop = FALSE]

  if (nrow(entrez_ids) == 0) {
    return(empty_pathway_result("No valid Entrez IDs were available for pathway enrichment."))
  }

  entrez_gene_list <- entrez_ids$ENTREZID
  pathway_enrichment <- NULL
  pathway_results <- data.frame()

  if (identical(index_s_pathway6, "KEGG")){
    kegg_organism <- switch(
      index_s_pathway5,
      "org.Hs.eg.db" = "hsa",
      "org.Mm.eg.db" = "mmu",
      "org.Rn.eg.db" = "rno",
      NULL
    )

    if (is.null(kegg_organism)) {
      return(empty_pathway_result("Selected organism is not supported for KEGG enrichment."))
    }

    pathway_enrichment <- tryCatch(
      enrichKEGG(
        gene = entrez_gene_list,
        organism = kegg_organism,
        pAdjustMethod = index_s_pathway7,
        pvalueCutoff = index_s_pathway8,
        qvalueCutoff = index_s_pathway9,
        minGSSize = index_s_pathway10,
        maxGSSize = index_s_pathway11
      ),
      error = function(e) NULL
    )

    if (!is.null(pathway_enrichment)) {
      pathway_enrichment@result <- reformat_kegg_gene_ids(pathway_enrichment@result, entrez_ids)
      pathway_results <- as.data.frame(pathway_enrichment@result)
      if (nrow(pathway_results) > 0 && all(c("Description", "geneID") %in% colnames(pathway_results))) {
        pathway_results <- pathway_results %>% dplyr::select(Description, geneID, everything())
      }
    }
  }
  else if (identical(index_s_pathway6, "Reactome")){
    reactome_organism <- switch(
      index_s_pathway5,
      "org.Hs.eg.db" = "human",
      "org.Mm.eg.db" = "mouse",
      "org.Rn.eg.db" = "rat",
      NULL
    )

    if (is.null(reactome_organism)) {
      return(empty_pathway_result("Selected organism is not supported for Reactome enrichment."))
    }

    pathway_enrichment <- tryCatch(
      enrichPathway(
        gene = entrez_gene_list,
        organism = reactome_organism,
        pAdjustMethod = index_s_pathway7,
        pvalueCutoff = index_s_pathway8,
        qvalueCutoff = index_s_pathway9,
        minGSSize = index_s_pathway10,
        maxGSSize = index_s_pathway11,
        readable = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(pathway_enrichment)) {
      pathway_results <- as.data.frame(pathway_enrichment@result)
    }
  }

  if (is.null(pathway_enrichment) || nrow(pathway_results) == 0) {
    return(empty_pathway_result("No enriched pathways were found for the selected settings."))
  }

  plots301 <- generate_pathway_plot(
    pathway_enrichment = pathway_enrichment,
    plot_type = index_s_pathway12,
    show_category = index_s_pathway13,
    fallback_message = "The pathway summary table is available, but this plot could not be generated for the selected result."
  )

  return(list(plot1 = plots301, data1 = pathway_results))
}
