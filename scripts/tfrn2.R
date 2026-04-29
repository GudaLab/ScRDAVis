library(magrittr)
datainput_single_multiple_sample_tfrn2<- function(index_multiple_sample_tfrn2_input, index_multiple_sample_tfrn2_input2,index_multiple_sample_tfrn2_input3,index_s_tfrn11, index_s_tfrn12, index_s_tfrn13, index_s_tfrn14){
  index_s_tfrn14 <- as.numeric(index_s_tfrn14)
  
  fallback_value <- function(x, y) {
    if (is.null(x) || length(x) == 0 || all(is.na(x))) {
      return(y)
    }
    x
  }
  
  empty_tfrn_plot <- function(title_text, body_text) {
    ggplot() +
      annotate("text", x = 0, y = 0.1, label = title_text, fontface = "bold", size = 5) +
      annotate("text", x = 0, y = -0.1, label = body_text, size = 4) +
      xlim(-1, 1) +
      ylim(-1, 1) +
      theme_void()
  }
  
  safe_matrix_column <- function(mat, column_name, fallback_value = 0) {
    if (!is.null(mat) && column_name %in% colnames(mat)) {
      return(mat[, column_name])
    }
    rep(fallback_value, nrow(single_multiple_sample_clustering@meta.data))
  }
  
  safe_tf_network_plot <- function(object, selected_tf, target_type = "both", edge_weight = "Cor", depth = 1, title_text = NULL) {
    tf_regulons <- tryCatch(GetTFRegulons(object), error = function(e) NULL)
    
    has_tf_regulon <- !is.null(tf_regulons) &&
      nrow(tf_regulons) > 0 &&
      "tf" %in% colnames(tf_regulons) &&
      selected_tf %in% tf_regulons$tf
    
    if (!has_tf_regulon) {
      return(empty_tfrn_plot(
        fallback_value(title_text, "TF network plot"),
        paste("No regulon network is available for", selected_tf)
      ))
    }
    
    tryCatch(
      {
        TFNetworkPlot(
          object,
          selected_tfs = selected_tf,
          target_type = target_type,
          edge_weight = edge_weight,
          depth = depth
        ) + ggtitle(fallback_value(title_text, "TF network plot"))
      },
      error = function(e) {
        empty_tfrn_plot(
          fallback_value(title_text, "TF network plot"),
          paste("Network plot not available:", conditionMessage(e))
        )
      }
    )
  }
  
  single_multiple_sample_clustering <- index_multiple_sample_tfrn2_input
  pos_regulon_scores <- index_multiple_sample_tfrn2_input2
  neg_regulon_scores <- index_multiple_sample_tfrn2_input3
  cur_tf <- index_s_tfrn11
  
  # add the regulon scores to the Seurat metadata
  single_multiple_sample_clustering$pos_regulon_score <- safe_matrix_column(pos_regulon_scores, cur_tf)
  single_multiple_sample_clustering$neg_regulon_score <- safe_matrix_column(neg_regulon_scores, cur_tf)
  
  
  # plot using FeaturePlot
  plots905 <- tryCatch(
    {
      FeaturePlot(single_multiple_sample_clustering, feature = cur_tf) + umap_theme()
    },
    error = function(e) {
      empty_tfrn_plot(
        paste("Feature plot:", cur_tf),
        paste("Feature plot not available:", conditionMessage(e))
      )
    }
  )
  
  
  #Network Visualization
  # select TF of interest
  #Visualize regulons
  plots906 <- tryCatch(
    {
      RegulonBarPlot(single_multiple_sample_clustering, selected_tf = cur_tf, top_n = index_s_tfrn12)
    },
    error = function(e) {
      empty_tfrn_plot(
        paste("Regulon bar plot:", cur_tf),
        paste("Regulon bar plot not available:", conditionMessage(e))
      )
    }
  )
  
  hub_df1 <- tryCatch(
    GetHubGenes(single_multiple_sample_clustering, n_hubs = 20),
    error = function(e) NULL
  )
  cur_mod <- character()
  cur_mod_genes <- character()
  if (!is.null(hub_df1) && nrow(hub_df1) > 0) {
    cur_mod <- subset(hub_df1, gene_name == cur_tf) %>% .$module %>% as.character
    cur_mod_genes <- subset(hub_df1, module == cur_mod) %>% .$gene_name
  }
  
  # plot with default settings
  plots907 <- safe_tf_network_plot(single_multiple_sample_clustering, selected_tf = cur_tf, target_type = "positive", edge_weight = index_s_tfrn13, depth = index_s_tfrn14, title_text = "Positive targets")
  plots908 <- safe_tf_network_plot(single_multiple_sample_clustering, selected_tf = cur_tf, target_type = "negative", edge_weight = index_s_tfrn13, depth = index_s_tfrn14, title_text = "Negative targets")
  plots909 <- safe_tf_network_plot(single_multiple_sample_clustering, selected_tf = cur_tf, target_type = "both", edge_weight = index_s_tfrn13, depth = index_s_tfrn14, title_text = "Pos & Neg targets")
  #plots910 <- TFNetworkPlot(single_multiple_sample_clustering, selected_tfs=cur_tf, label_TFs=0, label_genes=cur_mod_genes, depth=2) + ggtitle('Hub genes in the same module')
  # plots907 | plots908 | plots909
  
  return(list(plot1 = plots905, plot2 = plots906, plot3 = plots907, plot4 = plots908, plot5 = plots909, data1=single_multiple_sample_clustering))
}
