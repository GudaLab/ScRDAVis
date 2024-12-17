library(magrittr)
datainput_single_multiple_sample_tfrn2<- function(index_multiple_sample_tfrn2_input, index_multiple_sample_tfrn2_input2,index_multiple_sample_tfrn2_input3,index_s_tfrn11, index_s_tfrn12, index_s_tfrn13, index_s_tfrn14){
  index_s_tfrn14 <- as.numeric(index_s_tfrn14)
  
  single_multiple_sample_clustering <- index_multiple_sample_tfrn2_input
  pos_regulon_scores <- index_multiple_sample_tfrn2_input2
  neg_regulon_scores <- index_multiple_sample_tfrn2_input3
  cur_tf <- index_s_tfrn11
  
  # add the regulon scores to the Seurat metadata
  single_multiple_sample_clustering$pos_regulon_score <- pos_regulon_scores[,cur_tf]
  single_multiple_sample_clustering$neg_regulon_score <- neg_regulon_scores[,cur_tf]
  
  
  # plot using FeaturePlot
  plots905 <- FeaturePlot(single_multiple_sample_clustering, feature=cur_tf) + umap_theme()
  
  
  #Network Visualization
  # select TF of interest
  #Visualize regulons
  plots906 <- RegulonBarPlot(single_multiple_sample_clustering, selected_tf=cur_tf, top_n = index_s_tfrn12)
  
  # plot with default settings
  p <- TFNetworkPlot(single_multiple_sample_clustering, selected_tfs=cur_tf)
  
  
  hub_df1 <- GetHubGenes(single_multiple_sample_clustering, n_hubs = 20)
  cur_mod <- subset(hub_df1, gene_name == cur_tf) %>% .$module %>% as.character
  cur_mod_genes <- subset(hub_df1, module == cur_mod) %>% .$gene_name
  
  # plot with default settings
  plots907 <- TFNetworkPlot(single_multiple_sample_clustering, selected_tfs=cur_tf, target_type='positive', edge_weight = index_s_tfrn13, depth=index_s_tfrn14) + ggtitle("Positive targets")
  plots908 <- TFNetworkPlot(single_multiple_sample_clustering, selected_tfs=cur_tf, target_type ='negative', edge_weight = index_s_tfrn13, depth=index_s_tfrn14) + ggtitle("Negative targets")
  plots909 <- TFNetworkPlot(single_multiple_sample_clustering, selected_tfs=cur_tf, target_type ='both', edge_weight = index_s_tfrn13, depth=index_s_tfrn14) + ggtitle("Pos & Neg targets")
  #plots910 <- TFNetworkPlot(single_multiple_sample_clustering, selected_tfs=cur_tf, label_TFs=0, label_genes=cur_mod_genes, depth=2) + ggtitle('Hub genes in the same module')
  # plots907 | plots908 | plots909
  
  return(list(plot1 = plots905, plot2 = plots906, plot3 = plots907, plot4 = plots908, plot5 = plots909, data1=single_multiple_sample_clustering))
}