datainput_single_multiple_sample_cccn<- function(index_multiple_sample_cccn_input, index_subclustering_multiple_sample_cccn_input, index_multiple_sample_cccn_input2, index_subclustering_multiple_sample_cccn_input2, index_multiple_sample_normalization_method_cccn, index_subclustering_multiple_sample_normalization_method_cccn, index_s_cccn1, index_s_cccn2, index_s_cccn3){

  if (index_s_cccn1 == "multiple_sample" & index_s_cccn2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_multiple_sample_cccn_input 
      }
  else if (index_s_cccn1 == "multiple_sample_subclustering" & index_s_cccn2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_subclustering_multiple_cccn_sample_input  
  }
  else if (index_s_cccn1 == "multiple_sample" & index_s_cccn2 == "predicted"){
    if (index_multiple_sample_cccn_input2 == "sctype_classification"){
    single_multiple_sample_clustering <- index_multiple_sample_cccn_input
    Idents(single_multiple_sample_clustering) <- index_multiple_sample_cccn_input2
    }
    else if (index_multiple_sample_cccn_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_cccn_input2
    }
    else if (index_multiple_sample_cccn_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_cccn_input2
    }
    else if (index_multiple_sample_cccn_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_cccn_input2
    }
     }
  else if (index_s_cccn1 == "multiple_sample_subclustering" & index_s_cccn2 == "predicted"){
    single_multiple_sample_clustering <- index_subclustering_multiple_sample_cccn_input
    if (index_subclustering_multiple_sample_cccn_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_cccn_input2
    }
    else if (index_subclustering_multiple_sample_cccn_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_cccn_input2
    }
    else if (index_subclustering_multiple_sample_cccn_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_cccn_input2
    }
    else if (index_subclustering_multiple_sample_cccn_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_cccn_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_cccn_input2
    }
  }
  
  if (index_s_cccn1 == "multiple_sample" & index_multiple_sample_normalization_method_cccn == "LogNormalize"){
  #calculate gene-cluster specificity scores
  sg = genesorteR::sortGenes(single_multiple_sample_clustering@assays$RNA$scale.data, Idents(single_multiple_sample_clustering))
  }
  else if(index_s_cccn1 == "multiple_sample" & index_multiple_sample_normalization_method_cccn == "SCTransform"){
  sg = genesorteR::sortGenes(single_multiple_sample_clustering@assays$SCT$scale.data, Idents(single_multiple_sample_clustering))  
  }
  else if (index_s_cccn1 == "multiple_sample_subclustering" &index_subclustering_multiple_sample_normalization_method_cccn == "LogNormalize"){
  #calculate gene-cluster specificity scores
  sg = genesorteR::sortGenes(single_multiple_sample_clustering@assays$RNA$scale.data, Idents(single_multiple_sample_clustering))
  }
  else if(index_s_cccn1 == "multiple_sample_subclustering" &index_subclustering_multiple_sample_normalization_method_cccn == "SCTransform"){
  sg = genesorteR::sortGenes(single_multiple_sample_clustering@assays$SCT$scale.data, Idents(single_multiple_sample_clustering))  
  }
  #get highly variable genes between cell clusters
  pv = genesorteR::getPValues(sg)
  hvg = names(which(apply(pv$adjpval, 1, function(x) any(x < 0.05))))
  #calculate cluster correlations based on highly variable genes corMethod="pearson","spearman","kendall"
  pc <-  genesorteR::plotCorrelationHeat(sg, markers = hvg, corMethod=index_s_cccn3, outs = TRUE)$corr
  pc1 <-  as.data.frame(pc)
  plots901 <- pheatmap(pc, main = "Correlation Heatmap", display_numbers = TRUE, number_color = "black", fontsize_number = 8)
  plots903 <- as.ggplot(plots901)
  #set diagonal entries to 0 (network lingo: no self-loops)
  diag(pc) = 0
  #set negative entries to 0 (if two clusters' correlation is below 0.1, they are not connected. You can experiment with different cutoffs.)
  pc[pc < 0.1] = 0
  #Step Three: Visualize the Correlation Network
  #create an igraph network object from the weighted adjacency matrix stored in pc
  net = igraph::graph_from_adjacency_matrix(pc, weighted = TRUE)
  #remove multiple edges (meaning keep only one connection between each two cell clusters)
  net = igraph::simplify(net, edge.attr.comb = "max")
  #let's give the clusters names (identified in the Seurat tutorial)
  #cluster_name = single_multiple_sample_clustering$
  #Some information to use in our plots
  Correlation = E(net)$weight
  Percent_of_Cells = as.vector(sg$classProb*100)
  #plot network with ggraph
  set.seed(111)
  lay = ggraph::create_layout(net, layout = "fr")
  plots902 <-  ggraph(lay) + 
    geom_edge_link(aes(alpha = Correlation), edge_colour = "gray") + 
    geom_node_point(aes(size = Percent_of_Cells), colour = "blue") + 
    geom_node_text(aes(label = name), repel=TRUE) +
    theme(panel.background = element_blank())
  
  return(list(plot1 = plots903, plot2 = plots902, data1=pc1)) 
  
}
  
  