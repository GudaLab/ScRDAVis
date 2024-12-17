datainput_single_multiple_sample_trajectory4 <- function(index_trajectory4_input1, index_trajectory4_input2, index_trajectory4_input3, index_s_trajectory17, index_s_trajectory18){
  single_multiple_sample_clustering <- index_trajectory4_input1
  cds <- index_trajectory4_input2
  table_deg<-index_trajectory4_input3
  if(index_s_trajectory17 == "gene_name_list"){
    top_genes <- unlist(strsplit(index_s_trajectory18, ","))
    top_5 <- row.names(subset(fData(cds),gene_short_name %in% top_genes))
  }
  else{
    index_s_trajectory17 <- as.numeric(index_s_trajectory17)
  top_5 <- head(table_deg$gene_short_name, n=index_s_trajectory17)
  }
plots105 <- FeaturePlot(single_multiple_sample_clustering, features = top_5)

return(list(data1 = single_multiple_sample_clustering, data2 = cds, plot1 = plots105))

}