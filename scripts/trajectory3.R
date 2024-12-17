datainput_single_multiple_sample_trajectory3 <- function(index_trajectory3_input1, index_trajectory3_input2, index_s_trajectory15, index_s_trajectory16){
  single_multiple_sample_clustering <- index_trajectory3_input1
  cds <- index_trajectory3_input2
  
table_deg <- graph_test(cds, neighbor_graph = index_s_trajectory15)
table_deg <- table_deg %>% arrange(q_value) %>% dplyr::filter(status == "OK")
#Add pseudotime values into the seuratobject
single_multiple_sample_clustering$pseudotime <- pseudotime(cds)
plots104 <- FeaturePlot(single_multiple_sample_clustering, features = "pseudotime", label=index_s_trajectory16)
return(list(data1 = single_multiple_sample_clustering, data2 = cds, plot1 = plots104, data1 = table_deg))
}