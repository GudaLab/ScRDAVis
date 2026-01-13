library(motifmatchr)
library(TFBSTools)

datainput_single_multiple_sample_tfrn1<- function(index_multiple_sample_tfrn1_input, index_s_tfrn1, index_s_tfrn2, index_s_tfrn3, index_s_tfrn4, index_s_tfrn5, index_s_tfrn6, index_s_tfrn7, index_s_tfrn8){

  single_multiple_sample_clustering <- index_multiple_sample_tfrn1_input
  
  #pfm_core <- TFBSTools::getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))
  
  
  JASPAR2024 <- JASPAR2024()
  sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(JASPAR2024))

  pfm_core <- TFBSTools::getMatrixSet(x = sq24,  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

  #Construct TF Regulatory Network
  # run the motif scan
  if (index_s_tfrn1 == "EnsDb.Hsapiens.v86"){
    library("EnsDb.Hsapiens.v86")
    library("BSgenome.Hsapiens.UCSC.hg38")
  single_multiple_sample_clustering <- MotifScan(single_multiple_sample_clustering, species_genome = 'hg38', pfm = pfm_core, EnsDb = EnsDb.Hsapiens.v86)
  }
  
  else if(index_s_tfrn1 == "EnsDb.Mmusculus.v79"){
    library("EnsDb.Mmusculus.v79")
    library("BSgenome.Mmusculus.UCSC.mm10")
  single_multiple_sample_clustering <- MotifScan(single_multiple_sample_clustering, species_genome = 'mm10', pfm = pfm_core, EnsDb = EnsDb.Mmusculus.v79)
  }
  # get the motif df:
  motif_df <- GetMotifs(single_multiple_sample_clustering)
  
  
  # keep all TFs, and then remove all genes from the grey module
  tf_genes <- unique(motif_df$gene_name)
  modules <- GetModules(single_multiple_sample_clustering)
  nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
  genes_use <- c(tf_genes, nongrey_genes)
  
  # update the gene list and re-run SetDatExpr
  single_multiple_sample_clustering <- SetWGCNAGenes(single_multiple_sample_clustering, genes_use)
  #single_multiple_sample_clustering <- SetDatExpr(single_multiple_sample_clustering, group.by = cluster_types, group_name=cluster_number)
  
  
  # define model params:
  model_params <- list(objective = 'reg:squarederror',  max_depth = index_s_tfrn2,  eta = index_s_tfrn3,  nthread=64,  alpha=index_s_tfrn4)
  
  callbacks <- if (utils::packageVersion("xgboost") >= "2.0.0") {
    list(xgboost::xgb.cb.cv.predict(save_models = TRUE))
  } else {
    list(xgboost::cb.cv.predict(save_models = TRUE))
  }
  
  # construct the TF network
  single_multiple_sample_clustering <- ConstructTFNetwork(single_multiple_sample_clustering, model_params,
                                                          callbacks = callbacks)
  results <- GetTFNetwork(single_multiple_sample_clustering)
  head(results)
  
  
  
  #Define TF Regulons
  single_multiple_sample_clustering <- AssignTFRegulons(single_multiple_sample_clustering, strategy = "A",  reg_thresh = index_s_tfrn5,  n_tfs = index_s_tfrn6)
  #Calculate regulon expression signatures
  # positive regulons
  single_multiple_sample_clustering <- RegulonScores(single_multiple_sample_clustering, target_type = 'positive', cor_thresh = index_s_tfrn7, ncores=8)
  # negative regulons
  single_multiple_sample_clustering <- RegulonScores(single_multiple_sample_clustering, target_type = 'negative', cor_thresh = index_s_tfrn8, ncores=8)
  
  
  # access the results:
  pos_regulon_scores <- GetRegulonScores(single_multiple_sample_clustering, target_type='positive')
  neg_regulon_scores <- GetRegulonScores(single_multiple_sample_clustering, target_type='negative')
  #both_regulon_scores <- GetRegulonScores(single_multiple_sample_clustering, target_type='both')
  
  #Module regulatory networks
  #ModuleRegulatoryNetworkPlot
  plots901 <- ModuleRegulatoryNetworkPlot(single_multiple_sample_clustering, feature='positive', high_color='orange2')
  #circle, stress, eigen, linear, graphopt
  plots902 <- ModuleRegulatoryNetworkPlot(single_multiple_sample_clustering, feature='negative', high_color='dodgerblue')
  plots903 <- ModuleRegulatoryNetworkPlot(single_multiple_sample_clustering, feature='delta')
  # same plot with additional options
  plots904 <- ModuleRegulatoryNetworkPlot(single_multiple_sample_clustering, feature='delta', umap_background=TRUE, label_modules=FALSE)
  #plots901 | plots902 | plots903 | plots904
  
  
  tf_regulons <- GetTFRegulons(single_multiple_sample_clustering)
  
  # Step 1: Get column names from reference_table
  reference_genes <- colnames(pos_regulon_scores)
  
  # Step 2: Filter gene_table based on the column names of reference_table
  filtered_genes <- tf_regulons[tf_regulons$tf %in% reference_genes, ]
  filtered_genes <- filtered_genes[order(filtered_genes$Cor, decreasing = TRUE), ]
  
  filtered_tf<-filtered_genes$tf
  filtered_tf1<-as.data.frame(unique(filtered_tf))
  colnames(filtered_tf1) <-"TF"
  
  detach("package:TFBSTools", unload = TRUE)
  detach("package:motifmatchr", unload = TRUE)
return(list(plot1 = plots901, plot2 = plots902, plot3 = plots903, plot4 = plots904, data1=tf_regulons, data2=filtered_tf1, data3=single_multiple_sample_clustering, data4=pos_regulon_scores, data5=neg_regulon_scores))  
}
