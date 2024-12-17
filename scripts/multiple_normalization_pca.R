datainput_multiple_normalization_pca <- function(index_multiple_normalization_pca_input, index_multiple_sample_normalization_method, index_multiple_sample_scale_factor, index_multiple_sample_var_genes, index_multiple_sample_normalization_variable_genes, index_multiple_sample_pca_dim){
 
  if (index_multiple_sample_normalization_method == "LogNormalize")
  { 
  multiple_sample_normalized <- NormalizeData(index_multiple_normalization_pca_input, normalization.method = index_multiple_sample_normalization_method, scale.factor = index_multiple_sample_scale_factor)
  multiple_sample_normalized <- FindVariableFeatures(multiple_sample_normalized, selection.method = index_multiple_sample_normalization_variable_genes, nfeatures = index_multiple_sample_var_genes)
  multiple_sample_normalized <- ScaleData(multiple_sample_normalized, verbose = FALSE)
  multiple_sample_normalized <- RunPCA(multiple_sample_normalized, features = VariableFeatures(object = multiple_sample_normalized), npcs = index_multiple_sample_pca_dim)
  plots12 <- DimHeatmap(multiple_sample_normalized, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
  plots13 <- ElbowPlot(multiple_sample_normalized)
  plots14 <- DimPlot(multiple_sample_normalized, reduction = "pca")
  plots15 <- DimPlot(multiple_sample_normalized, reduction = "pca", group.by = "condition")
  }
  

  else if (index_multiple_sample_normalization_method == "SCTransform")
  {
  multiple_sample_normalized <- SCTransform(index_multiple_normalization_pca_input, vst.flavor="v2")
  multiple_sample_normalized <- ScaleData(multiple_sample_normalized, verbose = FALSE)
  multiple_sample_normalized <- RunPCA(multiple_sample_normalized, features = VariableFeatures(object = multiple_sample_normalized), npcs = index_multiple_sample_pca_dim)
  plots12 <- DimHeatmap(multiple_sample_normalized, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
  plots13 <- ElbowPlot(multiple_sample_normalized)
  plots14 <- DimPlot(multiple_sample_normalized, reduction = "pca")
  plots15 <- DimPlot(multiple_sample_normalized, reduction = "pca", group.by = "condition")
  }
  
  return(list(plot1 = plots12, plot2 = plots13, plot3 = plots14, plot4 = plots15, data1 = multiple_sample_normalized))
}