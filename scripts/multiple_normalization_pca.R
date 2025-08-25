# ================================
# scripts/multiple_normalization_pca.R
# ================================

# ---- MAIN FUNCTION (Normalization + HVGs + Scale + PCA + JackStraw if LogNormalize) ----
datainput_multiple_normalization_pca <- function(
    index_multiple_normalization_pca_input,
    index_multiple_sample_normalization_method,
    index_multiple_sample_scale_factor,
    index_multiple_sample_var_genes,
    index_multiple_sample_normalization_variable_genes,
    index_multiple_sample_pca_dim,       # Seurat PCA dimension (N)
    index_js_dims,                       # Max PCs to test for JackStraw/ScoreJackStraw
    index_js_numrep,                     # JackStraw num.replicate
    index_js_plot_max                    # Max PCs to display in JackStrawPlot
) {
  # 1) Normalize / HVGs / Scale
  if (identical(index_multiple_sample_normalization_method, "LogNormalize")) {
    obj <- NormalizeData(
      index_multiple_normalization_pca_input,
      normalization.method = index_multiple_sample_normalization_method,
      scale.factor = index_multiple_sample_scale_factor
    )
    obj <- FindVariableFeatures(
      obj,
      selection.method = index_multiple_sample_normalization_variable_genes,
      nfeatures = index_multiple_sample_var_genes
    )
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    
  } else if (identical(index_multiple_sample_normalization_method, "SCTransform")) {
    # SCTransform sets/uses variable features internally
    obj <- SCTransform(index_multiple_normalization_pca_input, vst.flavor = "v2", verbose = FALSE)
    DefaultAssay(obj) <- "SCT"
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
    
  } else {
    stop("Unknown normalization method: ", index_multiple_sample_normalization_method)
  }
  
  # 2) PCA (Seurat N)
  npcs_req <- max(2L, as.integer(index_multiple_sample_pca_dim))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = npcs_req, verbose = FALSE)
  
  # 3) Plots (PCA heatmap / elbow / scatter)
  plot1 <- DimHeatmap(obj, dims = 1, cells = 500, balanced = TRUE, fast = FALSE)
  plot2 <- ElbowPlot(obj, ndims = npcs_req)
  plot3 <- DimPlot(obj, reduction = "pca")
  
  group_col <- if ("condition" %in% colnames(obj@meta.data)) {
    "condition"
  } else if ("orig.ident" %in% colnames(obj@meta.data)) {
    "orig.ident"
  } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    "seurat_clusters"
  } else {
    NULL
  }
  plot4 <- if (!is.null(group_col)) {
    DimPlot(obj, reduction = "pca", group.by = group_col)
  } else {
    DimPlot(obj, reduction = "pca")
  }
  
  # 4) JackStraw only if LogNormalize
  js_plot <- NULL
  js_sig_text <- NULL
  js_sig_pcs <- integer(0)
  
  if (identical(index_multiple_sample_normalization_method, "LogNormalize")) {
    max_pcs_available <- ncol(Embeddings(obj, "pca"))
    js_dims   <- min(index_js_dims, 50L, max_pcs_available)  # cap at 50
    js_numrep <- max(20L, as.integer(index_js_numrep))
    plot_max  <- max(1L, min(as.integer(index_js_plot_max), max_pcs_available))
    
    obj <- JackStraw(obj, reduction = "pca",
                     dims = js_dims,
                     num.replicate = js_numrep,
                     verbose = FALSE)
    obj <- ScoreJackStraw(obj, dims = 1:js_dims)
    
    js_plot <- JackStrawPlot(obj, dims = 1:plot_max)
    
    # Significant PCs readout
    js_overall <- JS(obj[["pca"]], slot = "overall")
    if (is.numeric(js_overall)) {
      idx <- which(js_overall < 0.05)
      if (length(idx)) {
        nm <- names(js_overall); if (length(nm) == 0) nm <- paste0("PC", seq_along(js_overall))
        js_sig_pcs <- as.integer(sub("^PC", "", nm[idx]))
      }
    } else if (is.data.frame(js_overall) && all(c("PC","p.val") %in% names(js_overall))) {
      js_sig_pcs <- as.integer(js_overall$PC[js_overall$p.val < 0.05])
    }
    
    if (length(js_sig_pcs)) {
      js_sig_text <- paste0(
        "Significant PCs (p < 0.05): ",
        paste(js_sig_pcs, collapse = ", "),
        "\nTip: start with dims 1:", max(js_sig_pcs),
        " for neighbors/clustering and UMAP/t-SNE."
      )
    } else {
      js_sig_text <- "No PCs reached p < 0.05. Consider increasing num.replicate or dims, or revisit QC."
    }
    
  } else if (identical(index_multiple_sample_normalization_method, "SCTransform")) {
    js_sig_text <- "JackStraw cannot be run on SCTransform-normalized data. Please use ElbowPlot instead."
  }
  
  # 5) Annotate convenience fields
  obj@misc$pc_select_mode <- "seurat"
  obj@misc$pca_npcs_used  <- ncol(Embeddings(obj, "pca"))
  obj@misc$pca_dims       <- seq_len(obj@misc$pca_npcs_used)
  
  # 6) Return payload
  list(
    plot1 = plot1,                 # PCA heatmap
    plot2 = plot2,                 # Elbow
    plot3 = plot3,                 # PCA scatter (samples)
    plot4 = plot4,                 # PCA scatter (groups)
    data1 = obj,                   # Seurat object with PCA
    js_plot = js_plot,             # JackStraw plot (NULL if SCT)
    js_sig_text = js_sig_text,     # significant PCs (or SCT message)
    js_sig_pcs = js_sig_pcs,       # significant PCs (vector)
    
    used_dims = obj@misc$pca_dims
  )
}
