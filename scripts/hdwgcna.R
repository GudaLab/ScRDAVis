datainput_single_multiple_sample_hdwgcna<- function(index_multiple_sample_hdwgcna_input, index_subclustering_multiple_sample_hdwgcna_input, index_multiple_sample_hdwgcna_input2, index_subclustering_multiple_sample_hdwgcna_input2, index_multiple_sample_normalization_method_hdwgcna, index_subclustering_multiple_sample_normalization_method_hdwgcna, index_s_hdwgcna1, index_s_hdwgcna2, index_s_hdwgcna3, index_s_hdwgcna4, index_s_hdwgcna5, index_s_hdwgcna6, index_s_hdwgcna7, index_s_hdwgcna8, index_s_hdwgcna9, index_s_hdwgcna10, index_s_hdwgcna11, index_s_hdwgcna12, index_s_hdwgcna13, index_s_hdwgcna14){
  index_s_hdwgcna11 <-as.logical(index_s_hdwgcna11)
  index_s_hdwgcna14 <-as.logical(index_s_hdwgcna14)
  
  if (index_s_hdwgcna1 == "multiple_sample" & index_s_hdwgcna2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_multiple_sample_hdwgcna_input
    cluster_types <- "seurat_clusters"
    cluster_number <- index_s_hdwgcna3
    }
  else if (index_s_hdwgcna1 == "multiple_sample_subclustering" & index_s_hdwgcna2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_subclustering_multiple_sample_hdwgcna_input 
    cluster_types <- "seurat_clusters"
    cluster_number <- index_s_hdwgcna3
  }
 else if (index_s_hdwgcna1 == "multiple_sample" & index_s_hdwgcna2 == "predicted"){
    if (index_multiple_sample_hdwgcna_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_hdwgcna_input2
      cluster_types <- "sctype_classification"
      cluster_number <- index_s_hdwgcna3
    }
    else if (index_multiple_sample_hdwgcna_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_hdwgcna_input2
      cluster_types <- "singleR_labels"
      cluster_number <- index_s_hdwgcna3
    }
    else if (index_multiple_sample_hdwgcna_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_hdwgcna_input2
      cluster_types <- "GPTCelltype"
      cluster_number <- index_s_hdwgcna3
    }
    else if (index_multiple_sample_hdwgcna_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_multiple_sample_hdwgcna_input2
      cluster_types <- "cell_type"
      cluster_number <- index_s_hdwgcna3
     }
  }
  else if (index_s_hdwgcna1 == "multiple_sample_subclustering" & index_s_hdwgcna2 == "predicted"){
    if (index_subclustering_multiple_sample_hdwgcna_input2 == "sctype_classification"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_hdwgcna_input2
      cluster_types <- "sctype_classification"
      cluster_number <- index_s_hdwgcna3
      }
    else if (index_subclustering_multiple_sample_hdwgcna_input2 == "singleR_labels"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_hdwgcna_input2
      cluster_types <- "singleR_labels"
      cluster_number <- index_s_hdwgcna3
      }
    else if (index_subclustering_multiple_sample_hdwgcna_input2 == "GPTCelltype"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_hdwgcna_input2
      cluster_types <- "GPTCelltype"
      cluster_number <- index_s_hdwgcna3
      }
    else if (index_subclustering_multiple_sample_hdwgcna_input2 == "cell_type"){
      single_multiple_sample_clustering <- index_subclustering_multiple_sample_hdwgcna_input
      Idents(single_multiple_sample_clustering) <- index_subclustering_multiple_sample_hdwgcna_input2
      cluster_types <- "cell_type"
      cluster_number <- index_s_hdwgcna3
      }
  }
  
  enableWGCNAThreads(nThreads = 8)
  
  plots801 <- DimPlot(single_multiple_sample_clustering, group.by=cluster_types, label=TRUE) + NoLegend()
  
  single_multiple_sample_clustering <- SetupForWGCNA(
    single_multiple_sample_clustering,
    gene_select = "fraction", # the gene selection approach (variable, fraction)
    fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
    wgcna_name = "scrdavis" # the name of the hdWGCNA experiment
  )
    
  # construct metacells  in each group
  single_multiple_sample_clustering  <- MetacellsByGroups(
    seurat_obj = single_multiple_sample_clustering,
    group.by = cluster_types, # specify the columns in seurat_obj@meta.data to group by
    reduction = index_s_hdwgcna4, # select the dimensionality reduction to perform KNN on
    min_cells = index_s_hdwgcna6,
    k = index_s_hdwgcna5, # nearest-neighbors parameter
    max_shared = index_s_hdwgcna7, # maximum number of shared cells between two metacells
    ident.group = cluster_types, # set the Idents of the metacell seurat object
    target_metacells = index_s_hdwgcna8
  )
  # normalize metacell expression matrix:
  single_multiple_sample_clustering <- NormalizeMetacells(single_multiple_sample_clustering)
  
  if (index_s_hdwgcna1 == "multiple_sample" & index_multiple_sample_normalization_method_hdwgcna == "LogNormalize"){
    #Co-expression network analysis
    single_multiple_sample_clustering <- SetDatExpr(
      single_multiple_sample_clustering,
      group_name = cluster_number, # the name of the group of interest in the group.by column
      group.by=cluster_types, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
      assay = 'RNA', # using RNA assay
      slot = 'data' # using normalized data
    )
  }
  else if(index_s_hdwgcna1 == "multiple_sample" & index_multiple_sample_normalization_method_hdwgcna == "SCTransform"){
    #Co-expression network analysis
    single_multiple_sample_clustering <- SetDatExpr(
      single_multiple_sample_clustering,
      group_name = cluster_number, # the name of the group of interest in the group.by column
      group.by=cluster_types, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
      assay = 'SCT', # using RNA assay
      slot = 'data' # using normalized data
    )
  }
  else if (index_s_hdwgcna1 == "multiple_sample_subclustering" & index_subclustering_multiple_sample_normalization_method_hdwgcna == "LogNormalize"){
    #Co-expression network analysis
    single_multiple_sample_clustering <- SetDatExpr(
      single_multiple_sample_clustering,
      group_name = cluster_number, # the name of the group of interest in the group.by column
      group.by=cluster_types, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
      assay = 'RNA', # using RNA assay
      slot = 'data' # using normalized data
    )
  }
  else if(index_s_hdwgcna1 == "multiple_sample_subclustering" & index_subclustering_multiple_sample_normalization_method_hdwgcna == "SCTransform"){
    #Co-expression network analysis
    single_multiple_sample_clustering <- SetDatExpr(
      single_multiple_sample_clustering,
      group_name = cluster_number, # the name of the group of interest in the group.by column
      group.by=cluster_types, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
      assay = 'SCT', # using RNA assay
      slot = 'data' # using normalized data
    )
  }
  
  
  #Select soft-power threshold
  # Test different soft powers:
  single_multiple_sample_clustering <- TestSoftPowers(
    single_multiple_sample_clustering,
    networkType = index_s_hdwgcna9 # you can also use "unsigned" or "signed hybrid"
  )
  # plot the results:
  plot_list <- PlotSoftPowers(single_multiple_sample_clustering)
  # assemble with patchwork
  plots802 <-wrap_plots(plot_list, ncol=2)
  power_table <- GetPowerTable(single_multiple_sample_clustering)
  #head(power_table)
  # construct co-expression network:
  single_multiple_sample_clustering <- ConstructNetwork(
    single_multiple_sample_clustering,
    tom_name = cluster_number, # name of the topoligical overlap matrix written to disk
    overwrite_tom = TRUE
  )
  
  output_dir <- paste0(getwd(),"/www/")  # Replace with your desired directory
  output_file <- "PlotDendrogram.pdf"          # File name for the PDF
  pdf_path <- file.path(output_dir, output_file)
  pdf(pdf_path, width = 8, height = 6)  # Specify width and height if needed
  PlotDendrogram(single_multiple_sample_clustering, main='hdWGCNA Dendrogram')
  dev.off()
  #Optional: inspect the topoligcal overlap matrix (TOM)
  TOM <- GetTOM(single_multiple_sample_clustering)
  
  
  #Module Eigengenes and Connectivity
  # compute all MEs in the full single-cell dataset
  single_multiple_sample_clustering <- ModuleEigengenes(
    single_multiple_sample_clustering,
    scale.model.use= index_s_hdwgcna10, #"poisson", or "negbinom"
    #merge.cut.height = 0.25
    #group.by.vars="condition"
  )
  # harmonized module eigengenes:
  hMEs <- GetMEs(single_multiple_sample_clustering, harmonized=index_s_hdwgcna11)
  # module eigengenes:
  MEs <- GetMEs(single_multiple_sample_clustering, harmonized=FALSE)
  
  # compute eigengene-based connectivity (kME):
  single_multiple_sample_clustering <- ModuleConnectivity(
    single_multiple_sample_clustering,
    group.by = cluster_types, group_name = cluster_number
  )
  # rename the modules
  single_multiple_sample_clustering <- ResetModuleNames(
    single_multiple_sample_clustering,
    new_name = paste0("M-", cluster_number)
  )
  # plot genes ranked by kME for each module
  plots804 <-PlotKMEs(single_multiple_sample_clustering, ncol=3, n_hubs = index_s_hdwgcna13)
  # get the module assignment table:
  modules <- GetModules(single_multiple_sample_clustering) %>% subset(module != 'grey')
  mods <- levels(modules$module); mods <- mods[mods != 'grey']
  # get hub genes
  hub_df <- GetHubGenes(single_multiple_sample_clustering, n_hubs = index_s_hdwgcna12)
    
  # make a featureplot of hMEs for each module
  plot_list <- ModuleFeaturePlot(
    single_multiple_sample_clustering,
    features='hMEs', # plot the hMEs
    order=TRUE # order so the points with highest hMEs are on top
  )
  
  # stitch together with patchwork
  plots805 <-wrap_plots(plot_list, ncol=3)
  
  #correlation plot
  output_file1 <- "PlotModuleCorrelogram.pdf"          # File name for the PDF
  pdf_path1 <- file.path(output_dir, output_file1)
  pdf(pdf_path1, width = 8, height = 6) 
  ModuleCorrelogram(single_multiple_sample_clustering)
  dev.off()
  
  
  # add hMEs to Seurat meta-data:# get hMEs from seurat object
  MEs <- GetMEs(single_multiple_sample_clustering, harmonized=TRUE)
  modules1 <- GetModules(single_multiple_sample_clustering)
  mods <- levels(modules1$module); mods <- mods[mods != 'grey']
  
  # add hMEs to Seurat meta-data:
  single_multiple_sample_clustering@meta.data <- cbind(single_multiple_sample_clustering@meta.data, MEs)
  # plot with Seurat's DotPlot function
  plots807 <- DotPlot(single_multiple_sample_clustering, features=mods, group.by = cluster_types) +
    RotatedAxis() +
    scale_color_gradient2(high='red', mid='grey95', low='blue')
  
  ModuleNetworkPlot(
    single_multiple_sample_clustering
  )
  
  tempdir <-getwd()
  pdf_dir <- paste0(getwd(),"/ModuleNetworks")  # Replace with your directory
  
  # List all PDF files in the directory
  pdf_files <- list.files(pdf_dir, pattern = "\\.pdf$", full.names = TRUE)
  output_pdf <- "www/combined_output.pdf"
  # Combine the PDF files
  pdf_combine(pdf_files, output = output_pdf)
  # Remove the directory
  unlink(pdf_dir, recursive = TRUE)
  
  #Applying UMAP to co-expression networks
  single_multiple_sample_clustering <- RunModuleUMAP(
    single_multiple_sample_clustering,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors=15, # neighbors parameter for UMAP
    min_dist=0.1 # min distance between points in UMAP space
  )
  
  output_file2 <- "ModuleUMAPPlot.pdf"          # File name for the PDF
  pdf_path2 <- file.path(output_dir, output_file2)
  pdf(pdf_path2, width = 8, height = 6) 
  ModuleUMAPPlot(
    single_multiple_sample_clustering,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=index_s_hdwgcna13,# how many hub genes to plot per module?
    keep_grey_edges=index_s_hdwgcna14
  )
  dev.off()
  
  tom_dir <- paste0(getwd(),"/TOM")
  unlink(tom_dir, recursive = TRUE)
  
  return(list(plot801 = plots801, plot802 = plots802, plot804 = plots804, plot805 = plots805, plot807 = plots807, text_summary = tempdir, data1 = power_table, data2 = modules, data3 = hub_df, data4 = single_multiple_sample_clustering))
}
