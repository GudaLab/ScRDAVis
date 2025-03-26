datainput_subclustering_multiple_sample <- function(index_subclustering_multiple_sample_file, index_subclustering_multiple_sample_celltype, index_m_subclustering_4, index_m_subclustering_5, index_m_subclustering1, index_m_subclustering2, index_m_subclustering3){
  multiple_sample_clustering  <- index_subclustering_multiple_sample_file
  
  if (index_m_subclustering1 == "seurat_clusters"){
  Idents(multiple_sample_clustering) <- "seurat_clusters"
  subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering2)
  }
  else if (index_m_subclustering1 == "predicted"){
	if (index_subclustering_multiple_sample_celltype == "sctype_classification"){
    Idents(multiple_sample_clustering) <- "sctype_classification"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
   else if (index_subclustering_multiple_sample_celltype == "singleR_labels"){
    Idents(multiple_sample_clustering) <- "singleR_labels"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
   else if (index_subclustering_multiple_sample_celltype == "GPTCelltype"){
    Idents(multiple_sample_clustering) <- "GPTCelltype"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
   else if (index_subclustering_multiple_sample_celltype == "cell_type"){
    Idents(multiple_sample_clustering) <- "cell_type"
    subclustering_multiple_sample <-subset(multiple_sample_clustering, idents = index_m_subclustering3)
  }
  }
  else if (index_m_subclustering1 == "selected_gene"){
   #Idents(multiple_sample_clustering) <- "orig.ident"
   #gene_of_interest <- index_m_subclustering_4
   #selected_cells <- WhichCells(multiple_sample_clustering, expression =  get(gene_of_interest) > 1)
   #subclustering_multiple_sample <-subset(multiple_sample_clustering, cells = selected_cells)
   #Idents(multiple_sample_clustering) <- "orig.ident"  
   gene_of_interest <- unlist(strsplit(index_m_subclustering_4, ","))
#gene_of_interest <- index_m_subclustering_4
threshold <- 1

# First, verify genes exist
valid_genes <- intersect(gene_of_interest, rownames(multiple_sample_clustering))

if(length(valid_genes) == 0){
  stop("None of the genes in gene_of_interest are present in your Seurat object.")
}

# Access normalized data properly using GetAssayData
expr_data <- GetAssayData(multiple_sample_clustering, assay = "RNA", slot = "data")

# Select cells (OR condition: at least one gene above threshold)
selected_cells_or <- colnames(expr_data)[
  colSums(expr_data[valid_genes, , drop = FALSE] > threshold) == length(valid_genes)
]

# Create subset Seurat object
subclustering_multiple_sample <- subset(multiple_sample_clustering, cells = selected_cells_or)
  }
 else if (index_m_subclustering1 == "exclude_selected_gene"){
    
    gene_of_interest <- unlist(strsplit(index_m_subclustering_5, ","))
    threshold <- 1
    
    valid_genes <- intersect(gene_of_interest, rownames(multiple_sample_clustering))

if(length(valid_genes) == 0){
  stop("None of the genes in gene_of_interest are present in your Seurat object.")
}

# Access normalized data properly using GetAssayData
expr_data <- GetAssayData(multiple_sample_clustering, assay = "RNA", slot = "data")

# Select cells (OR condition: at least one gene above threshold)
selected_cells_or <- colnames(expr_data)[
  colSums(expr_data[valid_genes, , drop = FALSE] > threshold) == length(valid_genes)
]

# Create subset Seurat object
subclustering_multiple_sample <- subset(multiple_sample_clustering, cells = setdiff(colnames(expr_data), selected_cells_or))
  }
  table1 <- table(subclustering_multiple_sample$orig.ident) %>% as.data.frame 
  colnames(table1) <- c("Sample names", "Cell counts")
  
  subclustering_multiple_sample[["percent.mt"]] <- PercentageFeatureSet(subclustering_multiple_sample, pattern = "^MT-")


  plots1 <- VlnPlot(subclustering_multiple_sample, features = "nFeature_RNA", ncol = 1)
  plots2 <- VlnPlot(subclustering_multiple_sample, features = "nCount_RNA", ncol = 1)
  plots3 <- VlnPlot(subclustering_multiple_sample, features = "percent.mt", ncol = 1)
  
  return(list(plot = plots1+plots2+plots3, data1 = table1, data2 = subclustering_multiple_sample))
}