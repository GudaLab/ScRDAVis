datainput_multiple_sample <- function(index_multiple_sample_file, index_multiple_sample_file1, index_multiple_sample_file_names, index_multiple_sample_format, index_multiple_sample_name, index_multiple_sample_cell, index_multiple_sample_genes){
index_multiple_sample_format <- as.character(index_multiple_sample_format)

if (index_multiple_sample_format == "h5") {
    #remove extension
    #files <- list.files(pattern="*.h5")
    h5_files1 <- sub('\\.h5$', '', index_multiple_sample_file_names) 
    h5_read <- lapply(index_multiple_sample_file, Read10X_h5)
    
    
    #list file names
    names(h5_read) <- h5_files1
    multiple_list <- mapply(CreateSeuratObject, counts=h5_read, project=names(h5_read), min.cells = index_multiple_sample_cell, min.features = index_multiple_sample_genes)
    multiple_list_merged <- merge(multiple_list[[1]], y = multiple_list[-1], add.cell.ids = names(multiple_list),project="merged")
    
    Idents(multiple_list_merged) <- multiple_list_merged@meta.data$orig.ident
    table1 <- table(multiple_list_merged$orig.ident) %>% as.data.frame 
    colnames(table1) <- c("Sample names", "Cell counts")
  }

else if (index_multiple_sample_format == "MFB"){
#index_multiple_sample_file1 <- list.files(pattern="*.gz")
h5_files1 <- sub('\\_barcodes.tsv.gz$', '', index_multiple_sample_file_names) 
h5_files1 <- sub('\\_features.tsv.gz$', '', h5_files1) 
h5_files1 <- sub('\\_matrix.mtx.gz$', '', h5_files1) 
h5_files1 <- unique(h5_files1)
#Create a list of count matrices
multiple_list <- list()

path <- getwd()
setwd(index_multiple_sample_file1)
for(x in h5_files1){
  name <- gsub('','', x)
  h5_read <- ReadMtx(mtx = paste0(x,'_matrix.mtx.gz'), features = paste0(x,'_features.tsv.gz'),cells = paste0(x,'_barcodes.tsv.gz'))
  # create seurat objects
  test <- assign(name, CreateSeuratObject(counts =  h5_read, project = name, min.cells = index_multiple_sample_cell, min.features = index_multiple_sample_genes))
  multiple_list[x] =  test 
}
multiple_list_merged <- merge(multiple_list[[1]], y = multiple_list[-1], add.cell.ids = names(multiple_list),project="merged")

Idents(multiple_list_merged) <- multiple_list_merged@meta.data$orig.ident
table1 <- table(multiple_list_merged$orig.ident) %>% as.data.frame 
colnames(table1) <- c("Sample names", "Cell counts")
setwd(path)

}
else if (index_multiple_sample_format == "seurat_object")
{
  multiple_list_merged <- readRDS(index_multiple_sample_file)
  table1 <- table(multiple_list_merged$orig.ident) %>% as.data.frame 
  colnames(table1) <- c("Sample names", "Cell counts")
  multiple_list <- SplitObject(multiple_list_merged, split.by = "orig.ident")
}

else if (index_multiple_sample_format == "matrix_count")
{
  multiple_list_merged <- read.table(index_multiple_sample_file, header = TRUE, row.names = 1, sep = "\t")
  multiple_list_merged <- CreateSeuratObject(counts = multiple_list_merged, min.cells = index_multiple_sample_cell, min.features = index_multiple_sample_genes)
  table1 <- table(multiple_list_merged$orig.ident) %>% as.data.frame 
  colnames(table1) <- c("Sample names", "Cell counts")
  multiple_list <- SplitObject(multiple_list_merged, split.by = "orig.ident")
}

else if (index_multiple_sample_format == "exampledata")
{
  index_multiple_sample_file<- "www/example_data/example_data.RDS"
  index_multiple_sample_file_names <- "Example data to test the tool (C1_vs_P1 from GSE266873)"
  multiple_list_merged <- readRDS(index_multiple_sample_file)
  table1 <- table(multiple_list_merged$orig.ident) %>% as.data.frame 
  colnames(table1) <- c("Sample names", "Cell counts")
  multiple_list <- SplitObject(multiple_list_merged, split.by = "orig.ident")
}

multiple_list_merged[["percent.mt"]] <- PercentageFeatureSet(multiple_list_merged, pattern = "^MT-")

plots1 <- VlnPlot(multiple_list_merged, features = "nFeature_RNA", ncol = 1)
plots2 <- VlnPlot(multiple_list_merged, features = "nCount_RNA", ncol = 1)
plots3 <- VlnPlot(multiple_list_merged, features = "percent.mt", ncol = 1)
plots4 <- FeatureScatter(multiple_list_merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plots5 <- FeatureScatter(multiple_list_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

return(list(plot1 = plots1 + plots2 + plots3, data1 = table1, Plot2 = plots4 + plots5, text_summary = index_multiple_sample_file_names, data2 = multiple_list_merged, data3 = table1[,1], data4 = multiple_list))
}

