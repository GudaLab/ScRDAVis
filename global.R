## Check packages and install##
if (!require("shiny")) install.packages("shiny", dependencies = TRUE)
if (!require("remotes")) install.packages("remotes", dependencies = TRUE)
if (!require("DT")) install.packages("DT")
if (!require("shinythemes")) install.packages("shinythemes")
if (!require("shinyjs")) install.packages("shinyjs")
if (!require("shinyFiles")) install.packages("shinyFiles")
if (!require("shinyWidgets")) install.packages("shinyWidgets")
if (!require("shinycssloaders")) install.packages("shinycssloaders")
if (!require('devtools')) install.packages("devtools")
if (!require("ggplot2")) install.packages("ggplot2", dependencies = TRUE)
if (!require('ggplotify')) install.packages("ggplotify")
if (!require("data.table")) install.packages("data.table")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("shinydashboard")) install.packages("shinydashboard")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tibble")) install.packages("tibble") 
if (!require("HGNChelper")) install.packages("HGNChelper") 
if (!require("openai")) install.packages("openai")
if (!require("metap")) install.packages("metap")
if (!require("harmony")) install.packages("harmony")
if (!require("ggrepel")) install.packages("ggrepel")
if (!require("R.utils")) install.packages("R.utils")
if (!require("circlize")) install.packages("circlize")
if (!require("hdf5r")) install.packages("hdf5r")
if (!require("ggupset")) install.packages("ggupset")
if (!require("gridExtra")) install.packages("gridExtra")
if (!require("ggalluvial")) install.packages("ggalluvial")
if (!require("NMF")) install.packages("NMF")
if (!require("ggraph")) install.packages("ggraph")
if (!require("igraph")) install.packages("igraph")
if (!require("cowplot"))install.packages("cowplot")
if (!require("pdftools"))install.packages("pdftools")
if (!require("xgboost")) remotes::install_version("xgboost", version = "1.7.11.1", repos = "https://cloud.r-project.org")
if (!require("Seurat")) install.packages("Seurat")
if (!require("msigdbr"))install.packages("msigdbr")
if (!require('filelock')) install.packages("filelock")
#if (!require("msigdf"))install.packages("msigdf", repos = "https://cloud.r-project.org")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("patchwork")) install.packages("patchwork")
if (!require("SeuratObject")) install.packages("SeuratObject", update = FALSE)
if (!require("BiocManager")) install.packages("BiocManager", update = FALSE)
if (!require("sctransform")) BiocManager::install("sctransform", update = FALSE)
if (!require("celldex")) BiocManager::install("celldex", update = FALSE)
if (!require("SingleR")) BiocManager::install("SingleR", update = FALSE)
if (!require("scRNAseq")) BiocManager::install("scRNAseq", update = FALSE)
if (!require("GenomicRanges"))BiocManager::install("GenomicRanges", update = FALSE)
if (!require("glmGamPoi"))BiocManager::install("glmGamPoi", update = FALSE)
if (!require("scran"))BiocManager::install("scran", update = FALSE)
if (!require("EnhancedVolcano"))BiocManager::install("EnhancedVolcano", update = FALSE)
if (!require("ComplexHeatmap")) BiocManager::install("ComplexHeatmap", update = FALSE)
if (!require("clusterProfiler"))BiocManager::install("clusterProfiler", update = FALSE)
if (!require("org.Hs.eg.db"))BiocManager::install("org.Hs.eg.db", update = FALSE)
if (!require("org.Mm.eg.db"))BiocManager::install("org.Mm.eg.db", update = FALSE)
if (!require("org.Mmu.eg.db"))BiocManager::install("org.Mmu.eg.db", update = FALSE)
if (!require("org.Rn.eg.db"))BiocManager::install("org.Rn.eg.db", update = FALSE)
if (!require("org.Ss.eg.db"))BiocManager::install("org.Ss.eg.db", update = FALSE)
if (!require("ReactomePA"))BiocManager::install("ReactomePA", update = FALSE)
if (!require("msigdbr"))BiocManager::install("msigdbr", update = FALSE)
if (!require("fgsea"))BiocManager::install("fgsea", update = FALSE)
if (!require("enrichplot"))BiocManager::install("enrichplot", update = FALSE)
#if (!require("CellChat"))BiocManager::install("sqjin/CellChat", update = FALSE)
if (!require("multtest"))BiocManager::install("multtest", update = FALSE)
if (!require("WGCNA"))BiocManager::install("WGCNA", update = FALSE)
# if (!require("motifmatchr"))BiocManager::install("motifmatchr", update = FALSE)
# if (!require("TFBSTools"))BiocManager::install("TFBSTools", update = FALSE)
if (!require("JASPAR2020"))BiocManager::install("JASPAR2020", update = FALSE)
if (!require("JASPAR2024"))BiocManager::install("JASPAR2024", update = FALSE)
if (!require("EnsDb.Hsapiens.v86"))BiocManager::install("EnsDb.Hsapiens.v86", update = FALSE)
if (!require("EnsDb.Mmusculus.v79"))BiocManager::install("EnsDb.Mmusculus.v79", update = FALSE)
if (!require("BSgenome.Hsapiens.UCSC.hg38"))BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", update = FALSE)
if (!require("BSgenome.Mmusculus.UCSC.mm10"))BiocManager::install("BSgenome.Mmusculus.UCSC.mm10", update = FALSE)
source("scripts/PrctCellExpringGene.R")
options(shiny.maxRequestSize=2000*1024^2)
options(future.globals.maxSize= 925289600000)
Sys.setenv(OPENAI_API_KEY = '')  #Add your key here

optional_github_pkgs <- c(
  "DoubletFinder",
  "GPTCelltype",
  "presto",
  "monocle3",
  "SeuratWrappers",
  "SeuratDisk",
  "CellChat",
  "genesorteR",
  "enrichR",
  "hdWGCNA"
)

for (pkg in optional_github_pkgs) {
  suppressWarnings(suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE)))
}

if (!nzchar(Sys.getenv("OPENAI_API_KEY"))) {
  message("OPENAI_API_KEY is not set. GPTCelltype will remain unavailable until you set it in the R session or environment.")
}


# Determine OS-specific cache path
if (.Platform$OS.type == "windows") {
  # For Windows
  cache_path <- file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache", "R", "BiocFileCache")
} else {
  # For Linux/macOS
  cache_path <- file.path(Sys.getenv("HOME"), ".cache", "R", "BiocFileCache")
}

# Create directory if it doesn't exist
if (!dir.exists(cache_path)) {
  dir.create(cache_path, recursive = TRUE, showWarnings = FALSE)
}

# Set BiocFileCache directory environment variable
Sys.setenv("BIOCFILECACHE_DIR" = cache_path)

#views
count_file <- "view_counter.rds"
lock_file  <- "view_counter.lock"
if (!file.exists(count_file)) saveRDS(0L, count_file)

read_count <- function() {
  readRDS(count_file)
}
increment_count <- function() {
  lock <- lock(lock_file, timeout = 5000)   # wait up to 5s for the lock
  on.exit(unlock(lock), add = TRUE)
  n <- readRDS(count_file)
  n <- n + 1L
  saveRDS(n, count_file)
  n
}


#Inactivity timer
# 
# inactivity <- sprintf("function idleTimer() {
# var t = setTimeout(logout, %s);
# window.onmousemove = resetTimer; // catches mouse movements
# window.onmousedown = resetTimer; // catches mouse movements
# window.onclick = resetTimer;     // catches mouse clicks
# window.onscroll = resetTimer;    // catches scrolling
# window.onkeypress = resetTimer;  //catches keyboard actions
# 
# function logout() {
# Shiny.setInputValue('timeOut', '%ss')
# }
# 
# function resetTimer() {
# clearTimeout(t);
# t = setTimeout(logout, %s);  // time is in milliseconds (1000 is 1 second)
# }
# }
# idleTimer();", 300*1000, 300, 300*1000)
