# ScRDAVis

ScRDAVis is a browser-based, user-friendly R Shiny application designed for researchers without programming expertise to analyze and visualize single-cell RNA (scRNA) data. It supports single and multiple sample analyses as well as group comparisons, offering a range of functionalities for comprehensive data exploration.

## Features

ScRDAVis includes nine functional modules:
1. **Single or Multiple Samples Analysis**
   - Stats
   - Sample Groups and QC Filtering
   - Normalization and PCA Analysis
   - Clustering
   - Remove Doublets
   - Marker Identification
   - Cell Type Prediction
   - Cluster-Based Plots
   - Condition-Based Analysis
2. **Subclustering**
3. **Correlation Network Analysis**
4. **Genome Ontology (GO) Terms**
5. **Pathway Analysis**
6. **GSEA Analysis**
7. **Cell-Cell Communication**
8. **Trajectory and Pseudotime Analysis**
9. **Co-Expression and TF Analysis**
   - Co-Expression Network Analysis
   - Transcription Factor Regulatory Network Analysis

## Use ScRDAVis Online
ScRDAVis is deployed online and accessible at:  
**[https://www.gudalab-rtools.net/ScRDAVis](https://www.gudalab-rtools.net/ScRDAVis)**

## Launch ScRDAVis Locally

### Prerequisites
Ensure the following software is installed:
- **R** (>= 4.4.2): [Download R](https://www.r-project.org/)
- **RStudio** (>= 2024.09.1): [Download RStudio](https://posit.co/download/rstudio-desktop/)
- **Bioconductor** (>= 3.19)
- **Shiny** (>= 1.9.1)

**Note:** ScRDAVis has been tested with these versions. Using older versions of R may cause errors during package installation. Updating to the latest R version is recommended.

### Installation

Run the following commands in an R session to install required packages:

```R
if (!require("BiocManager")) install.packages("BiocManager", update = FALSE)
if (!require("devtools")) install.packages("devtools", update = FALSE)

# Function to check and install CRAN packages
install_cran_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
    }
  }
}

# Function to check and install Bioconductor packages
install_bioc_packages <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      BiocManager::install(pkg, update = FALSE, dependencies = TRUE)
    }
  }
}

# Function to check and install GitHub packages
install_github_packages <- function(packages) {
  for (pkg in names(packages)) {
    if (!require(pkg, character.only = TRUE)) {
      devtools::install_github(packages[[pkg]], upgrade = "never", dependencies = TRUE)
    }
  }
}

# List of CRAN packages
cran_packages <- c("shiny", "DT", "shinythemes", "shinyjs", "shinyFiles", "shinyWidgets", "shinycssloaders", "ggplot2", "data.table", "ggpubr", "shinydashboard", "dplyr", "tibble", "HGNChelper", "openai", "metap", "harmony", "ggrepel", "R.utils", "circlize", "hdf5r", "ggupset", "gridExtra", "ggalluvial", "NMF", "ggraph", "igraph", "cowplot", "pdftools", "xgboost")

# List of Bioconductor packages
bioc_packages <- c("Seurat", "SeuratObject", "sctransform", "celldex", "SingleR", "scRNAseq", "glmGamPoi", "scran", "EnhancedVolcano", "ComplexHeatmap", "clusterProfiler", "org.Hs.eg.db", "org.Mm.eg.db", "org.Mmu.eg.db", "org.Rn.eg.db", "org.Ss.eg.db", "ReactomePA", "msigdbr", "fgsea", "enrichplot", "multtest", "WGCNA", "hdWGCNA", "motifmatchr", "TFBSTools", "GenomicRanges", "JASPAR2020", "EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10")

# List of GitHub packages
github_packages <- list(
  "DoubletFinder" = "chris-mcginnis-ucsf/DoubletFinder",
  "GPTCelltype" = "Winnie09/GPTCelltype",
  "openxlsx" = "ycphs/openxlsx",
  "presto" = "immunogenomics/presto",
  "monocle3" = "cole-trapnell-lab/monocle3",
  "SeuratWrappers" = "satijalab/seurat-wrappers",
  "SeuratDisk" = "mojaveazure/seurat-disk",
  "patchwork" = "thomasp85/patchwork",
  "CellChat" = "sqjin/CellChat",
  "genesorteR" = "mahmoudibrahim/genesorteR",
  "enrichR" = "wjawaid/enrichR",
  "hdWGCNA" = "smorabit/hdWGCNA"
)

# Install all packages
install_cran_packages(cran_packages)
install_bioc_packages(bioc_packages)
install_github_packages(github_packages)

###Start the App

    Open an R session in RStudio.
    Run the following commands to start the application:

library(shiny)
shiny::runGitHub('ScRDAVis', 'GudaLab')

Alternatively, download the source code from GitHub and run:

library(shiny)
runApp('/path/to/the/ScRDAVis-master', launch.browser = TRUE)

###Usage

A detailed user manual is available under the "Manual" tab at:
https://www.gudalab-rtools.net/ScRDAVis
###Tested Platforms

This application was tested on:

    Linux (Red Hat and Ubuntu)
    Windows (10 and 11)
