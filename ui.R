library(shiny)
library(shinythemes)
library(shinyFiles)
library(DT)
library(shinyjs)


shinythemes::themeSelector()
shinyUI(
    navbarPage(id ="menu_tabs",
    theme = shinytheme("cerulean"),
    "",
    tabPanel(
      "ScRDAVis",
      mainPanel(
        h1("Single Cell RNA Data Analysis and Visualization (ScRDAVis)",align = "center"),
        hr(),
        h3("Introduction"),
        hr(),
    HTML("<p>ScRDAVis is a browser-based and user-friendly R Shiny application designed for researchers without programming proficiency to analyze and visualize single-cell RNA (scRNA) results. It supports single and multiple sample analyses as well as group comparisons. The application includes the following key functional analyses:</p>

    <h3>1. Single or Multiple Samples Analysis</h3>
    <p>This section offers various tabs to analyze one or more samples, which can be grouped into up to six groups.</p>
    
    <h5>1.1 Stats</h5>
    <p>Displays the QC plot and cell summary of the uploaded sample(s).</p>

    <h5>1.2 Sample Groups and QC Filtering</h5>
    <p>Assists in filtering QC metrics for the sample(s) for further analysis.</p>

    <h5>1.3 Normalization and PCA Analysis</h5>
    <p>Allows normalization of samples using multiple methods and generates PCA plots.</p>

    <h5>1.4 Clustering</h5>
    <p>Uses the Seurat clustering algorithm to group cells into clusters and visualizes them with UMAP or tSNE.</p>

    <h5>1.5 Remove Doublets</h5>
    <p>Employs DoubletFinder to detect doublet or singlet cells, allowing users to keep or remove doublets cells.</p>

    <h5>1.6 Marker Identification</h5>
    <p>Identifies markers for all clusters, a specific cluster, or between clusters and supports the identification of conserved markers.</p>

    <h5>1.7 Cell Type Prediction</h5>
    <p>Offers multiple options for cell type identification, including ScType, SingleR, GPTCelltype, or custom user-provided labels.</p>

    <h5>1.8 Cluster-Based Plots</h5>
    <p>Displays expressed genes in each cluster using Dot, Violin, Ridge, or Feature plots.</p>

    <h5>1.9 Condition-Based Analysis</h5>
    <p>Identifies expressed genes between two groups, with visualization options including Dot, Violin, Ridge, Feature, or Volcano plots.</p>

    <h3>2. Subclustering</h3>
    <p>Allows sub-clustering within one or more clusters from single or multiple sample analyses, following similar steps as in the primary analysis.</p>

    <h3>3. Correlation Network Analysis</h3>
    <p>Uses the genesorteR package to identify the correlation between cell clusters. Provides correlation summary tables and visualizations of correlation matrix and network plots.</p>

    <h3>4. Genome Ontology (GO) Terms</h3>
    <p>Uses the clusterProfiler package to identify biological processes, molecular functions, and cellular components for marker genes. Provides GO summary tables and visualizations in Dot, Bar, Net and UpSetplots.</p>

    <h3>5. Pathway Analysis</h3>
    <p>Employs the clusterProfiler and ReactomePA packages to identify pathways in single or multiple clusters, with results displayed in Dot, Bar, Net and UpSetplots.</p>

    <h3>6. GSEA Analysis</h3>
    <p>Performs Gene Set Enrichment Analysis (GSEA) using the fgsea and msigdb packages to identify enriched gene sets. Results are displayed in GSEA plots, Bar plots, and PlotGseaTables.</p>

    <h3>7. Cell-Cell Communication</h3>
    <p>Uses the Cellchat package to identify signaling communication between clusters, with receptor-ligand interactions visualized in Circular, Chord, Heatmap, Bubble, Bar, and Violin plots.</p>

    <h3>8. Trajectory and Pseudotime Analysis</h3>
    <p>Utilizes the Monocle3 package to order clusters in pseudotime and analyze gene function changes over time. Visualizations include trajectory and pseudotime plots, bar plots, and gene functional changes in pseudotime.</p>
    
    <h3>9. Co-Expression and TF analysis</h3>
    <h5>9.1 Co-Expression Network Analysis</h5>
    <p>Uses the hdWGCNA package to identify co-expression networks as undirected, weighted gene networks. These are visualized through co-expression networks with modules, soft power plots, module relationship plots, module network plots, and module UMAP plots</p>
    <h5>9.2 Transcription factor regulatory network analysis</h5>
    <p>Uses the hdWGCNA package to identify the transcription factor (TFs) within co-expression modules. These TFs play a key role in regulating gene expression networks in single-cell data. These TFs are visualized through bar plot, network plot and module UMAP plots</p>
    
    <h3>Outputs and Visualization</h3>
    <p>ScRDAVis provides publication-quality plots in seven formats: JPG, TIFF, PDF, SVG, BMP, EPS, and PS. Summary tables are also generated in .csv format for easy visualization and download.</p>
        <hr>
        <h3> use ScRDAVis online</h3>
        <p>ScRDAVis is deployed at: <a href='https://www.gudalab-rtools.net/ScRDAVis'>https://www.gudalab-rtools.net/ScRDAVis</a></p>
        <hr>
        <h3> Launch ScRDAVis using R and GitHub: </h3>
        <p> ScRDAVis were deposited under the GitHub repository: <a href='https://www.gudalab-rtools.net/ScRDAVis'>https://www.gudalab-rtools.net/ScRDAVis</a><br>
        Before running the app, users must have the following versions installed: R (>= 4.4.2), RStudio (>= 2024.12.0), Bioconductor (>= 3.20) and Shiny (>= 1.10.0) (Tested with this version).<br>
        Note: ScRDAVis has been tested with these versions. If users are running an older version of R, they may encounter errors during package installation. Therefore, it is recommended to update R to the latest version first.<br>
       Once R is open in the command line or in RStudio, users should run the following command in R to install the shiny package.<br><br></p>
          
<pre>install.packages('shiny')<br>
library(shiny)</pre>
          <hr>
          <h3>Start the app</h3>
          Start the R session using RStudio and run these lines:<br><br>
<pre>shiny::runGitHub('ScRDAVis','GudaLab')</pre>
or
Alternatively, download the source code from GitHub and run the following command in the R session using RStudio:
<pre>
library(shiny)
runApp('/path/to/the/ScRDAVis-master', launch.browser=TRUE)</pre>
<hr>
<h3> Developed and maintained by</h3>
<p>ScRDAVis was developed by Sankarasubramanian Jagadesan and Babu Guda. We share a passion for developing a user-friendly tool for biologists, particularly those who do not have access to bioinformaticians or programming expertise.
</p>
<hr>
"),        
      ),
    ),



##################################################################################################################################### 
##############################################################multiple samples#######################################################
#####################################################################################################################################
 tabPanel(
 "Single or Multiple Samples analysis",
 useShinyjs(),

 sidebarLayout(
   sidebarPanel(id="multiple_sidebar",
     selectInput("multiple_sample_format", label = "Select Input format", choices = list("Cell ranger h5 format" = "h5", "Cellranger Matrix, Feauture, Barcodes files" ="MFB", "Seurat object" = "seurat_object", "Matrix count file" = "matrix_count" ), selected = "h5"),
     
     #h3("Upload files"),
     fileInput("multiple_sample_file", label = "Upload multiple file(s) at once (eg: FILENAME1.h5 or FILENAME_1.h5)", multiple = T, accept =".h5"),
     fileInput("multiple_sample_file_mfb", label = "Upload multiple files at once (eg: FILENAME_matrix.mtx.gz, FILENAME_features.tsv, FILENAME_barcodes.tsv.gz)", multiple = T, accept =".gz"),
     fileInput("multiple_sample_file_rds", label = "Upload single seurat object file with single or multiple samples", multiple = F, accept=c(".rds", "RDS")),
     fileInput("multiple_sample_file_txt", label = "Upload single matrix file with single or multiple samples (eg: FILENAME1.txt or FILENAME1.tsv)", multiple = F, accept=".txt"),
     numericInput("multiple_sample_cell", label="Include minimum number of cells expressing in gene", min = 0, max = 1000000, value = "3"),
     numericInput("multiple_sample_genes", label="Include minimum number of genes a cell must express", min = 0, max = 1000000, value = "200"),
    actionBttn("multiple_sample_submit", "Submit",  style = "unite",color = "primary" ),
     width = 3),
   mainPanel(id="multiple_main_menu",
     tabsetPanel(
       type = "tabs", id="multiple_tabsets",
######################################################Tab1########################################################         
#####################multiple_samples QC_before filtering##########################
       
       tabPanel(
         "Stats",
         box(id = "m_bf_box1", width = 12,
         h3("Sample names of uploaded file(s)"),
         verbatimTextOutput("text_level"),
         ),
         fluidRow(
           box(id = "m_bf_box2", width = 10,
           column(
         h3("Number of cells in the given sample(s)"), 
         downloadBttn(outputId = "download_multiple_cell_table", label = "Download as csv"),
         withSpinner(dataTableOutput("multiple_cell_table")), width = 6),
           ),
           ),
       fluidRow(
           box(id = "m_bf_box3", 
               h3("QC Plot before filtering"), 
              actionBttn("download_m_qc_before_filtering", "Download plot",  style = "unite",color = "primary",icon = icon("download")), 
               status = "warning", 
               solidHeader = F,
               withSpinner(plotOutput("m_qc_before_filtering", width = "1000px", height = "1000px"))),
         ),
         fluidRow(
           box(id = "m_bf_box4", 
               h3("Feature-Feature relationships plot"), 
               actionBttn("download_m_ff_before_filtering", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
               status = "warning", 
               solidHeader = F,
               withSpinner(plotOutput("m_ff_before_filtering", width = "1000px", height = "500px"))),
         ),
       box(id = "m_bf_box5",  width = 12,
         hr(),
         downloadBttn(outputId = "m_so_before_filtering", label = "Download Seurat Object"),
         actionLink("link_m_qc_filtering", actionBttn("MNextButtonLoad", "NEXT STEP (Sample Groups and QC Filtering)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
         ),
       ),


######################################################Tab1.2########################################################  
#########################update groups##################################################
#####################multiple_samples QC_after filtering##########################       
       tabPanel(
         "Sample Groups and QC Filtering",
         selectInput("multiple_group_count", label="Select number of sample group(s)", multiple = F, choices = c("1"=1, "2" = 2, "3"=3, "4" = "4",  "5" = 5, "6" = "6" ), selected = 1),
         fluidRow(
           column(width =2, textInput("group1_name",label ="Type Group1 name", value = "Group1", width = NULL, placeholder = NULL)),
           column(width =2, uiOutput("group1_samples1")),
           column(width =2, textInput("group2_name",label ="Type Group2 name", value = "Group2", width = NULL, placeholder = NULL)),
           column(width =2, uiOutput("group2_samples2")),
         ), 
         fluidRow(
           column(width =2, textInput("group3_name",label ="Type Group3 name", value = "Group3", width = NULL, placeholder = NULL)),
           column(width =2, uiOutput("group3_samples3")),
           column(width =2, textInput("group4_name",label ="Type Group4 name", value = "Group4", width = NULL, placeholder = NULL)),
           column(width =2, uiOutput("group4_samples4")),
         ), 
         fluidRow(
           column(width =2, textInput("group5_name",label ="Type Group5 name", value = "Group5", width = NULL, placeholder = NULL)),
           column(width =2, uiOutput("group5_samples5")),
           column(width =2, textInput("group6_name",label ="Type Group6 name", value = "Group6", width = NULL, placeholder = NULL)),
           column(width =2, uiOutput("group6_samples6")),
         ),
         fluidRow(
           titlePanel("Define filtering parameters"),
           
           p("Exclude cells based on their number of expressed genes and the percentage of reads that map to the mitochondrial genome."),
           #   p("The parameters can be set on the right side of the plot and must be set using a higher value than the ones above. Otherwise, they will have no effect since the cells were already excluded."),
           
           
           column(width =3, numericInput("multiple_sample_min_count", label = "Keep cells that are expressed at least this number of genes", value = 0)),
           column(width =3, numericInput("multiple_sample_max_count", label = "Exclude cell that expressed more than this number of genes", value = 10000)),
           column(width =3, numericInput("multiple_sample_max_mito_perc", label="Filter cells that have more than this percentage mitochondrial counts",  value = 5)),
           br(),
           br(),
           column(width =3, actionBttn("multiple_sample_qc_filtering", "Update filtered data",  style = "unite",color = "primary", icon = icon("download"))),
           
         ),
         hr(),
         #withSpinner(verbatimTextOutput("multiple_cell_qc_text_level")),
         fluidRow(
           box(id = "m_qc_filter_box1",
           h3("Number of cells after QC"),
           column(
             h3("Sample(s) based"), 
             downloadBttn(outputId = "download_multiple_cell_table_after_qc", label = "Download as csv"),
             withSpinner(dataTableOutput("multiple_cell_table_after_qc")), width = 6),
           column(
             h3("Group(s) based"), 
             downloadBttn(outputId = "download_multiple_cell_table_after_qc2", label = "Download as csv"),
             withSpinner(dataTableOutput("multiple_cell_table_after_qc2")), width = 6),
         ),
         ),
         fluidRow(
           
         box(id = "m_qc_filter_box2", 
             h3("QC plot after filtering"),
             h3("Sample(s) based"),
         column(
         actionBttn("download_m_qc_after_filtering", "Download plot",  style = "unite",color = "primary", icon = icon("download")), 
             status = "warning", 
             solidHeader = F,
             withSpinner(plotOutput("m_qc_after_filtering", width = "900px", height = "1000px")), width = 7),
         ),
         
         box(id = "m_qc_filter_box3", 
             br(),
             br(),
             h3("Group(s) based"), 
             column(
             actionBttn("download_m_qc_after_filtering2", "Download plot",  style = "unite",color = "primary", icon = icon("download")), 
             status = "warning", 
             solidHeader = F,
             withSpinner(plotOutput("m_qc_after_filtering2", width = "500px", height = "1000px")), width = 5),
         ),
         ),
         fluidRow(
          box(id = "m_qc_filter_box4", 
                 column(
                  h3("Bar plots"),
                 h3("Sample(s) based"), 
                 actionBttn("download_m_qc_after_filtering3", "Download plot",  style = "unite",color = "primary", icon = icon("download")), 
                 status = "warning", 
                 solidHeader = F,
                 withSpinner(plotOutput("m_qc_after_filtering3", width = "900px", height = "800px")), width = 7),
          ),
           
             box(id = "m_qc_filter_box5", 
                 column(
                   br(),
                 h3("Group(s) based"), 
                 actionBttn("download_m_qc_after_filtering4", "Download plot",  style = "unite",color = "primary", icon = icon("download")), 
                 status = "warning", 
                 solidHeader = F,
                 withSpinner(plotOutput("m_qc_after_filtering4", width = "500px", height = "800px")),width = 5),
             ),
         ),
         fluidRow(
           box(id = "m_qc_filter_box6", width = 12,
           column(
             hr(),
         downloadBttn(outputId = "m_so_after_filtering", label = "Download Seurat Object"),
         actionLink("link_m_normalization", actionBttn("MNORMALIZATIONButtonLoad", "NEXT STEP (Normalization and PCA analysis)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
         width = 12),
         ),
       ),
       ),


#####################################################Tab1.3######################################   
#############################multiple_Normalization and PCA analysis#########################
tabPanel(
  "Normalization and PCA analysis",
  selectInput("multiple_sample_normalization_method", label = "Normalization method", choices = c("LogNormalize"="LogNormalize", "SCTransform"="SCTransform"), selected = "LogNormalize"),
  fluidRow( 
    #box(id = "boxLogNormalize",
    column(width =4, numericInput("multiple_sample_scale_factor", label = "Scale factor", value = 10000)),
    column(width =4, selectInput("multiple_sample_normalization_variable_genes", label = "variable genes detection", choices = c( "vst"="vst", "mean.var.plot (mvp)"="mean.var.plot", "dispersion (disp)"="dispersion"), selected = "vst")),  
    column(width =4, numericInput("multiple_sample_var_genes", label = "Number of top variable features", value = 2000)),
    #),
    column(width =4, numericInput("multiple_sample_pca_dim", label = "Number of dimensions (PCA)", value = 30)),
    br(),
    column(width =4, actionBttn("multiple_sample_normalization", "Submit",  style = "unite",color = "primary", icon = icon("download"))),
  ),
  hr(),
  
  fluidRow(
    box(id = "m_pca_box1", 
        h3("Dimension reduction heatmap  for PCA data"), 
        actionBttn("download_m_pca_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_pca_plot", height = 600))),
  
  box(id = "m_elbow_box", 
      h3("Elbow plot"), 
      actionBttn("download_m_elbow_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
      status = "warning", 
      solidHeader = F,
      withSpinner(plotOutput("m_elbow_plot", height = 600))),
  ),
  
  fluidRow(

    box(id = "m_pca_box2", 
        h3("PCA sample(s) based"), 
        actionBttn("download_m_pca2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_pca2_plot", height = 600))),
    
    box(id = "m_pca_box3", 
        h3("PCA group(s) based"), 
        actionBttn("download_m_pca3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_pca3_plot", height = 600))),
  ),
  box(id = "m_pca_box4", width = 12,
  hr(),
  downloadBttn(outputId = "m_normalization", label = "Download Seurat Object"),
  actionLink("link_m_clustering", actionBttn("MCLUSTERINGButtonLoad", "NEXT STEP (Clustering)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
), 
), 

################################Tab1.4######################################   
#############################multiple_Clustering#########################
tabPanel(
  "Clustering",
  fluidRow(
    h3("Nearest-neighbour graph construction"),
    column(4, numericInput("m_clustering1", label = "Number of dimensions", min = 0, max = 100, value = 30, step = 1)),
    column(4, numericInput("m_clustering2", label = "k.param", value = 20, step = 1)),
    column(4, numericInput("m_clustering3", label = "n.trees", value = 50, step = 1)),
  ),    
  fluidRow(		
    h3("Clustering parameters and integration method"),
    column(4, numericInput("m_clustering4", label = "Resolution", min = 0, max = 10, value = 0.5, step = 0.1)),
    column(4, selectInput("m_clustering5", label = "Clustering algorithm", choices = c("Louvain" = 1, "SLM" = 3, "Leiden" = 4), selected = 1)),
  # ),
  # fluidRow(
  #   h3("Select integration method"),
    column(4, selectInput("m_clustering13", label = "Integration method (except none, minimum two samples are required)", choices = c("None" = "None", "HarmonyIntegration" = "HarmonyIntegration", "CCAIntegration" = "CCAIntegration", "RPCAIntegration" = "RPCAIntegration", JointPCAIntegration = "JointPCAIntegration" ), selected = "None")),
  ),                 
  fluidRow( 
    h3("Dimension reduction"),
    column(3, radioButtons("m_clustering6", label = "Plot Options", choices = c("UMAP" = "umap", "t-SNE" = "tsne"), selected = "umap", inline=T)),
  ),                 
  fluidRow(   
    box(id = "m_umap_box",                      
        h3("UMAP parameters"), 
        column(3, numericInput("m_clustering7", label = "Number of dimensions", min = 0, max = 100, value = 30, step = 1)),
        column(3, numericInput("m_clustering8", label = "k-nearest-neighbours", min = 0, max = 50, value = 20, step = 1)),
        column(3, numericInput("m_clustering9", label = "min.dist", min = 0.001, max = 0.5, value = 0.3, step = 0.01)),
        column(3, selectInput("m_clustering10", label = "Show label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
    ),
  ),		
  fluidRow(                  
    box(id = "m_tsne_box",                      
        h3("t-SNE parameters"),                 
        column(3, numericInput("m_clustering11", label = "Number of dimensions",  min = 0, max = 100, value = 30, step = 1)),
        column(3, selectInput("m_clustering12", label = "Show label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
    ),
  ),   
  fluidRow(
    column(width =3, actionBttn("multiple_sample_clustering", "Submit",  style = "unite",color = "primary", icon = icon("download"))),
  ),
  hr(),
  
  fluidRow(
    box(id = "m_clustering_box1", 
        h3("UMAP / t-SNE cluster plot"), 
        actionBttn("download_m_umap_tsne1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_umap_tsne1_plot", height = 600))),
    
    box(id = "m_clustering_box2", 
        h3("Cluster based count bar plot"), 
        actionBttn("download_m_umap_tsne_bar1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_umap_tsne_bar1_plot", height = 600))),
  ),
  
  fluidRow(
    box(id = "m_clustering_box3", 
        h3("UMAP / t-SNE condition(s) based plot"), 
           actionBttn("download_m_umap_tsne2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_umap_tsne2_plot", height = 600))),
    
    box(id = "m_clustering_box4", 
        h3("Condition(s) based count bar plot"), 
        actionBttn("download_m_umap_tsne_bar2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_umap_tsne_bar2_plot", height = 600))),
  ),
  
  fluidRow(
    box(id = "m_clustering_box5", 
        h3("UMAP / t-SNE sample(s) based plot"), 
        actionBttn("download_m_umap_tsne3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_umap_tsne3_plot", height = 600))),
    
    box(id = "m_clustering_box6", 
        h3("Sample(s) based  count Bar plot"), 
        actionBttn("download_m_umap_tsne_bar3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_umap_tsne_bar3_plot", height = 600))),
  ),
  
  fluidRow(
    box(id = "m_clustering_box7", width = 5,
    column(
      h3("Number of cells in clusters", downloadBttn(outputId = "download_m_clustering_table1", label = "Download as csv")),
      withSpinner(dataTableOutput("m_clustering_table1")), width = 5),
    ),
    box(id = "m_clustering_box9", width = 7,
  column(
    h3("Number of cells in clusters based on condition(s)", downloadBttn(outputId = "download_m_clustering_table2", label = "Download as csv")),
    withSpinner(dataTableOutput("m_clustering_table2")), width = 7),
    ),
  box(id = "m_clustering_box10", width = 11,
 column(
  h3("Number of cells in clusters based on sample(s)", downloadBttn(outputId = "download_m_clustering_table3", label = "Download as csv")),
  withSpinner(dataTableOutput("m_clustering_table3")), width = 11),
  ),    
 ),

 box(id = "m_clustering_box8", width = 12,
  hr(),
  downloadBttn(outputId = "m_clustering", label = "Download Seurat Object"),
  actionLink("link_m_doublet", actionBttn("MDOUBLETButtonLoad", "NEXT STEP (Remove Doublets)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
  #actionLink("link_m_marker", actionBttn("MMARKERButtonLoad", "Skip Doublet, NEXT STEP (Markers Identification)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
 ), 
), 


######################################################Tab1.5########################################################   
####################################################Doublet finder#################################################
tabPanel(
  "Remove Doublets",
  fluidRow(   
    box(id = "m_doublet_box1",                      
        h3("Parameters to detect doublets"), 
        column(5, numericInput("m_doublet1", label = "Estimated percentage of doublets in dataset", min = 0.001, max = 1, value = 0.075)),
        column(4, selectInput("m_doublet3", label = "Dim plots label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        br(),
        column(width =3, actionBttn("multiple_sample_doublet", "Detect doublet",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ),	
  hr(),
  fluidRow(
    box(id = "m_doublet_box2", 
        h3("Singlet / doublets plot"), 
        actionBttn("download_m_doublet1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet1_plot", height = 600)),
    ),
    box(id = "m_doublet_box3", 
        h3("Singlet / doublets plot based on condition(s)"), 
        actionBttn("download_m_doublet2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet2_plot", height = 600)),
    ),
    box(id = "m_doublet_box4", 
        h3("Singlet / doublets plot based on sample(s)"), 
        actionBttn("download_m_doublet3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet3_plot", height = 600)),
    ),
    box(id = "m_doublet_box5", 
        h3("Singlet / doublets plot based on clusters"), 
        actionBttn("download_m_doublet4_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet4_plot", height = 600)),
    ),
    box(id = "m_doublet_box6",
    column(
      h4("Number of singlet / doublets in sample(s)"), 
      downloadBttn(outputId = "download_m_doublet1_table", label = "Download as csv"),
      withSpinner(dataTableOutput("m_doublet1_table")), 
      #withSpinner(verbatimTextOutput("m_doublet1_table")),
      width = 6),
  ),
  ),
  hr(),
  fluidRow(   
    box(id = "m_doublet_box7", width=12,                       
        h3("Parameters to keep or remove doublets"), 
        column(4, radioButtons("m_doublet2", label = "Keep or remove doublets", choices = c("Remove doublets" = "remove_doublet", "Keep doublets" = "keep_doublet"), selected = "remove_doublet", inline=T)),
        br(),
        column(width =4, actionBttn("multiple_sample_doublet2", "Keep or remove doublets",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ),
  
  fluidRow(
    box(id = "m_doublet_box8",status = "warning", 
        solidHeader = F,
        hr(),
        column(
          h3("Number of cell counts used for further analysis"), 
          downloadBttn(outputId = "download_m_doublet2_table", label = "Download as csv"),
          withSpinner(dataTableOutput("m_doublet2_table")), width = 8)
    ),
  ),
  
  fluidRow(
    box(id = "m_doublet_box9", width=10,
        h3("Singlet/doublet after keeping or removal"),
        actionBttn("download_m_doublet5_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet5_plot", height = 600)),
    ),
  ),
  fluidRow(
    box(id = "m_doublet_box10", width=10,
        h3("Based on condition(s)"),
        actionBttn("download_m_doublet6_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet6_plot", height = 900)),
    ),
  ),
  fluidRow(
    box(id = "m_doublet_box11", width=10,
        h3("Based on sample(s)"),
        actionBttn("download_m_doublet7_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet7_plot", height = 900)),
    ),
  ),
  fluidRow(
    box(id = "m_doublet_box12", width=10,
        h3("Based on clusters"),
        actionBttn("download_m_doublet8_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet8_plot", height = 900)),
    ),
  ),
  fluidRow(
    box(id = "m_doublet_box13", width=10,         
        h3("Clusters split by condition(s) and sample(s)"),
        actionBttn("download_m_doublet10_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet10_plot", height = 900)),
    ),
  ),
  fluidRow(
    box(id = "m_doublet_box16", width=10,         
        h3("Bar plots after keeping or removing doublet"),
        actionBttn("download_m_doublet9_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_doublet9_plot", height = 600)),
    ),
  ),
  fluidRow(
    box(id = "m_doublet_box14", width=5, 
    column(
      h3("Number of cells in clusters", downloadBttn(outputId = "download_m_doublet3_table", label = "Download as csv")),
      withSpinner(dataTableOutput("m_doublet3_table")), width = 5),
    ),
    box(id = "m_doublet_box17", width=7, 
    column(
      h3("Number of cells in clusters based on condition(s)", downloadBttn(outputId = "download_m_doublet4_table", label = "Download as csv")),
      withSpinner(dataTableOutput("m_doublet4_table")), width = 7),
    ),
    box(id = "m_doublet_box18", width=11, 
    column(
      h3("Number of cells in clusters based on sample(s)", downloadBttn(outputId = "download_m_doublet5_table", label = "Download as csv")),
      withSpinner(dataTableOutput("m_doublet5_table")), width = 11),
    ),
  ),
  box(id = "m_doublet_box15", width = 12,
  hr(),
  downloadBttn(outputId = "m_doublet", label = "Download Seurat Object"),
  actionLink("link_m_marker", actionBttn("MMARKERButtonLoad", "NEXT STEP (Markers Identification)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
),
),

######################################################Tab1.6########################################################   
#############################Identification of markers / Differential expression analysis#########################
tabPanel(
  "Markers Identification",
  fluidRow(   
    box(id = "m_marker_box1",
        h3("Markers identification or Differential expression analysis"),
        column(6, selectInput("m_marker1", label = "Select the analysis type", choices = c("Identify markers in all clusters" = 1, "Identify markers in one specific cluster" = 2, "Identify markers distinguishing a cluster from other cluster(s)" = 3, "Find conserved markers in one vs. all clusters" = 4, "Find conserved markers between two clusters" = 5), selected = 1)),
    ),
  ),
  fluidRow(   
    box(id = "m_marker_box2",                      
        h3("Gene expression markers parameters"), 
        column(3, numericInput("m_marker2", label = "Minimal percentage of cells", min = 0.01, max = 0.99, value = 0.25)),
        column(3, numericInput("m_marker3", label = "log fold change threshold", min = 0.01, max = 0.99, value = 0.25)),
        column(3, selectInput("m_marker4", label = "Statistical test", choices = c("wilcox"="wilcox", "wilcox_limma"="wilcox_limma", "bimod"="bimod", "roc"="roc", "t-test"="t", "LR"="LR", "MAST"="MAST"), selected = "wilcox")),
        column(3, selectInput("m_marker5", label = "Return only positive markers", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        
    ),
  ),	
  fluidRow(   
    box(id = "m_marker_box3",  
        column(6, uiOutput("m_marker_6")),
        column(6, uiOutput("m_marker_7")),
        column(6, uiOutput("m_marker_8")),
        column(6, uiOutput("m_marker_9")),
        column(3, selectInput("m_marker10", label = "group.var", choices = c("Condition" = "condition", "Samples" = "orig.ident"), selected = "condition")),
    ),
  ),
  fluidRow(   
    box(id = "m_marker_box4",                      
      column(width =3, actionBttn("multiple_sample_marker", "Detect marker genes",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  fluidRow(
    box(id = "m_marker_box5",status = "warning", width = 11,
        solidHeader = F,
        column(
          h3("Identified markers / differentially expressed genes", id = "m_marker11"), 
          h3("Conserved Markers genes", id = "m_marker12"),  
          downloadBttn(outputId = "download_m_marker1_table", label = "Download as csv"),
          withSpinner(dataTableOutput("m_marker1_table")), width = 11)
    ),
  ),
  fluidRow(
    box(id = "m_marker_box6", width=10,height = "1000px",
        br(),
        br(),
        h3("Heatmap for top 5 marker genes in cluster(s)"),
        actionBttn("download_m_marker1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_marker1_plot",height = "900px")),
        br(),
        br(),
    ),
  ),
  
  box(id = "m_marker_box7", width = 12,
      hr(),
      br(),
      downloadBttn(outputId = "m_marker", label = "Download Seurat Object"),
      actionLink("link_m_prediction", actionBttn("MANNOTATIONButtonLoad", "NEXT STEP (Cell Type Annotation)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
  ),
),
######################################################Tab1.7########################################################   
#############################################Cell type Prediction#################################################
tabPanel(
  "Cell Type Prediction",
  fluidRow(   
    box(id = "m_celltype_box1",
        h3("Predict Cell Type "),
        h5("Please make sure 'Identify markers in all clusters' were runned in the previous step, if you are using GPTCelltype"),
        column(6, selectInput("m_celltype1", label = "Cell type prediction method", choices = c("ScType" = 1, "SingleR" = 2, "GPTCelltype" = 3, "Use Own Labels" = 4), selected = 1)),
    ),
  ),
  
  fluidRow(   
    box(id = "m_celltype_box2",                      
        column(6, selectInput("m_celltype2", label = "Select reference data", choices = c("Adrenal"="Adrenal","Brain"="Brain","Eye"="Eye","Heart"="Heart","Immune system"="Immune system","Intestine"="Intestine","Kidney"="Kidney","Liver"="Liver","Lung"="Lung","Muscle"="Muscle","Pancreas"="Pancreas","Placenta"="Placenta","Spleen"="Spleen","Stomach"="Stomach","Thymus"="Thymus"))),  
    ),
  ),
  fluidRow(   
    box(id = "m_celltype_box3",                      
        column(6, selectInput("m_celltype3", label = "Select tissue", choices = c("Human primary cell atlas"="hpca","Blueprint/ENCODE"="blueprint_encode","Mouse RNA-seq"="mouse_rnaseq","Immunological Genome Project"="immgen","Database of Immune Cell Expression/eQTLs/Epigenomics"="dice","Novershtern hematopoietic data"="novershtern_hematopoietic","Monaco immune data"="monaco_immune"))),  
        column(6, selectInput("m_celltype4", label = "DE.method", choices = c("classic"="classic", "wilcox"="wilcox", "t"="t"))),  
    ),
  ),  
  fluidRow(   
    box(id = "m_celltype_box4",                      
        column(6, selectInput("m_celltype5", label = "Select model", choices = c("gpt-4"="gpt-4","gpt-4-turbo"="gpt-4-turbo","gpt-4o-mini"="gpt-4o-mini","gpt-4o"="gpt-4o","chatgpt-4o-latest"="chatgpt-4o-latest","gpt-3.5-turbo"="gpt-3.5-turbo"), selected = "gpt-4")),  
        column(6, numericInput("m_celltype6", label = "Top gene numbers to predict cell type",  min = 1, max = 25, value = 10)),  
    ),
  ),
  fluidRow(   
    box(id = "m_celltype_box5",
        uiOutput("m_celltype7"),
    ),
  ),
  fluidRow(   
    box(id = "m_celltype_box6", 
        column(4, selectInput("m_celltype8", label = "Dim plots label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        br(),
        column(3, actionBttn("multiple_sample_celltype", "Detect cell type",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  
  fluidRow(
    box(id = "m_celltype_box7", width=10,height = "1000px",
        br(),
        br(),
        h3("Dimplot of annotated clusters"),
        actionBttn("download_m_celltype1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_celltype1_plot",height = "900px")),
        #verbatimTextOutput("celltype_text_level"),
        br(),
        br(),
    ),
  ),
  fluidRow(
    box(id = "m_celltype_box8",status = "warning", width = 11,
        solidHeader = F,
        hr(),
        column(
          h3("ScType scores", id = "m_celltype9"), 
          h3("SingleR Scores", id = "m_celltype10"), 
          downloadBttn(outputId = "download_m_celltype1_table", label = "Download as csv"),
          withSpinner(dataTableOutput("m_celltype1_table")), width = 11)
    ),
  ),
  fluidRow(
    box(id = "m_celltype_box9",
        br(),
        h3("SingleR score heatmap"),
        actionBttn("download_m_celltype2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("m_celltype2_plot")),
        br(),
        h3("SingleR Delta distribution"),
        actionBttn("download_m_celltype3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
        withSpinner(plotOutput("m_celltype3_plot")),
        br(),
  ),
  ),
  box(id = "m_celltype_box10", width = 12,
      hr(),
      br(),
      downloadBttn(outputId = "m_celltype", label = "Download Seurat Object"),
      actionLink("link_m_clusterbased", actionBttn("MCPlotsButtonLoad", "NEXT STEP (Cluster Based Plots)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
  ),
),
######################################################Tab1.8########################################################   
#############################################Cluster-based plots#################################################
tabPanel(
  "Cluster-based plots",
  fluidRow(   
    box(id = "m_clusterbased_box1",                      
        h3("Select the plot type to display"), 
        h5("Please make sure 'Identify markers in all clusters' and the same 'cell prediction method' were runned in the previous steps."),
        column(6, selectInput("m_clusterbased1", label = "No. of features to display", choices = c("1" = 1, "2" = 2, "3"= 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10, "List of gene names" = "gene_name_list"), selected = 3)),
        column(6, textAreaInput("m_clusterbased2", label ="Enter your genes for ploting (eg: gene names separated by , )", value = "", placeholder = "KLK2,KLK3,CTSG,MS4A3,CLEC4OP,KDR", height = '400px')),
        column(3, selectInput("m_clusterbased3", label = "Plot type", choices = c("Dot Plot" = "Dot Plot", "Violin Plot" = "VlnPlot", "Ridge Plot"= "RidgePlot", "Feature Plot" = "FeaturePlot"), selected = "Dot Plot")),
        column(3, selectInput("m_clusterbased4", label = "group.by", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        column(3, selectInput("m_clusterbased5", label = "split.by", choices = c("None" = "NULL", "Condition" = "condition", "Samples" = "orig.ident"), selected = "NULL")),
        br(),
        br(),
        column(3, actionBttn("multiple_sample_clusterbased", "Generate plots",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  fluidRow(
    box(id = "m_clusterbased_box2", width=10,height = "2500px",

        br(),
        h3("Dot / Violin / Ridge / Feature plot"),
        actionBttn("download_m_clusterbased1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
        status = "warning",
        solidHeader = F,
        withSpinner(plotOutput("m_clusterbased1_plot",height = "2400px")),
    ),
  ),
  fluidRow(
    box(id = "m_clusterbased_box3",status = "warning", width = 11,
        solidHeader = F,
        hr(),
        column(
          h3("Top or selected genes, cell counts and proportion"), 
          downloadBttn(outputId = "download_m_clusterbased1_table", label = "Download as csv"),
          withSpinner(dataTableOutput("m_clusterbased1_table")), width = 11)
    ),
  ),
  box(id = "m_clusterbased_box4", width = 12,
      hr(),
      br(),
      downloadBttn(outputId = "m_clusterbased", label = "Download Seurat Object"),
      actionLink("link_m_conditionbased", actionBttn("MCLUSTERButtonLoad", "NEXT STEP (Condition based analysis)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
  ),
),
######################################################Tab1.9########################################################   
#############################################Condition-based analysis#################################################
tabPanel(
  "Condition based analysis",
  fluidRow(   
    box(id = "m_conditionbased_box1",     
        h3(" Differential expression analysis between two groups"),
        column(6, uiOutput("m_conditionbased_1")),
        column(6, uiOutput("m_conditionbased_2")),
        h3("Parameters to find the DEGs"), 
        column(3, numericInput("m_conditionbased3", label = "Minimal percentage of cells", min = 0.01, max = 0.99, value = 0.25)),
        column(3, numericInput("m_conditionbased4", label = "log fold change threshold", min = 0.01, max = 0.99, value = 0.25)),
        column(3, selectInput("m_conditionbased5", label = "Statistical test", choices = c("wilcox"="wilcox", "wilcox_limma"="wilcox_limma", "bimod"="bimod", "roc"="roc", "t-test"="t", "LR"="LR", "MAST"="MAST"), selected = "wilcox")),
        column(3, selectInput("m_conditionbased6", label = "Return only positive markers", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
    ),
  ),
  fluidRow( 
    box(id = "m_conditionbased_box2",  
        h3("Parameters for ploting"),
        column(3, selectInput("m_conditionbased7", label = "Plot type", choices = c("Dot Plot" = "Dot Plot", "Violin Plot" = "VlnPlot", "Ridge Plot"= "RidgePlot", "Feature Plot" = "FeaturePlot", "Volcano Plot" = "VolcanoPlot"), selected = "Dot Plot")),
        column(3, selectInput("m_conditionbased8", label = "group.by", choices = c("Condition" = "condition", "Samples" = "orig.ident"), selected = "condition")),
        column(6, numericInput("m_conditionbased9", label = "No. of features to display",  min = 1, max = 50, value = 15)),
        br(),
        column(3, actionBttn("multiple_sample_conditionbased", "Submit",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  fluidRow(
    box(id = "m_conditionbased_box3", width=10,height = "2000px",
        h3("Dot / Violin / Ridge / Feature / Volcano plot"),
        actionBttn("download_m_conditionbased1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
        status = "warning",
        solidHeader = F,
        withSpinner(plotOutput("m_conditionbased1_plot",height = "1900px")),
    ),
  ),
  fluidRow(
    box(id = "m_conditionbased_box4",status = "warning", width = 11,
        solidHeader = F,
        column(
          h3("Differentially expressed genes"), 
          downloadBttn(outputId = "download_m_conditionbased1_table", label = "Download as csv"),
          withSpinner(dataTableOutput("m_conditionbased1_table")), width = 11)
    ),
  ),
  box(id = "m_conditionbased_box5", width = 12,
      hr(),
      br(),
      downloadBttn(outputId = "m_conditionbased", label = "Download Seurat Object"),
      #actionLink("link_m_conditionbased", actionBttn("MCONDITIONButtonLoad", "NEXT STEP (Cell Type Prediction)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
  ),
),


       
     ),
   ),
 ),
),
########################################################################Tab2########################################################################################## 
#################################################################Multiple samples Sub clustering#####################################################################
################################################################################################################################################################## 
tabPanel(
  "Subclustering",
  useShinyjs(),
  sidebarLayout(
    sidebarPanel(id="subclustering_multiple_sidebar",
                 h3(id = "m_subclustering0", "Please run upto Single or Multiple samples Analysis upto Celltype prediction to begin this analysis"),
                 selectInput("m_subclustering1", label = "Select the cluster type for sub clustering", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters"),
                 uiOutput("m_subclustering_2"),
                 uiOutput("m_subclustering_3"),
                 # uiOutput("m_subclustering_4"),
                 # uiOutput("m_subclustering_5"),
                 # uiOutput("m_subclustering_6"),
                 actionBttn("subclustering_multiple_sample_submit", "Submit", style = "unite",color = "primary"),
                 width = 3),
    mainPanel(id="subclustering_multiple_main_menu",
              tabsetPanel(
                type = "tabs", id="subclustering_multiple_tabsets",
                
######################################################Tab2.1########################################################                 
######################################################STATS#######################################################
                tabPanel(
                  "Cell Stats", 
                  box(id = "m_subclustering_box1", width = 10,
                      # h4("Total number of cells in the selected clusters"),
                      # withSpinner(verbatimTextOutput("subclustering_multiple_cell_text_level")),
                      h3("Number of cells in the sample(s)"), 
                      downloadBttn(outputId = "download_subclustering_multiple_cell_table", label = "Download as csv"),
                      withSpinner(dataTableOutput("subclustering_multiple_cell_table"))
                  ),
                  fluidRow(
                    box(id = "m_subclustering_box2", 
                        h3("QC Stats for sleected sub clusters"), 
                        actionBttn("download_m_subclustering_qc", "Download plot",  style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_qc", width = "1000px", height = "500px"))),
                  ),
                  box(id = "m_subclustering_box3", width =12,
                      hr(),
                      downloadBttn(outputId = "m_subclustering_stats", label = "Download Seurat Object"),
                      actionLink("link_m_subclustering_normalization", actionBttn("MSUBQCButtonLoad", "NEXT STEP (Normalization and PCA analysis)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ),
                ),
                
#####################################################Tab2.2####################################  
#############################multiple_Normalization and PCA analysis#########################
                tabPanel(
                  "Normalization and PCA analysis",
                  h5("Please use the same normalization method used in single or multiple samples analysis"),
                  selectInput("subclustering_multiple_sample_normalization_method", label = "Normalization method", choices = c("LogNormalize"="LogNormalize", "SCTransform"="SCTransform"), selected = "LogNormalize"),
                  fluidRow( 
                    #box(id = "boxLogNormalize",
                    column(width =4, numericInput("subclustering_multiple_sample_scale_factor", label = "Scale factor", value = 10000)),
                    column(width =4, selectInput("subclustering_multiple_sample_normalization_variable_genes", label = "variable genes detection", choices = c( "vst"="vst", "mean.var.plot (mvp)"="mean.var.plot", "dispersion (disp)"="dispersion"), selected = "vst")),  
                    column(width =4, numericInput("subclustering_multiple_sample_var_genes", label = "Number of top variable features", value = 2000)),
                    #),
                    column(width =4, numericInput("subclustering_multiple_sample_pca_dim", label = "Number of dimensions (PCA)", value = 30)),
                    br(),
                    column(width =4, actionBttn("subclustering_multiple_sample_normalization", "Submit",  style = "unite",color = "primary", icon = icon("download"))),
                  ),
                  hr(),
                  
                  fluidRow(
                    box(id = "m_subclustering_pca_box1", 
                        h3("Dimension reduction heatmap  for PCA data"), 
                        actionBttn("download_m_subclustering_pca_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_pca_plot", height = 600))),
                    
                    box(id = "m_subclustering_elbow_box", 
                        h3("Elbow plot"), 
                        actionBttn("download_m_subclustering_elbow_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_elbow_plot", height = 600))),
                  ),
                  
                  fluidRow(
                    
                    box(id = "m_subclustering_pca_box2", 
                        h3("PCA sample(s) based"), 
                        actionBttn("download_m_subclustering_pca2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_pca2_plot", height = 600))),
                    
                    box(id = "m_subclustering_pca_box3", 
                        h3("PCA group(s) based"), 
                        actionBttn("download_m_subclustering_pca3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_pca3_plot", height = 600))),
                  ),
                  box(id = "m_subclustering_pca_box4", width = 12,
                      hr(),
                      downloadBttn(outputId = "m_subclustering_normalization", label = "Download Seurat Object"),
                      actionLink("link_m_subclustering_clustering", actionBttn("MSUBCLUSTERINGButtonLoad", "NEXT STEP (Clustering)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ), 
                ), 
                
###########################################Tab2.3########################################   
#############################subclustering_multiple_Clustering#########################
                tabPanel(
                  "Clustering",
                  fluidRow(
                    h3("Nearest-neighbour graph construction"),
                    column(3, numericInput("m_subclustering_clustering1", label = "Number of dimensions", min = 0, max = 100, value = 30, step = 1)),
                    column(3, numericInput("m_subclustering_clustering2", label = "k.param", value = 20, step = 1)),
                    column(3, numericInput("m_subclustering_clustering3", label = "n.trees", value = 50, step = 1)),
                  ),    
                  fluidRow(		
                    h3("Clustering parameters and integration method"),
                    column(3, numericInput("m_subclustering_clustering4", label = "Resolution", min = 0, max = 10, value = 0.5, step = 0.1)),
                    column(3, selectInput("m_subclustering_clustering5", label = "Clustering algorithm", choices = c("Louvain" = 1, "SLM" = 3, "Leiden" = 4), selected = 1)),
                    # ),
                    # fluidRow(
                    #   h3("Select integration method"),
                    column(3, selectInput("m_subclustering_clustering13", label = "Integration method (except none, minimum two samples are required)", choices = c("None" = "None", "HarmonyIntegration" = "HarmonyIntegration", "CCAIntegration" = "CCAIntegration", "RPCAIntegration" = "RPCAIntegration", JointPCAIntegration = "JointPCAIntegration" ), selected = "None")),
                  ),                 
                  fluidRow( 
                    h3("Dimension reduction"),
                    column(3, radioButtons("m_subclustering_clustering6", label = "Plot Options", choices = c("UMAP" = "umap", "t-SNE" = "tsne"), selected = "umap", inline=T)),
                  ),                 
                  fluidRow(   
                    box(id = "m_subclustering_umap_box",                      
                        h3("UMAP parameters"), 
                        column(3, numericInput("m_subclustering_clustering7", label = "Number of dimensions", min = 0, max = 100, value = 30, step = 1)),
                        column(3, numericInput("m_subclustering_clustering8", label = "k-nearest-neighbours", min = 0, max = 50, value = 20, step = 1)),
                        column(3, numericInput("m_subclustering_clustering9", label = "min.dist", min = 0.001, max = 0.5, value = 0.3, step = 0.01)),
                        column(3, selectInput("m_subclustering_clustering10", label = "Show label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
                    ),
                  ),		
                  fluidRow(                  
                    box(id = "m_subclustering_tsne_box",                      
                        h3("t-SNE parameters"),                 
                        column(3, numericInput("m_subclustering_clustering11", label = "Number of dimensions",  min = 0, max = 100, value = 30, step = 1)),
                        column(3, selectInput("m_subclustering_clustering12", label = "Show label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
                    ),
                  ),   
                  fluidRow(
                    column(width =3, actionBttn("subclustering_multiple_sample_clustering", "Submit",  style = "unite",color = "primary", icon = icon("download"))),
                  ),
                  hr(),
                  
                  fluidRow(
                    box(id = "m_subclustering_clustering_box1", 
                        h3("UMAP / t-SNE cluster plot"), 
                        actionBttn("download_m_subclustering_umap_tsne1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne1_plot", height = 600))),
                    
                    box(id = "m_subclustering_clustering_box2", 
                        h3("Cluster based count bar plot"), 
                        actionBttn("download_m_subclustering_umap_tsne_bar1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne_bar1_plot", height = 600))),
                  ),
                  
                  fluidRow(
                    box(id = "m_subclustering_clustering_box3", 
                        h3("UMAP / t-SNE condition(s) based plot"), 
                        actionBttn("download_m_subclustering_umap_tsne2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne2_plot", height = 600))),
                    
                    box(id = "m_subclustering_clustering_box4", 
                        h3("Condition(s) based count bar plot"), 
                        actionBttn("download_m_subclustering_umap_tsne_bar2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne_bar2_plot", height = 600))),
                  ),
                  
                  fluidRow(
                    box(id = "m_subclustering_clustering_box5", 
                        h3("UMAP / t-SNE sample(s) based plot"), 
                        actionBttn("download_m_subclustering_umap_tsne3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne3_plot", height = 600))),
                    
                    box(id = "m_subclustering_clustering_box6", 
                        h3("Sample(s) based  count Bar plot"), 
                        actionBttn("download_m_subclustering_umap_tsne_bar3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne_bar3_plot", height = 600))),
                  ),
                  fluidRow(
                    box(id = "m_subclustering_clustering_box11", width=10,         
                        h3("Clusters split by condition(s) and sample(s)"),
                        actionBttn("download_m_subclustering_umap_tsne4_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_umap_tsne4_plot", height = 900)),
                    ),
                  ),
                  
                  fluidRow(
                    box(id = "m_subclustering_clustering_box7", width = 5,
                        column(
                          h3("Number of cells in clusters", downloadBttn(outputId = "download_m_subclustering_clustering_table1", label = "Download as csv")),
                          withSpinner(dataTableOutput("m_subclustering_clustering_table1")), width = 5),
                    ),
                    box(id = "m_subclustering_clustering_box9", width = 7,
                        column(
                          h3("Number of cells in clusters based on condition(s)", downloadBttn(outputId = "download_m_subclustering_clustering_table2", label = "Download as csv")),
                          withSpinner(dataTableOutput("m_subclustering_clustering_table2")), width = 7),
                    ),
                    box(id = "m_subclustering_clustering_box10", width = 11,
                        column(
                          h3("Number of cells in clusters based on sample(s)", downloadBttn(outputId = "download_m_subclustering_clustering_table3", label = "Download as csv")),
                          withSpinner(dataTableOutput("m_subclustering_clustering_table3")), width = 11),
                    ),    
                  ),
                  
                  box(id = "m_subclustering_clustering_box8", width = 12,
                      hr(),
                      downloadBttn(outputId = "m_subclustering_clustering", label = "Download Seurat Object"),
                      actionLink("link_m_subclustering_marker", actionBttn("MSUBMARKERButtonLoad", "NEXT STEP (Markers Identification)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ), 
                ), 
                
                
######################################################Tab2.4########################################################   
#############################Identification of markers / Differential expression analysis#########################
                tabPanel(
                  "Markers Identification",
                  fluidRow(   
                    box(id = "m_subclustering_marker_box1",
                        h3("Markers identification or Differential expression analysis"),
                        column(6, selectInput("m_subclustering_marker1", label = "Select the analysis type", choices = c("Identify markers in all clusters" = 1, "Identify markers in one specific cluster" = 2, "Identify markers distinguishing a cluster from other cluster(s)" = 3, "Find conserved markers in one vs. all clusters" = 4, "Find conserved markers between two clusters" = 5), selected = 1)),
                    ),
                  ),
                  fluidRow(   
                    box(id = "m_subclustering_marker_box2",                      
                        h3("Gene expression markers parameters"), 
                        column(3, numericInput("m_subclustering_marker2", label = "Minimal percentage of cells", min = 0.01, max = 0.99, value = 0.25)),
                        column(3, numericInput("m_subclustering_marker3", label = "log fold change threshold", min = 0.01, max = 0.99, value = 0.25)),
                        column(3, selectInput("m_subclustering_marker4", label = "Statistical test", choices = c("wilcox"="wilcox", "wilcox_limma"="wilcox_limma", "bimod"="bimod", "roc"="roc", "t-test"="t", "LR"="LR", "MAST"="MAST"), selected = "wilcox")),
                        column(3, selectInput("m_subclustering_marker5", label = "Return only positive markers", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
                        
                    ),
                  ),	
                  fluidRow(   
                    box(id = "m_subclustering_marker_box3",  
                        column(6, uiOutput("m_subclustering_marker_6")),
                        column(6, uiOutput("m_subclustering_marker_7")),
                        column(6, uiOutput("m_subclustering_marker_8")),
                        column(6, uiOutput("m_subclustering_marker_9")),
                        column(3, selectInput("m_subclustering_marker10", label = "group.var", choices = c("Condition" = "condition", "Samples" = "orig.ident"), selected = "condition")),
                    ),
                  ),
                  fluidRow(   
                    box(id = "m_subclustering_marker_box4",                      
                        column(width =3, actionBttn("subclustering_multiple_sample_marker", "Detect marker genes",  style = "unite",color = "primary", icon = icon("download"))),
                    ),
                  ), 
                  hr(),
                  fluidRow(
                    box(id = "m_subclustering_marker_box5",status = "warning", width = 11,
                        solidHeader = F,
                        column(
                          h3("Identified markers / differentially expressed genes", id = "m_subclustering_marker11"), 
                          h3("Conserved Markers genes", id = "m_subclustering_marker12"),  
                          downloadBttn(outputId = "download_m_subclustering_marker1_table", label = "Download as csv"),
                          withSpinner(dataTableOutput("m_subclustering_marker1_table")), width = 11)
                    ),
                  ),
                  fluidRow(
                    box(id = "m_subclustering_marker_box6", width=10,height = "1000px",
                        br(),
                        br(),
                        h3("Heatmap for top 5 marker genes in cluster(s)"),
                        actionBttn("download_m_subclustering_marker1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_marker1_plot",height = "900px")),
                        br(),
                        br(),
                    ),
                  ),
                  
                  box(id = "m_subclustering_marker_box7", width = 12,
                      hr(),
                      br(),
                      downloadBttn(outputId = "m_subclustering_marker", label = "Download Seurat Object"),
                      actionLink("link_m_subclustering_prediction", actionBttn("MSUBANNOTATIONButtonLoad", "NEXT STEP (Cell Type Annotation)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ),
                ),
######################################################Tab2.5########################################################   
#############################################Cell type Prediction#################################################
                tabPanel(
                  "Cell Type Prediction",
                  fluidRow(   
                    box(id = "m_subclustering_celltype_box1",
                        h3("Predict Cell Type "),
                        h5("Please make sure 'Identify markers in all clusters' were runned in the previous step, if you are using GPTCelltype"),
                        column(6, selectInput("m_subclustering_celltype1", label = "Cell type prediction method", choices = c("ScType" = 1, "SingleR" = 2, "GPTCelltype" = 3, "Use Own Labels" = 4), selected = 1)),
                    ),
                  ),
                  
                  fluidRow(   
                    box(id = "m_subclustering_celltype_box2",                      
                        column(6, selectInput("m_subclustering_celltype2", label = "Select reference data", choices = c("Adrenal"="Adrenal","Brain"="Brain","Eye"="Eye","Heart"="Heart","Immune system"="Immune system","Intestine"="Intestine","Kidney"="Kidney","Liver"="Liver","Lung"="Lung","Muscle"="Muscle","Pancreas"="Pancreas","Placenta"="Placenta","Spleen"="Spleen","Stomach"="Stomach","Thymus"="Thymus"))),  
                    ),
                  ),
                  fluidRow(   
                    box(id = "m_subclustering_celltype_box3",                      
                        column(6, selectInput("m_subclustering_celltype3", label = "Select tissue", choices = c("Human primary cell atlas"="hpca","Blueprint/ENCODE"="blueprint_encode","Mouse RNA-seq"="mouse_rnaseq","Immunological Genome Project"="immgen","Database of Immune Cell Expression/eQTLs/Epigenomics"="dice","Novershtern hematopoietic data"="novershtern_hematopoietic","Monaco immune data"="monaco_immune"))),  
                        column(6, selectInput("m_subclustering_celltype4", label = "DE.method", choices = c("classic"="classic", "wilcox"="wilcox", "t"="t"))),  
                    ),
                  ),  
                  fluidRow(   
                    box(id = "m_subclustering_celltype_box4",                      
                        column(6, selectInput("m_subclustering_celltype5", label = "Select model", choices = c("gpt-4"="gpt-4","gpt-4-turbo"="gpt-4-turbo","gpt-4o-mini"="gpt-4o-mini","gpt-4o"="gpt-4o","chatgpt-4o-latest"="chatgpt-4o-latest","gpt-3.5-turbo"="gpt-3.5-turbo"), selected = "gpt-4")),  
                        column(6, numericInput("m_subclustering_celltype6", label = "Top gene numbers to predict cell type",  min = 1, max = 25, value = 10)),  
                    ),
                  ),
                  fluidRow(   
                    box(id = "m_subclustering_celltype_box5",
                        uiOutput("m_subclustering_celltype7"),
                    ),
                  ),
                  fluidRow(   
                    box(id = "m_subclustering_celltype_box6", 
                        column(4, selectInput("m_subclustering_celltype8", label = "Dim plots label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
                        br(),
                        column(3, actionBttn("subclustering_multiple_sample_celltype", "Detect cell type",  style = "unite",color = "primary", icon = icon("download"))),
                    ),
                  ), 
                  hr(),
                  
                  fluidRow(
                    box(id = "m_subclustering_celltype_box7", width=10,height = "1000px",
                        br(),
                        br(),
                        h3("Dimplot of annotated clusters"),
                        actionBttn("download_m_subclustering_celltype1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_celltype1_plot",height = "900px")),
                        #verbatimTextOutput("celltype_text_level"),
                        br(),
                        br(),
                    ),
                  ),
                  fluidRow(
                    box(id = "m_subclustering_celltype_box8",status = "warning", width = 11,
                        solidHeader = F,
                        hr(),
                        column(
                          h3("ScType scores", id = "m_subclustering_celltype9"), 
                          h3("SingleR Scores", id = "m_subclustering_celltype10"), 
                          downloadBttn(outputId = "download_m_subclustering_celltype1_table", label = "Download as csv"),
                          withSpinner(dataTableOutput("m_subclustering_celltype1_table")), width = 11)
                    ),
                  ),
                  fluidRow(
                    box(id = "m_subclustering_celltype_box9",
                        br(),
                        h3("SingleR score heatmap"),
                        actionBttn("download_m_subclustering_celltype2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        status = "warning", 
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_celltype2_plot")),
                        br(),
                        h3("SingleR Delta distribution"),
                        actionBttn("download_m_subclustering_celltype3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
                        withSpinner(plotOutput("m_subclustering_celltype3_plot")),
                        br(),
                    ),
                  ),
                  box(id = "m_subclustering_celltype_box10", width = 12,
                      hr(),
                      br(),
                      downloadBttn(outputId = "m_subclustering_celltype", label = "Download Seurat Object"),
                      actionLink("link_m_subclustering_clusterbased", actionBttn("MSUBCPlotsButtonLoad", "NEXT STEP (Cluster Based Plots)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ),
                ),
#####################################################Tab2.6########################################################   
#############################################Cluster-based plots#################################################
                tabPanel(
                  "Cluster-based plots",
                  fluidRow(   
                    box(id = "m_subclustering_clusterbased_box1",                      
                        h3("Select the plot type to display"), 
                        h5("Please make sure 'Identify markers in all clusters' and the same 'cell prediction method' were runned in the previous steps."),
                        column(6, selectInput("m_subclustering_clusterbased1", label = "No. of features to display", choices = c("1" = 1, "2" = 2, "3"= 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10, "List of gene names" = "gene_name_list"), selected = 3)),
                        column(6, textAreaInput("m_subclustering_clusterbased2", label ="Enter your genes for ploting (eg: gene names separated by , )", value = "", placeholder = "KLK2,KLK3,CTSG,MS4A3,CLEC4OP,KDR", height = '400px')),
                        column(3, selectInput("m_subclustering_clusterbased3", label = "Plot type", choices = c("Dot Plot" = "Dot Plot", "Violin Plot" = "VlnPlot", "Ridge Plot"= "RidgePlot", "Feature Plot" = "FeaturePlot"), selected = "Dot Plot")),
                        column(3, selectInput("m_subclustering_clusterbased4", label = "group.by", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
                        column(3, selectInput("m_subclustering_clusterbased5", label = "split.by", choices = c("None" = "NULL", "Condition" = "condition", "Samples" = "orig.ident"), selected = "NULL")),
                        br(),
                        br(),
                        column(3, actionBttn("subclustering_multiple_sample_clusterbased", "Generate plots",  style = "unite",color = "primary", icon = icon("download"))),
                    ),
                  ), 
                  hr(),
                  fluidRow(
                    box(id = "m_subclustering_clusterbased_box2", width=10,height = "2500px",

                        br(),
                        h3("Dot / Violin / Ridge / Feature plot"),
                        actionBttn("download_m_subclustering_clusterbased1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
                        status = "warning",
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_clusterbased1_plot",height = "2400px")),
                    ),
                  ),
                  fluidRow(
                    box(id = "m_subclustering_clusterbased_box3",status = "warning", width = 11,
                        solidHeader = F,
                        hr(),
                        column(
                          h3("Top or selected genes, cell counts and proportion"), 
                          downloadBttn(outputId = "download_m_subclustering_clusterbased1_table", label = "Download as csv"),
                          withSpinner(dataTableOutput("m_subclustering_clusterbased1_table")), width = 11)
                    ),
                  ),
                  box(id = "m_subclustering_clusterbased_box4", width = 12,
                      hr(),
                      br(),
                      downloadBttn(outputId = "m_subclustering_clusterbased", label = "Download Seurat Object"),
                      actionLink("link_m_subclustering_conditionbased", actionBttn("MSUBCLUSTERButtonLoad", "NEXT STEP (Condition based analysis)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ),
                ),
######################################################Tab2.7############################################################   
#############################################Condition-based analysis#################################################
                tabPanel(
                  "Condition based analysis",
                  fluidRow(   
                    box(id = "m_subclustering_conditionbased_box1",     
                        h3(" Differential expression analysis between two groups"),
                        column(6, uiOutput("m_subclustering_conditionbased_1")),
                        column(6, uiOutput("m_subclustering_conditionbased_2")),
                        h3("Parameters to find the DEGs"), 
                        column(3, numericInput("m_subclustering_conditionbased3", label = "Minimal percentage of cells", min = 0.01, max = 0.99, value = 0.25)),
                        column(3, numericInput("m_subclustering_conditionbased4", label = "log fold change threshold", min = 0.01, max = 0.99, value = 0.25)),
                        column(3, selectInput("m_subclustering_conditionbased5", label = "Statistical test", choices = c("wilcox"="wilcox", "wilcox_limma"="wilcox_limma", "bimod"="bimod", "roc"="roc", "t-test"="t", "LR"="LR", "MAST"="MAST"), selected = "wilcox")),
                        column(3, selectInput("m_subclustering_conditionbased6", label = "Return only positive markers", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
                    ),
                  ),

                  fluidRow( 
                    box(id = "m_subclustering_conditionbased_box2",  
                        h3("Parameters for ploting"),
                        column(3, selectInput("m_subclustering_conditionbased7", label = "Plot type", choices = c("Dot Plot" = "Dot Plot", "Violin Plot" = "VlnPlot", "Ridge Plot"= "RidgePlot", "Feature Plot" = "FeaturePlot", "Volcano Plot" = "VolcanoPlot"), selected = "Dot Plot")),
                        column(3, selectInput("m_subclustering_conditionbased8", label = "group.by", choices = c("Condition" = "condition", "Samples" = "orig.ident"), selected = "condition")),
                        column(6, numericInput("m_subclustering_conditionbased9", label = "No. of features to display",  min = 1, max = 50, value = 15)),
                        br(),
                        column(3, actionBttn("subclustering_multiple_sample_conditionbased", "Submit",  style = "unite",color = "primary", icon = icon("download"))),
                    ),
                  ), 
                  hr(),
                  fluidRow(
                    box(id = "m_subclustering_conditionbased_box3", width=10,height = "2000px",
                        br(),
                        h3("Dot / Violin / Ridge / Feature / Volcano plot"),
                        actionBttn("download_m_subclustering_conditionbased1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
                        status = "warning",
                        solidHeader = F,
                        withSpinner(plotOutput("m_subclustering_conditionbased1_plot",height = "1900px")),
                    ),
                  ),
                  fluidRow(
                    box(id = "m_subclustering_conditionbased_box4",status = "warning", width = 11,
                        solidHeader = F,
                        column(
                          h3("Differentially expressed genes"), 
                          downloadBttn(outputId = "download_m_subclustering_conditionbased1_table", label = "Download as csv"),
                          withSpinner(dataTableOutput("m_subclustering_conditionbased1_table")), width = 11)
                    ),
                  ),
                  box(id = "m_subclustering_conditionbased_box5", width = 12,
                      hr(),
                      br(),
                      downloadBttn(outputId = "m_subclustering_conditionbased", label = "Download Seurat Object"),
                      #actionLink("link_m_subclustering_conditionbased", actionBttn("MSUBCONDITIONButtonLoad", "NEXT STEP (Cell Type Prediction)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                  ),
                ),
                
                
                
              ),
    ),
  ),
),
######################################################Menu3#####################################################################   
#########################################Cell Cluster Correlation Network####################################################################
tabPanel(
  "Correlation network",
  useShinyjs(),
  fluidRow( 
    box(id = "s_cccn_box0",  width=12,
        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction and Marker Identification step."),
    ),
  ),
  
  fluidRow(   
    box(id = "s_cccn_box1",
        h3("Cell cluster correlation network analysis"),
        h3("Select the input data and celltype method for analysis"),
        column(3, selectInput("s_cccn1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering"), selected = "multiple_sample")),
        column(3, selectInput("s_cccn2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        column(3, selectInput("s_cccn3", label = "Correlation method", choices = c("pearson"="pearson","spearman"="spearman","kendall"="kendall"), selected = "spearman")),
        column(3, actionBttn("single_multiple_sample_cccn", "Cluster correlation network",  style = "unite",color = "primary", icon = icon("download"))),		
    ),
  ),
  hr(),
  fluidRow(
    box(id = "s_cccn_box2",status = "warning", width=12,
        
        h3("Cluster-based correlation matrix plot"),
        # actionBttn("download_s_cccn1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cccn1_plot", height = "600px")),
        h3("Cluster-based Correlation Network plot"),
        actionBttn("download_s_cccn2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cccn2_plot", height = "600px")),
        br(),
        h3("Cluster-based correlation table"), 
        downloadBttn(outputId = "download_s_cccn1_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_cccn1_table")),
        hr(),
    ),
  ),
),
             
######################################################Menu4#####################################################################   
####################################################GO terms####################################################################
tabPanel(
  "GO terms",
  useShinyjs(),
  fluidRow( 
    box(id = "s_go_box0",  width=12,
        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction and Marker Identification step."),
    ),
  ),
  
  fluidRow(   
    box(id = "s_go_box1",
        h3("Select the input data and cluster(s) for analysis"),
        column(6, selectInput("s_go1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering", "List of gene names" = "gene_name_list"), selected = "multiple_sample")),
        column(6, textAreaInput("s_go14", label ="Enter your genes (eg: gene names separated by , )", value = "", placeholder = "KLK2,KLK3,CTSG,MS4A3,CLEC4OP,KDR", height = '400px')),
        column(6, selectInput("s_go2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        column(6, uiOutput("s_go_3")),
        column(6, numericInput("s_go4", label = "p_val_adj", value = 0.05)),
         
    ),
        box(id = "s_go_box2",
        h3("GO term parameters"),
        column(4, selectInput("s_go5", label = "Organism", choices = c("Human"="org.Hs.eg.db", "Mouse"="org.Mm.eg.db", "Rat"="org.Rn.eg.db", "Pig"="org.Ss.eg.db", "Rhesus"="org.Mmu.eg.db"), selected = "org.Hs.eg.db")),
        column(4, selectInput("s_go6", label = "Ontology", choices = c("BP"="BP", "MF"="MF", "CC"="CC", "ALL"="ALL"), selected = "BP")),
        column(4, selectInput("s_go7", label = "pAdjustMethod", choices = c("holm"="holm", "hochberg"="hochberg", "hommel"="hommel", "bonferroni"="bonferroni", "BH"="BH", "BY"="BY", "fdr"="fdr", "none"="none"), selected = "BH")),
        column(4, numericInput("s_go8", label = "pvalueCutoff", value = 0.05)),
        column(4, numericInput("s_go9", label = "qvalueCutoff", value = 0.2)),
        column(4, numericInput("s_go10", label = "Minimal size of genes", value = 10)),
        column(4, numericInput("s_go11", label = "Maximal size of genes", value = 500)),
        column(4, selectInput("s_go12", label = "Plot type", choices = c("Dot plot"="dotplot", "Bar plot"="barplot", "Net plot"="cnetplot", "upset plot"="upsetplot"), selected = "dotplot")),
        column(4, numericInput("s_go13", label = "No. of category to plot", value = 10, min = 1, max = 50)),
        column(4, actionBttn("single_multiple_sample_go", "Go term",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  fluidRow(
    box(id = "s_go_box3",status = "warning", width = 11,
        
        h3("Go term plot"),
        actionBttn("download_s_go1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_go1_plot", height = "600px")),
        br(),
        h3("Summary table"), 
        downloadBttn(outputId = "download_s_go1_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_go1_table")),
        hr(),
    ),
  ),  

  ),
######################################################Menu5#####################################################################   
####################################################Pathway analysis####################################################################
tabPanel(
  "Pathway analysis",
  useShinyjs(),
  fluidRow( 
    box(id = "s_pathway_box0",  width=12,
        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction and Marker Identification step."),
    ),
  ),
  
  fluidRow(   
    box(id = "s_pathway_box1",
        h3("Select the input data and cluster(s) for analysis"),
        column(4, selectInput("s_pathway6", label = "Pathway analysis type", choices = c("KEGG"="KEGG", "Reactome"="Reactome"), selected = "KEGG")),
        column(4, selectInput("s_pathway1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering", "List of gene names" = "gene_name_list"), selected = "multiple_sample")),
        column(4, textAreaInput("s_pathway14", label ="Enter your genes (eg: gene names separated by , )", value = "", placeholder = "KLK2,KLK3,CTSG,MS4A3,CLEC4OP,KDR", height = '400px')),
        column(4, selectInput("s_pathway2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        column(4, uiOutput("s_pathway_3")),
        column(4, numericInput("s_pathway4", label = "p_val_adj", value = 0.05)),
    ),
      box(id = "s_pathway_box2",
        h3("Pathway parameters"), 
        column(4, selectInput("s_pathway5", label = "Organism", choices = c("Human"="org.Hs.eg.db", "Mouse"="org.Mm.eg.db", "Rat"="org.Rn.eg.db"), selected = "org.Hs.eg.db")),
        column(4, selectInput("s_pathway7", label = "pAdjustMethod", choices = c("holm"="holm", "hochberg"="hochberg", "hommel"="hommel", "bonferroni"="bonferroni", "BH"="BH", "BY"="BY", "fdr"="fdr", "none"="none"), selected = "BH")),
        column(4, numericInput("s_pathway8", label = "pvalueCutoff", value = 0.05)),
        column(4, numericInput("s_pathway9", label = "qvalueCutoff", value = 0.2)),
        column(4, numericInput("s_pathway10", label = "Minimal size of genes", value = 10)),
        column(4, numericInput("s_pathway11", label = "Maximal size of genes", value = 500)),
        column(4, selectInput("s_pathway12", label = "Plot type", choices = c("Dot plot"="dotplot", "Bar plot"="barplot", "Net plot"="cnetplot", "upset plot"="upsetplot"), selected = "dotplot")),
        column(4, numericInput("s_pathway13", label = "No. of category to plot", value = 10, min = 1, max = 50)),
        column(4, actionBttn("single_multiple_sample_pathway", "Pathway analysis",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  fluidRow(
    box(id = "s_pathway_box3",status = "warning", width = 11,

        h3("Pathway plot"),
        actionBttn("download_s_pathway1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_pathway1_plot", height = "600px")),
        br(),
        h3("Summary table"), 
        downloadBttn(outputId = "download_s_pathway1_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_pathway1_table")),
        hr(),
    ),
  ),  
  
),
######################################################Menu6#####################################################################   
####################################################GSEA analysis####################################################################
tabPanel(
  "GSEA analysis",
  useShinyjs(),
  fluidRow( 
    box(id = "s_gsea_box0",  width=12,
        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction and Marker Identification step."),
    ),
  ),
  
  fluidRow(   
    box(id = "s_gsea_box1",
        h5("Run the single or multiple samples analysis before start this analysis"),
        h3("Select the input data and cluster(s) for analysis"),
        column(6, selectInput("s_gsea1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering"), selected = "multiple_sample")),
        #column(6, textAreaInput("s_gsea13", label ="Enter your genes (eg: gene names separated by , )", value = "", placeholder = "KLK2,KLK3,CTSG,MS4A3,CLEC4OP,KDR", height = '400px')),
        column(6, selectInput("s_gsea2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        column(6, uiOutput("s_gsea_3")),
        column(6, numericInput("s_gsea4", label = "p_val_adj", value = 0.05)),
    ),
    box(id = "s_gsea_box2",
        h3("GSEA parameters"), 
        column(4, selectInput("s_gsea5", label = "Organism", choices = c("Human"="Homo sapiens", "Mouse"="Mus musculus"), selected = "Homo sapiens")),
        column(4, selectInput("s_gsea6", label = "Category (from MSigDB)", choices = c("Hallmark gene sets (H)"="H", "Positional gene sets (C1)"="C1", "Curated gene sets (C2)"="C2", "Regulatory target gene sets (C3)"="C3", "Computational gene sets (C4)"="C4", "Ontology gene sets (C5)" ="C5", "Oncogenic signature gene sets (C6)"="C6", "Immunologic signature gene sets (C7)"="C7", "Cell type signature gene sets (C8)"="C8"), selected = "H")),
        column(4, selectInput("s_gsea7", label = "ScoreType", choices = c("std"="std", "pos"="pos", "neg"="neg"), selected = "std")),
        column(4, numericInput("s_gsea8", label = "Minimal size of genes", value = 15)),
        column(4, numericInput("s_gsea9", label = "Maximal size of genes", value = 50)),
        column(4, numericInput("s_gsea10", label = "Number of permutations", value = 100)),
        column(4, selectInput("s_gsea11", label = "Plot type", choices = c("GSEA plot"="GSEA_plot", "plotGseaTable"="plotGseaTable", "Bar plot"="barplot"), selected = "GSEA_plot")),
        column(4, numericInput("s_gsea12", label = "No. of significance to plot", value = 10, min = 2, max = 40)),
        column(4, actionBttn("single_multiple_sample_gsea", "GSEA analysis",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ), 
  hr(),
  fluidRow(
    box(id = "s_gsea_box3",status = "warning", width = 12,

        h3("GSEA plot"),
        actionBttn("download_s_gsea1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_gsea1_plot", width="1000px", height="1000px")),
        br(),
        h3("Summary table"), 
        downloadBttn(outputId = "download_s_gsea1_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_gsea1_table")),
        hr(),
    ),
  ),  
  ),
######################################################Menu7#####################################################################   
#########################################Cell-Cell Communication####################################################################
tabPanel(
  "Cell-Cell Communication",
  useShinyjs(),
  fluidRow( 
    box(id = "s_cellchat_box0",  width=12,
        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction and Marker Identification step."),
    ),
  ),
  
  fluidRow(   
    box(id = "s_cellchat_box1",
        h3("Select the input data and celltype method for analysis"),
        column(6, selectInput("s_cellchat1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering"), selected = "multiple_sample")),
        column(6, selectInput("s_cellchat2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        h3("Cell-cell communication parameters  (CellChat)"),
        column(6, selectInput("s_cellchat3", label = "Organism", choices = c("Human"="PPI.human", "Mouse"="PPI.mouse"), selected = "PPI.human")),
        column(6, numericInput("s_cellchat4", label = "Threshold of the percent of cells expressed", value = 0)),
        column(6, numericInput("s_cellchat5", label = "Threshold of Log Fold Change", value = 0)),
        column(6, numericInput("s_cellchat6", label = "Threshold of p-values", value = 0.05)),
        column(6, selectInput("s_cellchat7", label = "Methods for computing the average gene expression per cell group", choices = c("triMean"="triMean", "truncatedMean"="truncatedMean", "thresholdedMean"="thresholdedMean", "median"="median"), selected = "triMean")),
        column(6, numericInput("s_cellchat8", label = "Minmum number of cells required in each cell group for cell-cell communication", value = 10)),
        column(6, numericInput("s_cellchat9", label = "Communication pattern k-value", value = 2)),
        column(6, selectInput("s_cellchat10", label = "Show label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(6, actionBttn("single_multiple_sample_cellchat1", "Cell-Cell communication analysis",  style = "unite",color = "primary", icon = icon("download"))),
        ),

  ),
  hr(),
  fluidRow(
    br(),
    box(id = "s_cellchat_box2",status = "warning", width = 6,
        h3("Interactions plot with counts"),
        #actionBttn("download_s_cellchat1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat1_plot", width="600px",height = "600px")),
        ),
    box(id = "s_cellchat_box3",status = "warning", width = 6,
        h3("Interactions plot with weights/strength"),
        #actionBttn("download_s_cellchat2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat2_plot", width="600px",height = "600px")),
        br(),
    ),
    box(id = "s_cellchat_box4",status = "warning", width = 12,
        h3("Interaction heatmap"),
        #actionBttn("download_s_cellchat3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat3_plot", height = "800px")),
        br(),
    ),
    box(id = "s_cellchat_box5",status = "warning", width = 12,
        h3("Incoming and outgoing signaling patterns"),
        #actionBttn("download_s_cellchat4_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat4_plot", height = "1000px")),
        br(),
        h3("Incoming and Outgoing communication pattern of target and secreting cells"),
        #actionBttn("download_s_cellchat12_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
        withSpinner(plotOutput("s_cellchat12_plot", height = "1000px")),
        br(),
        h3("Interaction table"), 
        downloadBttn(outputId = "download_s_cellchat1_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_cellchat1_table")),
        hr(),
    ),
  ), 
  fluidRow(   
    box(id = "s_cellchat_box6", width = 12,  
        h3("Show all the significant interactions associated with certain signaling pathways"), 
        column(4, uiOutput("s_cellchat_11")),
        column(4, selectInput("s_cellchat12", label = "Show label", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        br(),
        column(4, actionBttn("single_multiple_sample_cellchat2", "Vizualize selected pathway",  style = "unite",color = "primary", icon = icon("download"))),
        br(),
    ),
    hr(),
    br(),
    box(id = "s_cellchat_box7",status = "warning", width = 6,
        h3("Interactions plot (Circle)"),
        #actionBttn("download_s_cellchat5_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat5_plot", width="600px",height = "600px")),
        ),
    box(id = "s_cellchat_box8",status = "warning", width = 6,
        h3("Interactions plot (Chord)"),
        #actionBttn("download_s_cellchat6_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat6_plot", width="600px",height = "600px")),
    ),
    box(id = "s_cellchat_box9",status = "warning", width = 12,
        h3("Interaction heatmap"),
        #actionBttn("download_s_cellchat7_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat7_plot", height = "800px")),
    ),
    box(id = "s_cellchat_box10",status = "warning", width = 12,
        h3("Hierachy plot"),
        #actionBttn("download_s_cellchat11_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat11_plot", height = "800px")),
        h3("Bubble plot"),
        actionBttn("download_s_cellchat8_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat8_plot", height = "800px")),
        ),
    box(id = "s_cellchat_box11",status = "warning", width = 12,
        h3("Network analysis contribution bar plot"),
        actionBttn("download_s_cellchat9_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat9_plot", height = "800px")),
        br(),
        h3("Gene expression plot"),
        actionBttn("download_s_cellchat10_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_cellchat10_plot", height = "1000px")),
        br(),
        h3("Interaction table"), 
        downloadBttn(outputId = "download_s_cellchat2_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_cellchat2_table")),
    ),
  ),
),
######################################################Menu8#####################################################################   
#############################################Trajectory and Pseudotime analysis#################################################
tabPanel(
  "Trajectory and Pseudotime analysis",
  useShinyjs(),
  fluidRow( 
    box(id = "s_trajectory_box0",  width=12,
        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction and Marker Identification step."),
    ),
  ),
  
  fluidRow(   
    box(id = "s_trajectory_box1", 
        h3("Select the input data and annotation method"),
        column(6, selectInput("s_trajectory1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering"), selected = "multiple_sample")),
        column(6, selectInput("s_trajectory2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
        h3("Parameters to Learn Trajectory"), 
        h5("Please make sure you have used UMAP in the clustering steps"),
        column(3, selectInput("s_trajectory3", label = "use_partition", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
        column(3, selectInput("s_trajectory4", label = "close_loop", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(3, selectInput("s_trajectory5", label = "label_groups_by_cluster", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
        column(3, selectInput("s_trajectory6", label = "label_branch_points", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(3, selectInput("s_trajectory7", label = "label_roots", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(3, selectInput("s_trajectory8", label = "label_leaves", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
        br(),
        column(3, actionBttn("single_multiple_sample_trajectory1", "Learn Trajectory",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ),
  hr(),
  fluidRow(
    box(id = "s_trajectory_box2", width=12,
        
        br(),
        h3("Trajectory plot"),
        actionBttn("download_s_trajectory1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("s_trajectory1_plot", width="800px",height = "600px")),
        br(),
        hr(),
    ),
  ),
  fluidRow(   
    box(id = "s_trajectory_box3",                      
        h3("Order cells in pseudotime"), 
        #column(3, selectInput("s_trajectory9", label = "use_partition", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
        column(3, uiOutput("s_trajectory_10")),
        column(3, selectInput("s_trajectory11", label = "label_groups_by_cluster", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(3, selectInput("s_trajectory12", label = "label_branch_points", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(3, selectInput("s_trajectory13", label = "label_roots", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
        column(3, selectInput("s_trajectory14", label = "label_leaves", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
        br(),
        column(3, actionBttn("single_multiple_sample_trajectory2", "Pseudotime analysis",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ),
  fluidRow(
    box(id = "s_trajectory_box4", width=12,
        hr(),
        br(),
        h3("Cells plotted in pseudotime"),
        actionBttn("download_s_trajectory2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("s_trajectory2_plot", width="800px",height = "600px")),
        br(),
        h3("Cells ordered by Seurat cluster and Monocle3 pseudotime"),
        actionBttn("download_s_trajectory3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_trajectory3_plot", height = "600px")),
        br(),
        hr(),
    ),
  ),
  fluidRow(   
    box(id = "s_trajectory_box5",                      
        h3("Find genes that changes function during the pseudotime"), 
        column(3, selectInput("s_trajectory15", label = "neighbor_graph", choices = c("Principal graph" = "principal_graph", "KNN" = "knn"), selected = "principal_graph")),
        column(3, selectInput("s_trajectory16", label = "label_groups_by_cluster", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
        column(3, actionBttn("single_multiple_sample_trajectory3", "Find genes",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ),
  fluidRow(
    box(id = "s_trajectory_box6",status = "warning", width = 11,
        hr(),
        h3("Pseudotime plot"),
        actionBttn("download_s_trajectory4_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        withSpinner(plotOutput("s_trajectory4_plot", width="600px",height = "600px")),
        br(),
        h3("List of genes that changes function during the pseudotime"), 
        downloadBttn(outputId = "download_s_trajectory1_table", label = "Download as csv"),
        withSpinner(dataTableOutput("s_trajectory1_table")),
        hr(),
    ),
  ),
  fluidRow(   
    box(id = "s_trajectory_box7",                      
        h3("Plot the top or user listed genes to see the changes in pseudotime"), 
        column(6, selectInput("s_trajectory17", label = "No. of genes to display", choices = c("1" = 1, "2" = 2, "3"= 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10, "List of gene names" = "gene_name_list"), selected = 5)),
        column(6, textAreaInput("s_trajectory18", label ="Enter your genes for ploting (eg: gene names separated by , )", value = "", placeholder = "KLK2,KLK3,CTSG,MS4A3,CLEC4OP,KDR", height = '400px')),
        column(3, actionBttn("single_multiple_sample_trajectory4", "Plot genes",  style = "unite",color = "primary", icon = icon("download"))),
    ),
  ),
  fluidRow(
    box(id = "s_trajectory_box8", width=12,
        hr(),
        br(),
        h3("Pseudotime plot for the top selected genes"),
        actionBttn("download_s_trajectory5_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
        status = "warning", 
        solidHeader = F,
        withSpinner(plotOutput("s_trajectory5_plot", height = "1000px")),
        br(),
        hr(),
    ),
  ),
), 

######################################################Menu9.1#####################################################################   
####################################################hdWGCNA####################################################################
tabPanel(
  "Co-expression and TF analysis",
  useShinyjs(),
  # sidebarLayout(
  #   sidebarPanel(width=1),
  #   mainPanel(id = "hdWGCNA_tfs_main_menu",
              tabsetPanel(
                type = "tabs", id = "Coexpression_tabsets",
                tabPanel(
                  "Co-expression network analysis",
                  useShinyjs(),
                  fluidRow( 
                    box(id = "s_hdwgcna_box0",  width=12,
                        h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction step."),
                    ),
                  ),
                  
                  fluidRow(   
                    box(id = "s_hdwgcna_box1",
                        h3("Co-expression network analysis using hdWGCNA"),
                        h3("Select the input data and cluster(s) for analysis"),
                        column(6, selectInput("s_hdwgcna1", label = "Input data", choices = c("Output from single or multiple samples" = "multiple_sample", "Output from subclustering" = "multiple_sample_subclustering"), selected = "multiple_sample")),
                        column(6, selectInput("s_hdwgcna2", label = "Select the celltype method", choices = c("Seurat clusters" = "seurat_clusters", "Predicted or own label from previous methods" = "predicted"), selected = "seurat_clusters")),
                        column(6, uiOutput("s_hdwgcna_3")),
                        column(6, selectInput("s_hdwgcna4", label = " select the reduction type", choices = c("umap"="umap", "PCA"="PCA"), selected = "umap")),
                    ),
                    box(id = "s_hdwgcna_box2",
                        h3("Construct metacells"),
                        column(6, numericInput("s_hdwgcna5", label = "Nearest-neighbors parameter (k)", value = 10, min = 1, max = 100)),
                        column(6, numericInput("s_hdwgcna6", label = "Minimum number of cells in a particular grouping to construct metacells", value = 10, min = 1, max = 100)),
                        column(6, numericInput("s_hdwgcna7", label = "Maximum number of shared cells between two metacells", value = 15, min = 1, max = 100)),
                        column(6, numericInput("s_hdwgcna8", label = "Maximum target number of metacells to construct", value = 1000)),
                    ),
                    box(id = "s_hdwgcna_box3",
                        h3("Select soft-power"),
                        column(6, selectInput("s_hdwgcna9", label = "Network Type", choices = c("signed"="signed", "unsigned"="unsigned", "signed hybrid"="signed hybrid"), selected = "signed")),
                    ),
                    box(id = "s_hdwgcna_box4",
                        h3("Module eigengenes and connectivity"),
                        column(6, selectInput("s_hdwgcna10", label = "Scale model", choices = c("linear"="linear", "poisson"="poisson", "negbinom"="negbinom"), selected = "linear")),
                        column(6, selectInput("s_hdwgcna11", label = "Harmonized module eigengenes", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "TRUE")),
                        column(6, numericInput("s_hdwgcna12", label = "Show top N hub genes", value = 10, min = 1, max = 25)),
                    ),
                    box(id = "s_hdwgcna_box5",
                        h3("Module based UMAP Plot"),
                        column(6, numericInput("s_hdwgcna13", label = "No. of hub genes to label in each module", value = 5, min = 1, max = 10)),
                        column(6, selectInput("s_hdwgcna14", label = "Show edges between genes in different modules (grey edges)", choices = c("Yes" = "TRUE", "No" = "FALSE"), selected = "FALSE")),
                        column(6, actionBttn("single_multiple_sample_hdwgcna", "WGCNA analysis",  style = "unite",color = "primary", icon = icon("download"))),
                    ),
                    hr(),
                  ),
                  
                  fluidRow(
                    box(id = "s_hdwgcna_box6",
                        h3("UMAP plot to check the loaded the data is correct"),
                        actionBttn("download_s_hdwgcna1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_hdwgcna1_plot", height = "600px")),
                        h3("Soft power threshold plots"),
                        actionBttn("download_s_hdwgcna2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_hdwgcna2_plot", height = "600px")),
                        h3("Co-expression network plot"),
                        #verbatimTextOutput("text_level_test"),
                        withSpinner(uiOutput("s_hdwgcna3_plot")),
                        h3("Module ranked by eigengene-based connectivity kME"),
                        actionBttn("download_s_hdwgcna4_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_hdwgcna4_plot", height = "600px")),
                        h3("Module feature plots"),
                        actionBttn("download_s_hdwgcna5_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_hdwgcna5_plot", height = "600px")),
                        h3("Module correlagram plot"),
                        withSpinner(uiOutput("s_hdwgcna6_plot")),
                        h3("Module with Seurat’s dot plot"),
                        actionBttn("download_s_hdwgcna7_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_hdwgcna7_plot", height = "600px")),
                        h3("Individual module network plots"),
                        withSpinner(uiOutput("s_hdwgcna8_plot")),
                        h3("UMAP plot for co-expression networks"),
                        withSpinner(uiOutput("s_hdwgcna9_plot")),
                    ),
                    box(id = "s_hdwgcna_box7", width=8,
                        h3("Soft-power threshold table"), 
                        downloadBttn(outputId = "download_s_hdwgcna1_table", label = "Download as csv"),
                        withSpinner(dataTableOutput("s_hdwgcna1_table")),
                    ),
                    box(id = "s_hdwgcna_box8", width=12,
                        h3("Module assignment table"), 
                        downloadBttn(outputId = "download_s_hdwgcna2_table", label = "Download as csv"),
                        withSpinner(dataTableOutput("s_hdwgcna2_table")),
                    ),
                    box(id = "s_hdwgcna_box9", width=5,
                        h3("Top N hub genes"), 
                        downloadBttn(outputId = "download_s_hdwgcna3_table", label = "Download as csv"),
                        withSpinner(dataTableOutput("s_hdwgcna3_table")),
                        br(),
                    ),
                    hr(),
                    box(id = "s_hdwgcna_box10", width=12,
                        downloadBttn(outputId = "s_hdwgcna", label = "Download Object file"),
                        actionLink("link_s_tfrn", actionBttn("SNextButtonLoad", "NEXT STEP (Transcription factor regulatory network analysis)", style = "unite",color = "primary",icon = icon("arrow-right"), size = "md")),
                    ),
                  ),
                ),
######################################################Menu9.2#####################################################################   
####################################################TFs####################################################################
                tabPanel(
                  "Transcription factor regulatory network analysis",
                  # fluidRow( 
                  #   box(id = "s_tfrn_box0",  width=12,
                  #       h3("To begin this analysis, please complete Single or Multiple samples or subclustering analysis until Cell Type Prediction step."),
                  #   ),
                  # ),
                  
                  fluidRow(   
                    box(id = "s_tfrn_box1",
                        h3("Transcription factor regulatory network analysis"),
                        column(6, selectInput("s_tfrn1", label = "Organism", choices = c("Human"="EnsDb.Hsapiens.v86", "Mouse"="EnsDb.Mmusculus.v79"), selected = "EnsDb.Hsapiens.v86")),
                    ),
                    box(id = "s_tfrn_box2",
                        h3("Identify TFs in promoter regions (uses JASPAR 2020 database, Motif scan and XGBOost)"),
                        column(6, numericInput("s_tfrn2", label = "max_depth", value = 1)),
                        column(6, numericInput("s_tfrn3", label = "eta", value = 0.1)),
                        column(6, numericInput("s_tfrn4", label = "alpha", value = 0.5)),
                    ),
                    box(id = "s_tfrn_box3",
                        h3("Define TF Regulons"),
                        column(6, numericInput("s_tfrn5", label = "Threshold for regulatory score", value = 0.01, min=0, max=1)),
                        column(6, numericInput("s_tfrn6", label = "The number of top TFs to keep for each gene", value = 10, min=1, max=20)),
                    ),
                    box(id = "s_tfrn_box4",
                        h3("Calculate regulon expression signatures"),
                        column(6, numericInput("s_tfrn7", label = "Positive regulon score thresold", value = 0.05, min = 0, max = 1)),
                        column(6, numericInput("s_tfrn8", label = "Negative regulon score thresold", value = -0.05, min = 0, max = -1)),
                        column(6, actionBttn("single_multiple_sample_tfrn1", "Transcription factor analysis",  style = "unite",color = "primary", icon = icon("download"))),
                    ),
                  ),
                  hr(),
                  fluidRow(
                    box(id = "s_tfrn_box5", width=12,
                        h3("Module regulatory network plot (Positive)"),
                        actionBttn("download_s_tfrn1_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn1_plot", width="800px", height = "600px")),
                        
                        h3("Module regulatory network plot (Negative)"),
                        actionBttn("download_s_tfrn2_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn2_plot", width="800px", height = "600px")),
                        
                        h3("Module regulatory network plot (Both)"),
                        actionBttn("download_s_tfrn3_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
                        withSpinner(plotOutput("s_tfrn3_plot", width="800px",height = "600px")),
                        
                        h3("Module regulatory network plot (Module UMAP)"),
                        actionBttn("download_s_tfrn4_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn4_plot", width="800px", height = "600px")),
                        
                        h3("TF network table"), 
                        downloadBttn(outputId = "download_s_tfrn1_table", label = "Download as csv"),
                        withSpinner(dataTableOutput("s_tfrn1_table")),
                        hr(),
                    ),
                  ),
#####################submenu#############
                  fluidRow(   
                    box(id = "s_tfrn_box6", width = 12,  
                        h3("Select a TF of interest"), 
                        column(6, uiOutput("s_tfrn_11")),
                        h3("Bar plot parameter"), 
                        column(6, numericInput("s_tfrn12", label = "Number of top and bottom target genes", value = 10, min = 1 , max = 25)),
                        h3("Network plot parameter"), 
                        column(4, selectInput("s_tfrn13", label = "Attribute to color the network edges", choices = c("Cor" = "Cor", "Gain" = "Gain"), selected = "Cor")),
                        column(4, selectInput("s_tfrn14", label = "Number of layers to extend the TF network", choices = c("Primary" = 1, "Primary and secondary" = 2, "Primary, secondary and tertiary" = 3), selected = 2)),
                        column(4, actionBttn("single_multiple_sample_tfrn2", "Plot networks",  style = "unite",color = "primary", icon = icon("download"))),
                        br(),
                       
                    ),
                  ),
                   hr(),
                  fluidRow(
                    box(id = "s_tfrn_box7", width=12,
                        h3("Feature plot of selected TF"),
                        actionBttn("download_s_tfrn11_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn11_plot", width="800px", height = "600px")),
                        
                        h3("Top target genes within TF regulons"),
                        actionBttn("download_s_tfrn12_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn12_plot", width="800px", height = "600px")),
                        
                        h3("TF network plot (Positive)"),
                        actionBttn("download_s_tfrn13_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")),
                        withSpinner(plotOutput("s_tfrn13_plot", width="800px",height = "600px")),
                        
                        h3("TF network plot (Negative)"),
                        actionBttn("download_s_tfrn14_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn14_plot", width="800px", height = "600px")),
                        
                        h3("TF network plot (Both)"),
                        actionBttn("download_s_tfrn15_plot", "Download plot", style = "unite",color = "primary", icon = icon("download")), 
                        withSpinner(plotOutput("s_tfrn15_plot", width="800px", height = "800px")),
                        hr(),
                        br(),
                        downloadBttn(outputId = "s_tfrn", label = "Download Object file"),
                    ),
                  ),
                ),
              ),
  #   ),
  # ),
),
######################################################Menu10#####################################################################   
###################################################Manual####################################################################
tabPanel(
  "Manual",
  mainPanel(
    HTML("<h1 align='center'>Single Cell RNA Data Analysis and Visualization</h1>
<hr>
<h3>This section will introduce how to prepare input files:</h3>
<p><strong>Supported Input Formats:</strong></p>
<ol>
  <li><strong>H5 Files (Cell Ranger Output)</strong></li>
<ul>
  <li>Cell Ranger file: filtered_feature_bc_matrix.h5.</li>
  <li>Rename it to SAMPLE_NAME.h5 for proper identification. </li>
</ul>

  <li><strong>Cell Ranger Matrix Files</strong></li>
<ul>
  <li>Cell Ranger files: matrix.mtx.gz, feature.tsv.gz, barcode.tsv.gz.</li>
  <li>Rename as SAMPLE_NAME_matrix.mtx.gz, SAMPLE_NAME_features.tsv.gz, and SAMPLE_NAME_barcodes.tsv.gz, and upload together as a set.</li>
</ul>
  <li><strong>Seurat Objects</strong></li>
<ul>
<li>Format: filename.rds (Seurat object). The orig.ident attribute should match the sample name(s).</li>
</ul>
  <li><strong>Matrix count file</strong></li>
<ul>
<li>Format: Filename.txt with rows as genes and columns as sample_cellID.</li>
</ul>
</ol>
<p><strong>Data Size and Handling:</strong></p>
<ul>
  <li>The tool can handle scRNA-seq data up to 3GB in the specified formats.</li>
  <li>Supports analysis of single or multiple samples, including up to six sample groups.</li>
  <li>After data upload, users can proceed with the analysis through a step-by-step workflow for the 1st Module, with the 'Next Step' button guiding users through each tab in the process.</li>
  <li>Once the single or multiple analysis is completed, users can analysis as per their need, there is no steps involved further</li>
</ul>
<p><strong>Output and Visualizations:</strong></p>
<ul>
  <li><strong>High-Quality Plot Download:</strong> Users can download plots in seven formats: JPG, TIFF, PDF, SVG, BMP, EPS, and PS. However, a few specific plots, such as those requiring exceptionally high detail or complex rendering (e.g., network graphs or high-resolution heatmaps), are only available as PDF files to preserve their quality and detail.</li>
  <li><strong>Summary Tables:</strong> Tables are displayed using the DT package. Users can visualize up to 100 rows (default is 10) and download the entire table as a CSV file.</li>
  <li><strong>Download Seurat Object:</strong> In single or multiple sample analyses, users can download the processed results as an RDS file (Seurat Object). </li>
</ul>
<p><strong>Example Datasets:</strong></p>
<p>To ensure seamless analysis and reproducibility, ScRDAVis includes one reference dataset for each input format, sourced from NCBI, which has been pre-tested with the tool. These datasets allow users to explore the tool's functionalities and understand the analysis workflow effectively.</p>
<ul>
  <li><strong>H5 File:</strong> <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE271107'>GSE271107</a></li>
  <li><strong>Matrix Files:</strong> <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277476'>GSE277476</a></li>
  <li><strong>Seurat Object:</strong> <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE250488'>GSE250488</a></li>
  <li><strong>Matrix count file:</strong> <a href='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155953'>GSE155953</a></li>
</ul>

<hr>
<h3>Step-by-Step Approach for User Interaction</h3>
<h3>1. Single or Multiple samples analysis</h3>
	<h3>1.1 Stats</h3>
<ul>
    <li><strong>Upload and Adjust Parameters:</strong>
        <ul>
            <li>Minimum cell expression per gene: Define the minimum number of cells that should express each gene.</li>
            <li>Minimum gene expression per cell: Set the minimum number of genes each cell should express.</li>
        </ul>
		
    </li>
	<li><strong>Execution:</strong>
<ul>
<li>Click the Submit button to run the analysis based on selected parameters.</li>
</ul>
</li>
    <li><strong>Output:</strong>
        <ul>
            <li>QC plots: Quality metrics before filtering <strong>(Fig. 1.1a)</strong>.</li>
            <li>Feature-Feature Relationships Plot <strong>(Fig. 1.1b)</strong> </li>
          <li>Sample(s) cell counts table <strong>(Fig. 1.1c)</strong></li>
            <img src='images/1.1.jpg' width='800' height='500' alt=''/>
        </ul>
    </li>
</ul>

<h3>1.2. Sample Groups and QC Filtering</h3>
<ul>
    <li><strong>Assign sample(s) group(s):</strong><br>
        
            Choose the number of group(s) based on sample grouping
		<ul>
            <li>If only one sample: Select 1 group.</li>
<li>For multiple samples in a single condition: Select 1 group.</li>
            <li>For multiple samples in different conditions: Select up to 6 groups.</li>
        </ul>
    </li>
    <li><strong>Define Filtering Parameters:</strong>
        <ul>
            <li>Set thresholds to remove low-quality cells and genes, and optionally filter out cells with high mitochondrial gene content.</li>
			 </ul>
	</li>
			<li><strong>Execution:</strong>
<ul>
	<li>Update Filtered Data: Click to apply filters and update the dataset.</li>
</ul>
</li>
            
        
    <li><strong>Outputs after Filtering:</strong>
        <ul>
            <li>QC matrics and Bar Plot <strong>(Fig. 1.2a-d)</strong></li>
            <li>Summary Table: Filtered cell count data <strong>(Fig. 1.2e,f)</strong></li>
			<img src='images/1.2.jpg' width='800' height='650' alt=''/>
        </ul>
    </li>
</ul>

<h3>1.3. Normalization and PCA Analysis</h3>
<ul>
    <li><strong>Normalization Methods:</strong>
        <ul>
            <li>LogNormalize: Adjusts for sequencing depth or read count differences.
                <ul>
                    <li>Additional options: Scale Factor, Variable Genes Detection methods (vst, mvp, disp), Top Variable Features to use.</li>
                </ul>
            </li>
            <li>SCT (SCTransform): Uses regularized negative binomial regression for clustering and differential expression.</li>
        </ul>
    </li>
    <li><strong>PCA Settings:</strong>
        <ul>
            <li>Choose the PCA Dimensions (typically between 1-50) for analysis.</li>
        </ul>
    </li>
	<li><strong>Execution:</strong>
<ul>
	<li>Click  Submit button to start the analysis.</li>
</ul>
</li>
    <li><strong>Outputs:</strong>
        <ul>
            <li>PCA Heatmap <strong>(Fig. 1.3a)</strong></li>
            <li>Elbow Plot <strong>(Fig. 1.3b)</strong></li>
            <li>PCA Plot (sample-wise or group-wise) <strong>(Fig. 1.3c,d)</strong></li>
			<img src='images/1.3.jpg' width='800' height='600' alt=''/>
        </ul>
    </li>
</ul>

<h3>1.4. Clustering</h3>
<ul>
    <li><strong>Clustering Step:</strong>
        <ul>
            <li>Find Neighbors: The users selects the dimensions to use (PCA, integrated dimensions, etc.) and k-nearest neighbors.</li>
            <li>Clustering Algorithm: Select between Louvain, SLM or Leiden algorithms for clustering.</li>
            <li>Resolution Control: The users can adjust the resolution (0.1 to 1) parameter to control the granularity of clusters.</li>
        </ul>
    </li>
    <li><strong>Dimension Reduction:</strong>
        <ul>
            <li>o Choose between UMAP or t-SNE for dimensionality reduction.
                <ul>
                    <li>For UMAP: Users can adjust parameters like min.dist, k-nearest-neighbours, and the number of dimensions.</li>
                    <li>For t-SNE: Users can adjust the number of dimensions</li>
                </ul>
            </li>
        </ul>
    </li>
    <li><strong>Integration Method:</strong>
		 <ul>
    <li>Select integration method(s): CCA, RPCA, Harmony, or JointPCA. These methods will allow users to handle dataset complexities and integrate data from multiple samples.</li>
	  </ul>
  </li>
	<li><strong>Execution:</strong>
<ul>
	<li>Click  Submit button to start the analysis.</li>
</ul>
</li>
<li><strong>Visualize and Compare:</strong>
      <ul>
          <li>Display UMAP or t-SNE plots with clustering labels and sample/condition overlays <strong>(Fig. 1.4a,c,e)</strong>.</li>
          <li>Bar charts (Fig. 4b,d,f) and tables show cell counts per cluster and per sample/condition <strong>(Fig. 1.4g-i)</strong>.</li>
		  <img src='images/1.4.jpg' width='800' height='700' alt=''/>
      </ul>
  </li>
</ul>

<h3>1.5. Remove Doublets</h3>
<ul>
    <li><strong>Doublet Detection:</strong> 
<ul>
    <li>Uses DoubletFinder to identify potential doublets, helping identify cells that may contain RNA from more than one original cell.</li>
    <li>Start with an estimated doublet rate of 7.5%-10%&nbsp; (0.075 to 0.1) of total cell count.</li>
	</ul>
    </li>
	<li><strong>Execution:</strong>
<ul>
	<li>Click Detect Doublet button to start the analysis.</li>
</ul>
</li>
    <li><strong>Outputs:</strong>
        <ul>
            <li>UMAP or t-SNE Plot: Shows singlet and doublet cells, color-coded by cluster, sample, or group <strong>(Fig. 1.5a-d)</strong>.</li>
            <li>Summary Table: Counts of singlet and doublet cells <strong>(Fig. 1.5e)</strong>.</li>
        </ul>
    </li>
  <li><strong>Remove or keep doublets:</strong>
	  <ul>
    <li>Users can choose to Keep or Remove Doublets. Updated plots and tables reflect the selection <strong>(Fig. 1.5f-o)</strong>.</li>
		<img src='images/1.5.jpg' width='800' height='700' alt=''/>
	</ul>
  </li>
</ul>

<h3>1.6. Marker Identification</h3>
<ul>
    <li><strong>Identify markers in all clusters (FindAllMarkers):</strong>
        <ul>
            <li>Customizable parameters:
				<ul>
          <li>Minimum cell percentage (min.pct) to specify the minimum fraction of cells in which a gene is expressed.</li>
          <li>Log fold-change threshold (logfc.threshold) to filter markers based on expression magnitude.</li>
          <li>Statistical test options (test.use), including Wilcoxon rank sum (wilcox), Wilcoxon-Limma hybrid (wilcox_limma), binomial (bimod), ROC, t-test, likelihood ratio test (LR), and MAST.</li>
          <li>Positive markers only (only.pos), with options for yes or no, to focus on upregulated genes in the target cluster.</li>
            <li>When using SCTransform for normalization, out tool uses PrepSCTFindMarkers preps the data for accurate differential testing by adjusting the SCT assay, making results more reliable for FindMarkers and FindAllMarkers.</li>
				</ul>
			</li>
<li>Output: Heatmap of the top 5 genes per cluster <strong>(Fig. 1.6a)</strong>, helping users visualize the distinguishing genes for each cluster and Summary table of markers or expressed genes <strong>(Fig. 1.6b)</strong>.</li>
	<img src='images/1.6.jpg' width='800' height='750' alt=''/>
        </ul>
    </li>
    <li><strong>Marker identification in one specific cluster or between two clusters (FindMarkers):</strong> 
  <ul>
    <li>Identifies markers for one cluster against another or against all other clusters.</li>
    <li>Includes all the customizable parameters noted above, enabling targeted cluster comparison with refined criteria.</li>
    <li>Output: A table format displaying the expressed genes for the specified clusters, ideal for in-depth comparisons.</li>
	  </ul>
</li>
    <li><strong>Conserved marker identification for one vs. all cluster or between two clusters (FindConservedMarkers):</strong>
        <ul>
            <li>Finds markers conserved across groups (e.g., conditions) for a cluster, or conserved markers between two specific clusters.</li>
          <li>Utilizes the same customizable parameters for consistency across comparisons.</li>
            <li>Output: A table format with expressed genes, providing insights into markers consistently expressed across groups or clusters.</li>
        </ul>
    </li>
</ul>

<h3>1.7. Cell Type Prediction</h3>
<ul>
  <li><strong>ScType:</strong></li>
       <ul>
         <li>Predefined Tissue Types: Users can select from 15 tissue types, including: Adrenal, Brain, Eye, Heart, Immune, Intestine, Kidney, Liver, Lung, Muscle, Pancreas, Placenta, Spleen, Stomach, Thymus.</li>
         <li>Tissue Classification: Automatically classifies the cells based on the selected tissue type.</li>
   </ul>
    <li><strong>SingleR:</strong></li>
      <ul>
         <li>Reference Datasets: Users can use reference datasets such as: Human Primary Cell Atlas, Blueprint/ENCODE, Mouse RNA-seq, Immunological Genome Project, Database of Immune Cell Expression/eQTLs/Epigenomics, Novershtern Hematopoietic data, Monaco immune data.</li>
         <li>Prediction: Predicts the cell types based on these well-known reference datasets.</li>
     </ul>
    <li><strong>GPTCelltype:</strong></li>
      <ul>
         <li>GPT Models: Utilizes various GPT models, including: GPT-4, GPT-4-turbo, GPT-4o-mini, GPT-4o, ChatGPT-4o-latest, GPT-3.5-turbo, GPT-3.5-turbo.</li>
         <li>Gene Requirements: Requires a minimum number of top genes for accurate prediction.</li>
         <li>Availability: Available via the web platform. To use it locally, users need to update their API key by setting Sys.setenv(OPENAI_API_KEY = 'your_openai_API_key') in the global.R file.</li>
      </ul>
  <li><strong>Own Cell Labels:</strong></li>
      <ul>
         <li>User-Defined Labels: Users can manually input their own cell type labels for each cluster.</li>
         <li>Cluster Grouping: If multiple clusters need the same label, users should provide the same label name for those clusters.</li>
      </ul>
      
    <li><strong>UMAP/t-SNE Labels:</strong>
		<ul>
         <li>Label display options: Users can choose to show or hide cell type labels in the UMAP or t-SNE plots.</li>
      </ul>
	</li>
	<li><strong>Execution:</strong>
<ul>
	<li>Click Detect cell type button to start the analysis.</li>
</ul>
</li>
   <li><strong>Output:</strong></li>
      <ul>
         <li>Plot: Generates an image plot showing the predicted cell types <strong>(Fig. 1.7a,b)</strong>.</li>
         <li>Summary Table: Provides a summary table with the predicted cell types and associated scores <strong>(Fig. 1.7c)</strong>.</li>
		  <img src='images/1.7.jpg' width='800' height='650' alt=''/>
      </ul>
</ul>

<h3>1.8. Cluster-based Plots</h3>
<ul>
    <li><strong>Gene Selection:</strong>
        <ul>
            <li>Top Genes: Users can select the top features or genes (from 2 to 10).</li>
            <li>Custom Genes: Users may also input custom gene names by selecting them from a drop-down menu (list of genes) and entering the desired gene names as a comma-separated list.</li>
        </ul>
    </li>
    <li><strong>Plot Types:</strong>
		<ul>
    <li>Multiple visualization formats are available, including: Dot Plot, Violin Plot, Ridge Plot, Feature Plot.</li>
		</ul>
  </li>
<li><strong>Grouping and Splitting:</strong>
	<ul>
	<li>Group by: Users can organize the data by Seurat clusters or labels generated from previous cell type prediction steps.</li>
	<li>Split by: If multiple samples are present, plots can be split by condition or sample to compare expression patterns across groups.</li>
	</ul>
  </li>
	<li><strong>Execution:</strong>
<ul>
	<li>Click Generate plots button to start the analysis.</li>
</ul>
</li>
    <li><strong>Output:</strong>
	<ul>
	<li>Plot: The user receives one of the chosen plot formats (violin plot, dot plot, feature plot or ridge plot) <strong>(Fig. 1.8a-d)</strong>.</li>
	<li>Summary Tables: The tool generates tables showing marker gene cell counts and cell proportions, providing an additional layer of quantitative insight <strong>(Fig. 1.8e)</strong>.</li>
	  <img src='images/1.8.jpg' width='800' height='600' alt=''/>
		</ul>
	</li>
</ul>

<h3>1.9. Condition-based Analysis</h3>
<ul>
    <li><strong>Group Selection:</strong> 
	<ul>
	<li>Users can compare gene expression between two conditions by selecting one group per dropdown menu.</li>
	</ul>
	</li>
    <li><strong>Customizable Parameters:</strong><ul>
	<li>Minimum Cell Percentage (min.pct): Sets the minimum fraction of cells in which a gene must be expressed.</li>
	<li>Log Fold-Change Threshold (logfc.threshold): Filters markers by expression magnitude.</li>
	<li>Statistical Tests (test.use): Users can choose from various methods, including Wilcoxon rank sum, Wilcoxon-Limma hybrid, binomial, ROC, t-test, likelihood ratio test, and MAST.</li>
	<li>Positive Markers Only (only.pos): Option to display only upregulated genes in the target cluster.</li>
	</ul>
	</li>
    <li><strong>Visualization Options:</strong>
	<ul>
	<li>Multiple formats are available, including: Dot Plot, Violin Plot, Ridge Plot, Feature Plot, Volcano Plot</li>
	</ul>
	</li>
    <li><strong>Grouping:</strong><ul>
	<li>Group By: Users can group data by Seurat clusters or predicted cell type labels.</li>
	<li>Number of Features: Allows display of a specific number of up- and down-regulated genes (e.g., 15).</li>
	</ul>
	</li>
	<li><strong>Execution:</strong>
<ul>
	<li>Click Submit button to start the analysis.</li>
</ul>
</li>
    <li><strong>Output:</strong><ul>
	<li>Plot: The users receives the chosen plot type, providing visual comparison <strong>(Fig. 1.9a-e)</strong>.</li>
	<li>Summary Tables: Table contains the differentially expressed genes between the slected groups <strong>(Fig. 1.9f)</strong>.</li>
	<li>This setup enables users to conduct detailed comparisons between conditions, facilitating insights into differential gene expression and cellular responses.</li>
	<img src='images/1.9.jpg' width='800' height='700' alt=''/>
    </ul>
	</li>
</ul>
	<hr>
	<h3>Subclustering</h3>
	In ScRDAVis, users can further explore specific clusters of interest by performing subclustering analysis. This feature allows for a more granular examination of cell populations within one or multiple clusters, based on the user’s selection. Similar output were generated as like as above for the selected cluster(s) or cell type(s).

<h3>2.1. Cluster Selection:</h3>
<ul>
<li>Users can choose one or multiple clusters for subclustering.</li>
<li>Clusters can be selected based on Seurat clusters or previously predicted annotation labels.</li>
</ul>
	
<h3>2.2. Subclustering Analysis Steps:</h3>
<ul>
<li>The subclustering process mirrors the main workflow, with dedicated tabs for each stage, allowing users to perform the following analyses on the selected cluster(s):
<ul>
<li><strong>Cell Stats:</strong> Overview of cell metrics within the selected clusters, including minimum gene and cell expression thresholds.</li>
<li><strong>Normalization and PCA Analysis: </strong>Options to normalize and vizualize the  PCA data for secific cluster(s). (use the same method used in the above menu)</li>
<li><strong>Clustering:</strong> Allows users to re-cluster cells within the subclusters, providing insights into finer subpopulations.</li>
<li><strong>Marker Identification:</strong> Users can identify markers specific to subclusters, with options to customize parameters for marker detection.</li>
<li><strong>Cell Type Prediction:</strong> Provides options to predict cell types within the selected subclusters using ScType, SingleR, GPTCelltype, or custom labels.</li>
<li><strong>Cluster-Based Plots:</strong> Users can visualize gene expression within subclusters through Dot, Violin, Ridge, or Feature plots.</li>
<li><strong>Condition-Based Analysis:</strong> Enables differential expression comparisons within subclusters, providing insights into condition-specific gene expression patterns.</li>
</ul>
</li>
</ul>
<hr>
<h3>3. Correlation Network Analysis</h3>
ScRDAVis includes Cluster-Based Correlation Analysis using the genesorteR package. This feature helps users explore relationships and interactions among genes within specific clusters by calculating pairwise correlations.
<ul>
<li><strong>Prerequisites:</strong>
<ul>
<li>Correlation Network Analysis becomes available after completing single or multiple samples analysis or subclustering analysis up to cell type prediction.</li>
<li>Users can choose to conduct  analysis on: Seurat clusters or predicted cell type labels from single, multiple, or subcluster analyses.&nbsp;</li>
</ul>
</li>
<li><strong>Correlation Methods:</strong>
<ul>
<li>Pearson</li>
<li>Spearman</li>
<li>Kendall</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
<li>Click the Cluster correlation network button to run the analysis based on selected parameters.</li>
</ul>
</li>
<li><strong>Output:</strong>
<ul>
<li>Correlation Heatmap: Displays the correlation values between genes within clusters in a matrix format <strong>(Fig. 3a)</strong>.</li>
<li>Correlation Network Plot: Depicts the relationships between genes as a network, highlighting strongly correlated pairs <strong>(Fig. 3b)</strong>.</li>
<li>Summary Table: With the complete correlation matrix for detailed analysis <strong>(Fig. 3c)</strong></li>
	<img src='images/3.1.jpg' width='800' height='600' alt=''/>
</ul>
</li>
</ul>
This analysis provides a deeper understanding of gene co-expression and interaction patterns within clusters, aiding in the identification of significant biological relationships.
<hr>
<h3>4. GO Term Analysis</h3>
	ScRDAVis provides integrated Gene Ontology (GO) term analysis using the clusterProfiler package, enabling users to explore biological functions, molecular mechanisms, and cellular components related to gene expression patterns in single, multiple, or subcluster analyses. Here’s how users can conduct GO analysis:.

<ul>
<li><strong>Prerequisites:</strong>
  <ul>
    <li>GO analysis becomes available after completing single or multiple samples analysis or subclustering analysis up to cell type prediction.</li>
  </ul>
</li>
<li><strong>Input Options:</strong><br>
Users can choose to conduct GO analysis on:
<ul>
<li>Seurat clusters or predicted cell type labels from single, multiple, or subcluster analyses.</li>
<li>Users can use one or multiple clusters at a time, with an adjustable parameter (p_val_adj < 0.05) for significant results.</li>
<li>A custom list of genes: Users can manually enter gene names (comma-separated) to investigate GO terms for genes of specific interest.</li>
</ul>
</li>
<li><strong>Organisms Supported for GO Term Mapping:</strong><br>
ScRDAVis supports GO analysis for five organisms, mapping gene IDs to gene symbols:
<ul>
<li>Human: org.Hs.eg.db</li>
<li>Mouse: org.Mm.eg.db</li>
<li>Rat: org.Mmu.eg.db</li>
<li>Pig: org.Ss.eg.db</li>
<li>Rhesus: org.Rn.eg.db</li>
</ul>
</li>
<li><strong>GO Term Analysis Parameters:</strong><br>
Ontology Method: Users can choose to focus on specific biological aspects or all three:
<ul>

<li>Biological Process (BP)</li>
<li>Molecular Function (MF)</li>
<li>Cellular Component (CC)</li>
<li>All: To analyze across all three categories.</li>
</ul>
</li>
<li><strong>Adjustable Parameters:</strong>
<ul>
<li>pAdjustMethod: Select the method to adjust for multiple testing.</li>
<li>pvalueCutoff: Set a cutoff for p-values.</li>
<li>qvalueCutoff: Define a q-value threshold for significance.</li>
<li>Minimum Size of Genes: Minimum number of genes required in a GO term.</li>
<li>Maximum Size of Genes: Maximum number of genes in a GO term.</li>
<li>Plot Type: Choose a visualization format (Dot Plot, Bar Plot, Net and UpSetPlot).</li>
<li>Number of Categories to Plot: Select the number of categories to display (1 to 50).</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
<li>Click the GO Term button to run the analysis based on selected parameters.</li>
</ul>
</li>
<li><strong>Output</strong>
<ul>
The GO Term analysis provides:
<li>Plots: Dot plot, Bar plot, UpSet plot and Network plot for the selected ontology categories <strong>(Fig. 4a-d)</strong>.</li>
<li>Summary Table: A downloadable table summarizing the GO terms, adjusted p-values, and other relevant metrics, allowing users to interpret and visualize biological insights <strong>(Fig. 4e)</strong>.</li>
	<img src='images/4.1.jpg' width='800' height='600' alt=''/>
</ul>
</li>
</ul>
This GO term analysis feature in ScRDAVis provides users with an accessible, visually informative, and comprehensive view of gene functionality across clusters and conditions, enabling enhanced biological interpretation of scRNA-seq data.
<hr>
<h3>5. Pathway Analysis</h3>
ScRDAVis offers pathway analysis through KEGG and Reactome databases using the clusterProfiler and ReactomePA packages. Users can gain insights into biological pathways associated with specific gene expression profiles from single, multiple, or subcluster analyses.
<ul>
<li><strong>Prerequisites:</strong>
<ul>
<li>Pathway analysis becomes available after completing single or multiple samples analysis or subclustering analysis up to cell type prediction.</li>
</ul>
</li>
<li><strong>Input Options:</strong><br>
Pathway analysis can be performed on:
<ul>
<li>Seurat clusters or predicted cell type labels from single, multiple, or subcluster analyses.</li>
<li>One or multiple clusters simultaneously, with results filtered by an adjusted p-value (p_val_adj < 0.05).</li>
<li>A custom list of genes: Users can input specific gene names (comma-separated) to focus on pathways for genes of interest.</li>
</ul>
</li>
<li><strong>Organisms Supported for Pathway Mapping:</strong><br>
ScRDAVis enables pathway mapping for multiple organisms:
<ul>
<li>KEGG Pathways: Supports human (org.Hs.eg.db), mouse (org.Mm.eg.db), and rat (org.Mmu.eg.db) for mapping gene IDs to symbols.</li>
<li>Reactome Pathways: Available for human, mouse, and rat.</li>
</ul>
</li>
<li><strong>Pathway Analysis Parameters:</strong>
<ul>
<li>pAdjustMethod: Choose a method for multiple testing correction.</li>
<li>pvalueCutoff: Set a threshold for p-values.</li>
<li>qvalueCutoff: Define a q-value cutoff for pathway significance.</li>
<li>Minimum Size of Genes: Minimum gene count per pathway.</li>
<li>Maximum Size of Genes: Maximum gene count per pathway.</li>
<li>Plot Type: Choose visualization format (Dot Plot, Bar Plot, Net and UpSetPlot Plot).</li>
<li>Number of Pathways to Plot: Select the number of pathways to display (1 to 50).</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
<li>Click the Pathway Analysis button to run the analysis with the selected parameters.</li>
</ul>
</li>
<li><strong>Output</strong><br>
Pathway analysis results include:
<ul>
<li>Visualizations: Dot plot, Bar plot, UpSet plot and Network plot, showcasing significant pathways <strong>(Fig. 5a-d)</strong>.</li>
<li>Summary Table: A downloadable table with pathway details, adjusted p-values, and other metrics for further exploration and interpretation <strong>(Fig. 5e)</strong>.</li>
	<img src='images/5.1.jpg' width='800' height='600' alt=''/>
</ul>
</li>
</ul>
The pathway analysis functionality in ScRDAVis helps users understand the biological processes and signaling pathways linked to gene expression profiles across clusters and conditions, providing a deep functional understanding of their ScRNA-seq data.
<hr>
<h3>6. GSEA Analysis</h3>
The Gene Set Enrichment Analysis (GSEA) feature in ScRDAVis leverages the fgsea package to identify enriched pathways using ranked gene lists, such as those generated from differential expression analysis. This allows users to assess pathway-level expression changes and gain insights into functional changes across clusters or conditions.
<ul>
<li><strong>Prerequisites</strong>
<ul>
<li>GSEA analysis can be conducted following single or multiple samples analysis or subclustering analysis  up to cell type prediction.</li>
</ul>
</li>
<li><strong>Input Options</strong><br>
GSEA analysis can be performed on:
<ul>
<li>Seurat clusters or predicted cell type labels derived from single, multiple, or subcluster analyses.</li>
<li>Single or multiple clusters simultaneously, with results filtered by an adjusted p-value (p_val_adj < 0.05).</li>
</ul>
</li>
<li><strong>Organisms and Gene Sets Supported</strong>
<ul>
<li>Organisms: Human and mouse gene ID mapping.</li>
<li>Gene Set Categories: Using the msigdbr package, which provides gene sets compatible with fgsea from the Molecular Signatures Database (MSigDB). Available categories include:
<ul>
<li>Hallmark Gene Sets (H)</li>
<li>Positional Gene Sets (C1)</li>
<li>Curated Gene Sets (C2)</li>
<li>Regulatory Target Gene Sets (C3)</li>
<li>Computational Gene Sets (C4)</li>
<li>Ontology Gene Sets (C5)</li>
<li>Oncogenic Signature Gene Sets (C6)</li>
<li>Immunologic Signature Gene Sets (C7)</li>
<li>Cell Type Signature Gene Sets (C8)</li>
</ul>
</li>
</ul>
</li>
<li><strong>GSEA Analysis Parameters:</strong>
<ul>
<li>scoreType: Define the scoring method for pathway enrichment.</li>
<li>Minimal Size of Genes: Minimum number of genes in a gene set.</li>
<li>Maximal Size of Genes: Maximum number of genes in a gene set.</li>
<li>Number of Permutations: Control the precision of p-value calculations.</li>
<li>Plot Type: Choose visualization format (GSEA Plot, PlotGseaTable, Bar Plot).</li>
<li>Number of Significant Pathways to Plot: Select the number of pathways to display (1 to 40).</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
<li>Click the GSEA Analysis button to run the analysis with the selected parameters.</li>
</ul>
</li>
<li><strong>Output</strong><br>
The GSEA analysis provides:
<ul>
<li>Visualizations:
<ul>
<li>GSEA Plot: Displays the enrichment score curve <strong>(Fig. 6a)</strong>.</li>
<li>PlotGseaTable: Shows enriched pathways and their enrichment scores <strong>(Fig. 6b)</strong>.</li>
<li>Bar Plot: Highlights top significant pathways <strong>(Fig. 6c)</strong>.</li>
</ul>
</li>
<li>Summary Table: A downloadable table of enriched pathways, adjusted p-values, and scores. If the users selects the top 10 significant pathways, the tool displays the top 5 upregulated and top 5 downregulated pathways <strong>(Fig. 6d)</strong>.</li>
	<img src='images/6.1.jpg' width='800' height='550' alt=''/>
</ul>
</li>
</ul>
GSEA analysis in ScRDAVis offers a powerful method for understanding pathway-level dynamics, supporting biological interpretation of ScRNA-seq data through visual and quantitative assessments of enriched pathways.
<hr>
<h3>7. Cell-Cell Communication Analysis</h3>	
ScRDAVis integrates CellChat to enable users to analyze cell-cell communication within single or multiple samples, as well as for subclusters. This analysis identifies potential ligand-receptor interactions, allowing users to explore how different cell types or clusters communicate based on gene expression patterns.
<ul>
<li><strong>Input Options</strong>
<ul>
<li>Source of Input: Users can analyze cell-cell communication using Seurat clusters or predicted cell type labels generated from single, multiple, or subcluster analysis.</li>
<li>Organisms Supported: Human and mouse datasets are available for ligand-receptor interaction mapping.</li>
</ul>
</li>
</ul>
<h3>7.1. Parameters for Cell-Cell Communication</h3>
<ul>
<li><strong>Identify Over-Expressed Genes:</strong>
<ul>
<li>Threshold of Cell Expression Percentage: Minimum percentage of cells expressing the genes.</li>
<li>Log Fold Change Threshold: Minimum log fold-change required for genes to be considered over-expressed.</li>
<li>p-Value Threshold: Statistical significance threshold.</li>
</ul>
</li>
<li><strong>Compute Communication Probability:</strong>
<ul>
<li>Expression Method: Choose how to compute the average expression per cell group (options: triMean, truncatedMean, thresholdedMean, median).</li>
</ul>
</li>
<li><strong>Filter Communication:</strong>
<ul>
<li>Minimum Cell Requirement: Minimum number of cells needed in each cell group to analyze cell-cell communication.</li>
</ul>
</li>
<li><strong>Communication Pattern Identification:</strong>
<ul>
<li>Pattern k-Value: Defines the number of communication patterns to identify.</li>
</ul>
</li>
<li><strong>Label Option:</strong>
<ul>
<li>Show or hide labels in plots.</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
	<li>Click to Cell-Cell communication analysis button to start the analysis.</li>
</ul>
</li>
<li><strong>Output for Cell-Cell Communication Analysis</strong><br>
The analysis generates the following visual outputs:
<ul>
<li>Interaction Plots:
<ul>
<li>Counts and Weights/Strength: Displays the frequency and intensity of interactions among cell groups <strong>(Fig. 7.1a,b)</strong>.</li>
<li>Interaction Heatmap: Shows interaction strengths across all clusters or cell types <strong>(Fig. 7.1c)</strong>.</li>
<li>Incoming and Outgoing Signaling Patterns: Visualizes communication patterns for target and secreting cells <strong>(Fig. 7.1d,e)</strong>.</li>
</ul>
</li>
<li>Interaction Table: Includes source and target cell types, ligand-receptor pairs, and interaction scores <strong>(Fig. 7.1f)</strong>.</li>
</ul>
</li>
<img src='images/7.1.jpg' width='800' height='700' alt=''/>
</ul>
<h3>7.2. Analyzing Specific Signaling Pathways</h3>
For a more focused analysis, users can select a specific signaling pathway from a drop-down menu, enabling detailed visualization of the chosen pathway.
<ul>
<li><strong>Outputs for Specific Signaling Pathway:</strong>
<ul>
<li>Circle Plot: Visualizes interactions among cell groups by counts <strong>(Fig. 7.2a)</strong>.</li>
<li>Chord Plot: Depicts connections between cell types via ligand-receptor pairs <strong>(Fig. 7.2b)</strong>.</li>
<li>Interaction Heatmap: Interaction strengths among clusters for the specific pathway <strong>(Fig. 7.2c)</strong>.</li>
<li>Bubble Plot and Bar Plot: Display interaction intensity for the selected pathway <strong>(Fig. 7.2d)</strong>.</li>
<li>Hierarchy Plot: Shows the hierarchical organization of cell types and their interactions <strong>(Fig. 7.2e)</strong>.</li>
<li>Bar Plot: Shows the network analysis contribution in bar plot <strong>(Fig. 7.2f)</strong>.</li>	
<li>Violin Plot: Shows expression of pathway-associated genes <strong>(Fig. 7.2g)</strong>.</li>
<li>Signaling Pathway Table: Contains source, target, ligand, receptor, and interaction details for the specific pathway <strong>(Fig. 7.2h)</strong>.</li>
<img src='images/7.2.jpg' width='800' height='650' alt=''/>
</ul>
</li>
</ul>
This suite of tools and visualizations enables detailed exploration of cell communication, allowing users to interpret inter-cellular signaling dynamics in ScRNA-seq datasets with biological relevance.
<hr>
<h3>8. Trajectory and Pseudotime Analysis</h3>	
ScRDAVis integrates Monocle3 for trajectory and pseudotime analysis, allowing users to study the dynamic progression of cells over pseudotime and identify genes with functional changes along this trajectory.
<ul>
<li><strong>Preparing for Trajectory and Pseudotime Analysis</strong>
<ul>
<li>Prerequisites: Users must complete analysis up to the cell type prediction step in either single or multiple sample analysis, or subclustering analysis.</li>
<li>Input Format: The tool automatically converts the Seurat object to Monocle3 format, and users can choose between Seurat clusters or predicted cell type labels as input.</li>
<li>UMAP Requirement: UMAP should be used in clustering steps for compatibility with Monocle3.</li>
</ul>
</li>
</ul>
<h3>8.1. Parameters for Learning Trajectory</h3>
<ul>
<li><strong>Partitioning Options:</strong>
<ul>
<li>use_partition: Toggle to specify partitions for different groups.</li>
<li>close_loop: Set to close or open the trajectory loop.</li>
<li>label_groups_by_cluster: Labels cell groups by cluster.</li>
<li>label_branch_points, label_roots, label_leaves: Allows labeling of key points on the trajectory (branches, roots, leaves).</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
<li>Once parameters are set, users can click the Learn Trajectory button to generate the trajectory plot.</li>
</ul>
</li>
<li><strong>Output:</strong>
<ul>
<li>Trajectory Plot: Displays cell progression in trajectory space, providing insight into the cellular development path <strong>(Fig. 8a)</strong>.</li>
</ul>
</li>
</ul>
<h3>8.2. Pseudotime Ordering of Cells</h3>
<ul>
<li><strong>Parameters: </strong>
  <ul>
    <li>Root Cluster Selection: Users must select one cluster to serve as the root cluster, marking the starting point of pseudotime.</li>
    <li>Labeling Options: Parameters include options to label groups by clusters, as well as marking branch points, roots, and leaves.</li>
  </ul>
</li>
	<li><strong>Execution:</strong>
<ul>
<li>Click to Submit button to start the analysis.</li>
</ul>
</li>
<li><strong>Output:</strong>
<ul>
<li>Pseudotime Plot: Cells are arranged by pseudotime, showing the developmental trajectory <strong>(Fig. 8b)</strong>.</li>
<li>Bar Chart: Cells are ordered based on both Seurat clusters and Monocle3 pseudotime <strong>(Fig. 8c)</strong>.</li>
</ul>
</li>
</ul>
<h3>8.3. Identifying Genes with Functional Changes in Pseudotime</h3>
To explore gene expression dynamics along the pseudotime trajectory, users can analyze gene expression changes:
<ul>
<li><strong>Parameters:</strong>
<ul>
<li>Neighbor Graph Selection: Users can select between Principal Graph or K-Nearest Neighbor (KNN) to model gene expression changes.</li>
</ul>
</li>
	<li><strong>Execution:</strong>
<ul>
<li>Click Find Genes Button: Begins the identification of genes whose functions vary along pseudotime.</li>
</ul>
</li>
<li><strong>Output:</strong>
<ul>
<li>Pseudotime Plot of Cells: Visual representation of cells in pseudotime with associated gene expression <strong>(Fig. 8d)</strong>.</li>
<li>Summary Table: Lists genes with dynamic functional changes along pseudotime <strong>(Fig. 8e)</strong>.</li>
</ul>
</li>
</ul>
<h3>8.4. Plotting Gene Expression in Pseudotime</h3>
Users can visualize specific genes to observe their expression patterns over pseudotime:
<ul>
<li><strong>Gene Selection: </strong>
  <ul>
    <li>Top Genes: By default, the tool plots the top 5 genes with dynamic changes, adjustable between 1 to 10 genes.</li>
    <li>Custom Genes: Users can specify a custom list of genes (comma-separated) to plot in pseudotime.</li>
	 </ul>
</li>
    <li><strong>Output: </strong>
      <ul>
        <li>Creates a feature plot to display gene expression across cells in pseudotime <strong>(Fig. 8f)</strong>.</li>
		 <img src='images/8.1.jpg' width='800' height='650' alt=''/>
      </ul>
  </li>
</ul>
This functionality helps users analyze and visualize gene dynamics, offering insights into cellular progression and identifying key genes in developmental pathways.
<hr>
<h3> 9. Co-Expression and TF Analysis</h3>
<h3>9.1. Co-Expression Network Analysis</h3>
ScRDAVis incorporates co-expression network analysis for ScRNA-seq data using the hdWGCNA package. This feature enables users to identify gene modules and their relationships in Seurat clusters or predicted cell type labels. 
<ul>
<li><strong>Prerequisites:</strong>
<ul>
<li>Co-expression network analysis becomes available after completing single or multiple samples analysis or subclustering analysis up to cell type prediction.</li>
<li>User can use one cluster at a time.</li>
</ul>
</li>
<li><strong>Metacell Construction:</strong><br>
Aggregates small groups of similar cells from the same biological sample. Uses the k-Nearest Neighbors (KNN) algorithm to group similar cells and compute a metacell gene expression matrix.
<ul>
<li><strong>Parameters:</strong>
<ul>
<li>k: Number of nearest neighbors for aggregation.</li>
<li>min_cells: Minimum number of cells in a group to construct metacells.</li>
<li>max_shared: Maximum number of cells shared across two metacells.</li>
<li>target_metacells: Maximum number of target metacells to construct.</li>
</ul>
</li>
</ul>
</li>
<li><strong>Co-Expression Network Construction:</strong><br>
Builds networks with customizable parameters:
<ul>
<li>softpower: Determines the scale-free topology for constructing networks.</li>
<li>networkType: Options include signed, unsigned, or signed hybrid.</li>
</ul>
</li>
<li><strong>Module Eigengenes and Connectivity:</strong>
<ul>
<li>Scales data using selectable models: linear, poisson, or negbinom.</li>
<li>Allows Harmony batch correction for harmonized module eigengenes (hMEs), selectable by the users.</li>
</ul>
</li>
<li><strong>Hub Gene Extraction:</strong>
<ul>
<li>Extracts the top N hub genes for selected modules, aiding in the identification of key regulators.</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
	<li>Click the WGCNA Analysis button initiates co-expression network analysis.</li>
</ul>
</li>
<li><strong>Outputs:</strong><br>
Few plots were not available in image files format so we have provided those as pdf files.
<ul>
<li>Soft Power Plots: Visualizes the selection of the optimal soft power parameter for network construction <strong>(Fig. 9.1a)</strong>.</li>
<li>Co-Expression Network Visualization: Displays modules with distinct colors representing gene clusters <strong>(Fig. 9.1b)</strong>.</li>
<li>Ranked Genes in Modules: Provides a list of genes ranked by module membership (kME) <strong>(Fig. 9.1c)</strong>.</li>
<li>Feature Plots: Highlights the expression of modules or specific genes <strong>(Fig. 9.1d)</strong>.</li>
<li>Module Relationships Plots: Correlation between modules based on harmonized module eigengenes (hMEs) <strong>(Fig. 9.1e)</strong>.</li>
<li>Seurat DotPlot with Modules: Displays module-specific gene expression across clusters <strong>(Fig. 9.1f)</strong>.</li>
<li>Individual Module Network Plots: Visualizes the gene network for specific modules <strong>(Fig. 9.1g)</strong>.</li>
<li>Module UMAP Plots: Maps modules onto UMAP visualizations for spatial context <strong>(Fig. 9.1h)</strong>.</li>
<li>Summary Table: Soft Power Table: Lists optimal soft power values <strong>(Fig. 9.1i)</strong>. Module Assignment Table: Details gene-module relationships with colors <strong>(Fig. 9.1j)</strong>. Hub Genes Table: Identifies top hub genes per module <strong>(Fig. 9.1k)</strong>.</li>
<img src='images/9.1.jpg' width='800' height='800' alt=''/>
</ul>
</li>
</ul>
This functionality provides a robust framework for uncovering intricate co-expression patterns and identifying key drivers in single-cell datasets.<br>
<h3>9.2. Transcription Factor Regulatory Network Analysis</h3>
Transcription Factor (TF) Regulatory Network Analysis in ScRDAVis employs the hdWGCNA package to construct and analyze  TF regulatory networks based on ScRNA-seq data. This feature allows users to identify gene modules and investigate TF-mediated regulation within clusters or predicted cell type labels.
<ul>
<li><strong>Prerequisites:</strong>
<ul>
<li>Complete single or multiple sample analysis or subclustering analysis, including cell type prediction.</li>
<li>Analysis is performed one cluster at a time.</li>
</ul>
</li>
<li><strong>TF Regulatory Network Construction:</strong>
<ul>
<li><strong>TF Binding Motif Information:</strong>
<ul>
<li>Human: EnsDb.Hsapiens.v86, BSgenome.Hsapiens.UCSC.hg38.</li>
<li>Mouse: EnsDb.Mmusculus.v79, BSgenome.Mmusculus.UCSC.mm10.</li>
<li>Motifs from the JASPAR 2020 database for multiple species.</li>
</ul>
</li>
<li><strong>Machine Learning Model:</strong>
<ul>
<li>XGBoost used to model TF regulation for each gene with:</li>
<li>max_depth : Maximum depth of a tree</li>
<li>eta : Step size shrinkage used in update to prevent overfitting</li>
<li>alpha: L1 regularization term on weights</li>
</ul>
</li>
<li><strong>TF Regulon Strategy:</strong>
<ul>
<li>Strategy A selects the top TFs for each gene by default</li>
<li>reg_thresh : Threshold for regulatory score)</li>
<li>n_tfs : The number of top TFs to keep for each gene</li>
</ul>
</li>
<li><strong>Regulon Expression Signatures:</strong>
<ul>
<li>Positive correlation: cor_thresh = 0.05. Threshold for TF-gene correlation for genes to be included in the positive regulon score</li>
<li>Negative correlation: cor_thresh = -0.05. threshold for TF-gene correlation for genes to be included in the negative regulon score</li>
</ul>
</li>
</ul>
</li>
<li><strong>Execution:</strong>
<ul>
<li>Click Transcription factor analysis button to start the analysis.</li>
</ul>
</li>	
<li><strong>Output and Visualization:</strong>
<ul>
<li>Module Regulatory Network Plots: Positive, negative, and combined regulatory network plots. Visualize TF-to-target relationships categorized by regulatory effects <strong>(Fig. 9.2.1a-d)</strong>.</li>
<li>Regulated Scores Table: Comprehensive list of TFs and their downstream targets <strong>(Fig. 9.2.1e)</strong>.</li>
	<img src='images/9.2.1.jpg' width='800' height='800' alt=''/>
</ul>
</li>
<li><strong>TF-Specific Visualizations:</strong><br>
Unravel regulatory mechanisms governing gene expression in cellular contexts. Identify key transcription factors and their target genes for hypothesis generation and validation. Explore positive and negative regulatory effects within gene modules.
<ul>
<li>Select a TF from a dropdown menu to generate specific plots:</li>
</ul>
</li>
<li><strong>Outputs:</strong>
<ul>
<li>UMAP Plots: Spatial distribution of the TF <strong>(Fig. 9.2.2a)</strong>.</li>
<li>Bar Plots: Contribution of the TF across modules <strong>(Fig. 9.2.2b)</strong>.</li>
<li>Network Plots: Positive, negative, and combined networks, with primary, secondary and tertiary targets <strong>(Fig. 9.2.2c-e)</strong>.</li>
<img src='images/9.2.2.jpg' width='800' height='800' alt=''/>
</ul>
</li>
</ul>
This functionality provides a comprehensive view of transcriptional regulation in ScRNA-seq data, enabling detailed exploration of TF-driven cellular processes.

	
<hr>

         ") 
  ),
),





)
)

