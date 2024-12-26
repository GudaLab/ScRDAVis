source("global.R")
server <- function(input, output, session) {
  
 
############################################################################################################################################################
                                                                  ##    Multiple Input     ##
############################################################################################################################################################
  
  ##########multiple sidebar hide##########
  observeEvent(input[["multiple_tabsets"]], {
    if(input[["multiple_tabsets"]] == "Stats"){
      showElement(selector = "#multiple_sidebar")
      removeCssClass("multiple_main_menu", "col-sm-12")
      addCssClass("multiple_main_menu", "col-sm-8")
    }else{
      hideElement(selector = "#multiple_sidebar")
      removeCssClass("multiple_main_menu", "col-sm-8")
      addCssClass("multiple_main_menu", "col-sm-12")
    }
  })
  
  
  ################multiple hide tabs########
  hideTab(inputId = "multiple_tabsets", target = "Sample Groups and QC Filtering")
  hideTab(inputId = "multiple_tabsets", target = "Normalization and PCA analysis")
  hideTab(inputId = "multiple_tabsets", target = "Clustering")
  hideTab(inputId = "multiple_tabsets", target = "Remove Doublets")
  hideTab(inputId = "multiple_tabsets", target = "Markers Identification")
  hideTab(inputId = "multiple_tabsets", target = "Cell Type Prediction")
  hideTab(inputId = "multiple_tabsets", target = "Cluster-based plots")
  hideTab(inputId = "multiple_tabsets", target = "Condition based analysis")
  
  ################showtabbutton#############     
  observeEvent(input$link_m_qc_filtering, {
    showTab(inputId = "multiple_tabsets", target = "Sample Groups and QC Filtering")
  })
  
  observeEvent(input$link_m_normalization, {
    showTab(inputId = "multiple_tabsets", target = "Normalization and PCA analysis")
  })
  
  observeEvent(input$link_m_clustering, {
    showTab(inputId = "multiple_tabsets", target = "Clustering")
  })
  observeEvent(input$link_m_doublet, {
    showTab(inputId = "multiple_tabsets", target = "Remove Doublets")
  })
  
  observeEvent(input$link_m_marker, {
    showTab(inputId = "multiple_tabsets", target = "Markers Identification")
  })
  
  observeEvent(input$link_m_prediction, {
    showTab(inputId = "multiple_tabsets", target = "Cell Type Prediction")
  })
  
  observeEvent(input$link_m_clusterbased, {
    showTab(inputId = "multiple_tabsets", target = "Cluster-based plots")
  })
  
  observeEvent(input$link_m_conditionbased, {
    showTab(inputId = "multiple_tabsets", target = "Condition based analysis")
  })
  
  
  ################ multiple input###########   
  observe({
    if (input$multiple_sample_format == "h5") {
      shinyjs::show("multiple_sample_file")
      shinyjs::hide("multiple_sample_file_mfb")
      shinyjs::hide("multiple_sample_file_rds")
      shinyjs::hide("multiple_sample_file_txt")
    }
    else if (input$multiple_sample_format == "MFB") {
      shinyjs::hide("multiple_sample_file")
      shinyjs::show("multiple_sample_file_mfb")
      shinyjs::hide("multiple_sample_file_rds")
      shinyjs::hide("multiple_sample_file_txt")
    }
    else if (input$multiple_sample_format == "seurat_object") {
      shinyjs::hide("multiple_sample_file")
      shinyjs::hide("multiple_sample_file_mfb")
      shinyjs::show("multiple_sample_file_rds")
      shinyjs::hide("multiple_sample_file_txt")
    }
    else if (input$multiple_sample_format == "matrix_count") {
      shinyjs::hide("multiple_sample_file")
      shinyjs::hide("multiple_sample_file_mfb")
      shinyjs::hide("multiple_sample_file_rds")
      shinyjs::show("multiple_sample_file_txt")
    }
    
  })
  
  
  ########multiple hide qc filtering######## 
  shinyjs::hide("m_bf_box1")
  shinyjs::hide("m_bf_box2")
  shinyjs::hide("m_bf_box3")
  shinyjs::hide("m_bf_box4")
  shinyjs::hide("m_bf_box5")
  
  observeEvent(input$multiple_sample_submit,{
    shinyjs::show("m_bf_box1")
    shinyjs::show("m_bf_box2")
    shinyjs::show("m_bf_box3")
    shinyjs::show("m_bf_box4")
    shinyjs::show("m_bf_box5")
  })
  
  
  ########multiple hide qc filtering########  
  shinyjs::hide("m_qc_filter_box1")
  shinyjs::hide("m_qc_filter_box2")
  shinyjs::hide("m_qc_filter_box3")
  shinyjs::hide("m_qc_filter_box4")
  shinyjs::hide("m_qc_filter_box5")
  shinyjs::hide("m_qc_filter_box6")
  
  observeEvent(input$multiple_sample_qc_filtering,{
    shinyjs::show("m_qc_filter_box1")
    shinyjs::show("m_qc_filter_box2")
    shinyjs::show("m_qc_filter_box3")
    shinyjs::show("m_qc_filter_box4")
    shinyjs::show("m_qc_filter_box5")
    shinyjs::show("m_qc_filter_box6")
  })
  
  observe({
    if (input$multiple_group_count == 1) {
      shinyjs::show("group1_name")
      shinyjs::show("group1_samples")
      shinyjs::hide("group2_name")
      shinyjs::hide("group2_samples")
      shinyjs::hide("group3_name")
      shinyjs::hide("group3_samples")
      shinyjs::hide("group4_name")
      shinyjs::hide("group4_samples")
      shinyjs::hide("group5_name")
      shinyjs::hide("group5_samples")
      shinyjs::hide("group6_name")
      shinyjs::hide("group6_samples")
      shinyjs::hide("group2_samples2")
      shinyjs::hide("group3_samples3")
      shinyjs::hide("group4_samples4")
      shinyjs::hide("group5_samples5")
      shinyjs::hide("group6_samples6")
      hideTab(inputId = "multiple_tabsets", target = "Condition based analysis")
      hideTab(inputId = "subclustering_multiple_tabsets", target = "Condition based analysis")
      shinyjs::hide("link_m_conditionbased")
      shinyjs::hide("link_m_subclustering_conditionbased")
    }
   else if (input$multiple_group_count == 2) {
      shinyjs::show("group1_name")
      shinyjs::show("group1_samples")
      shinyjs::show("group2_name")
      shinyjs::show("group2_samples")
      shinyjs::hide("group3_name")
      shinyjs::hide("group3_samples")
      shinyjs::hide("group4_name")
      shinyjs::hide("group4_samples")
      shinyjs::hide("group5_name")
      shinyjs::hide("group5_samples")
      shinyjs::hide("group6_name")
      shinyjs::hide("group6_samples")
      shinyjs::show("group2_samples2")
      shinyjs::hide("group3_samples3")
      shinyjs::hide("group4_samples4")
      shinyjs::hide("group5_samples5")
      shinyjs::hide("group6_samples6")
      #showTab(inputId = "multiple_tabsets", target = "Condition based analysis")
      #showTab(inputId = "subclustering_multiple_tabsets", target = "Condition based analysis")
      shinyjs::show("link_m_conditionbased")
      shinyjs::show("link_m_subclustering_conditionbased")
    }
    else if (input$multiple_group_count == 3) {
      shinyjs::show("group1_name")
      shinyjs::show("group1_samples")
      shinyjs::show("group2_name")
      shinyjs::show("group2_samples")
      shinyjs::show("group3_name")
      shinyjs::show("group3_samples")
      shinyjs::hide("group4_name")
      shinyjs::hide("group4_samples")
      shinyjs::hide("group5_name")
      shinyjs::hide("group5_samples")
      shinyjs::hide("group6_name")
      shinyjs::hide("group6_samples")
      shinyjs::show("group2_samples2")
      shinyjs::show("group3_samples3")
      shinyjs::hide("group4_samples4")
      shinyjs::hide("group5_samples5")
      shinyjs::hide("group6_samples6")
      shinyjs::show("link_m_conditionbased")
      shinyjs::show("link_m_subclustering_conditionbased")
    }
    else if (input$multiple_group_count == 4) {
      shinyjs::show("group1_name")
      shinyjs::show("group1_samples")
      shinyjs::show("group2_name")
      shinyjs::show("group2_samples")
      shinyjs::show("group3_name")
      shinyjs::show("group3_samples")
      shinyjs::show("group4_name")
      shinyjs::show("group4_samples")
      shinyjs::hide("group5_name")
      shinyjs::hide("group5_samples")
      shinyjs::hide("group6_name")
      shinyjs::hide("group6_samples")
      shinyjs::show("group2_samples2")
      shinyjs::show("group3_samples3")
      shinyjs::show("group4_samples4")
      shinyjs::hide("group5_samples5")
      shinyjs::hide("group6_samples6")
      shinyjs::show("link_m_conditionbased")
      shinyjs::show("link_m_subclustering_conditionbased")
    }
    else if (input$multiple_group_count == 5) {
      shinyjs::show("group1_name")
      shinyjs::show("group1_samples")
      shinyjs::show("group2_name")
      shinyjs::show("group2_samples")
      shinyjs::show("group3_name")
      shinyjs::show("group3_samples")
      shinyjs::show("group4_name")
      shinyjs::show("group4_samples")
      shinyjs::show("group5_name")
      shinyjs::show("group5_samples")
      shinyjs::hide("group6_name")
      shinyjs::hide("group6_samples")
      shinyjs::show("group2_samples2")
      shinyjs::show("group3_samples3")
      shinyjs::show("group4_samples4")
      shinyjs::show("group5_samples5")
      shinyjs::hide("group6_samples6")
      shinyjs::show("link_m_conditionbased")
      shinyjs::show("link_m_subclustering_conditionbased")
    }
    else if (input$multiple_group_count == 6) {
      shinyjs::show("group1_name")
      shinyjs::show("group1_samples")
      shinyjs::show("group2_name")
      shinyjs::show("group2_samples")
      shinyjs::show("group3_name")
      shinyjs::show("group3_samples")
      shinyjs::show("group4_name")
      shinyjs::show("group4_samples")
      shinyjs::show("group5_name")
      shinyjs::show("group5_samples")
      shinyjs::show("group6_name")
      shinyjs::show("group6_samples")
      shinyjs::show("group2_samples2")
      shinyjs::show("group3_samples3")
      shinyjs::show("group4_samples4")
      shinyjs::show("group5_samples5")
      shinyjs::show("group6_samples6")
      shinyjs::show("link_m_conditionbased")
      shinyjs::show("link_m_subclustering_conditionbased")
    }
  })
  
  
  ########multiple hide normalization########     
  shinyjs::hide("m_pca_box1")
  shinyjs::hide("m_elbow_box")
  shinyjs::hide("m_pca_box2")
  shinyjs::hide("m_pca_box3")
  shinyjs::hide("m_pca_box4")
  
  observeEvent(input$multiple_sample_normalization,{
    shinyjs::show("m_pca_box1")
    shinyjs::show("m_elbow_box")
    shinyjs::show("m_pca_box2")
    shinyjs::show("m_pca_box3")
    shinyjs::show("m_pca_box4")
  })
  
  observe({
    if (input$multiple_sample_normalization_method == "LogNormalize") {
      shinyjs::show("multiple_sample_scale_factor")
      shinyjs::show("multiple_sample_normalization_variable_genes")
      shinyjs::show("multiple_sample_var_genes")
    }
    else if (input$multiple_sample_normalization_method  == "SCTransform") {
      shinyjs::hide("multiple_sample_scale_factor")
      shinyjs::hide("multiple_sample_normalization_variable_genes")
      shinyjs::hide("multiple_sample_var_genes")
    }
  })
  
  ########multiple hide clustering########     
  shinyjs::hide("m_clustering_box1")
  shinyjs::hide("m_clustering_box2")
  shinyjs::hide("m_clustering_box3")
  shinyjs::hide("m_clustering_box4")
  shinyjs::hide("m_clustering_box5")
  shinyjs::hide("m_clustering_box6") 
  shinyjs::hide("m_clustering_box7") 
  shinyjs::hide("m_clustering_box8")   
  shinyjs::hide("m_clustering_box9") 
  shinyjs::hide("m_clustering_box10") 
  
  observeEvent(input$multiple_sample_clustering,{
    shinyjs::show("m_clustering_box1")
    shinyjs::show("m_clustering_box2")
    shinyjs::show("m_clustering_box3")
    shinyjs::show("m_clustering_box4")
    shinyjs::show("m_clustering_box5")
    shinyjs::show("m_clustering_box6") 
    shinyjs::show("m_clustering_box7") 
    shinyjs::show("m_clustering_box8") 
    shinyjs::show("m_clustering_box9") 
    shinyjs::show("m_clustering_box10") 
  })
  
  
  observe({
    if (input$m_clustering6 == "umap") {
      shinyjs::show("m_umap_box")
      shinyjs::hide("m_tsne_box")
    }
    else if (input$m_clustering6  == "tsne") {
      shinyjs::hide("m_umap_box")
      shinyjs::show("m_tsne_box")
    }
  })
  
  
  ########multiple hide doublet boxes########
  shinyjs::hide("m_doublet_box2")
  shinyjs::hide("m_doublet_box3")
  shinyjs::hide("m_doublet_box4")
  shinyjs::hide("m_doublet_box5")
  shinyjs::hide("m_doublet_box6")
  shinyjs::hide("m_doublet_box7")
  shinyjs::hide("m_doublet_box8")
  shinyjs::hide("m_doublet_box9")
  shinyjs::hide("m_doublet_box10")
  shinyjs::hide("m_doublet_box11")
  shinyjs::hide("m_doublet_box12")
  shinyjs::hide("m_doublet_box13")
  shinyjs::hide("m_doublet_box14")
  shinyjs::hide("m_doublet_box15")
  shinyjs::hide("m_doublet_box16")
  shinyjs::hide("m_doublet_box17")
  shinyjs::hide("m_doublet_box18")
  
  observeEvent(input$multiple_sample_doublet,{
    shinyjs::show("m_doublet_box2")
    shinyjs::show("m_doublet_box3")
    shinyjs::show("m_doublet_box4")
    shinyjs::show("m_doublet_box5")
    shinyjs::show("m_doublet_box6")
    shinyjs::show("m_doublet_box7")
  })
  observeEvent(input$multiple_sample_doublet2,{
    
    shinyjs::show("m_doublet_box8")
    shinyjs::show("m_doublet_box9")
    shinyjs::show("m_doublet_box10")
    shinyjs::show("m_doublet_box11")
    shinyjs::show("m_doublet_box12")
    shinyjs::show("m_doublet_box13")
    shinyjs::show("m_doublet_box14")
    shinyjs::show("m_doublet_box15")
    shinyjs::show("m_doublet_box16")
    shinyjs::show("m_doublet_box17")
    shinyjs::show("m_doublet_box18")
    
  })
  
  
  ########multiple hide markers box########  
  shinyjs::hide("m_marker_box5")
  shinyjs::hide("m_marker_box6")
  shinyjs::hide("m_marker_box7")
  shinyjs::hide("m_marker10")
  shinyjs::hide("m_marker11")
  shinyjs::hide("m_marker12")
  
  
  observe({
    if (input$m_marker1 == 1) {
      shinyjs::hide("m_marker_6")
      shinyjs::hide("m_marker_7")
      shinyjs::hide("m_marker_8")
      shinyjs::hide("m_marker_9")
      shinyjs::hide("m_marker6")
      shinyjs::hide("m_marker7")
      shinyjs::hide("m_marker8")
      shinyjs::hide("m_marker9")
      shinyjs::hide("m_marker10")
      shinyjs::show("m_marker11")
      shinyjs::hide("m_marker12")
    }
    else if (input$m_marker1 == 2) {
      shinyjs::show("m_marker_6")
      shinyjs::hide("m_marker_7")
      shinyjs::hide("m_marker_8")
      shinyjs::hide("m_marker_9")
      shinyjs::show("m_marker6")
      shinyjs::hide("m_marker7")
      shinyjs::hide("m_marker8")
      shinyjs::hide("m_marker9")
      shinyjs::hide("m_marker_box6")
      shinyjs::hide("m_marker10")
      shinyjs::show("m_marker11")
      shinyjs::hide("m_marker12")
      
    }
    else if (input$m_marker1 == 3) {
      shinyjs::show("m_marker_6")
      shinyjs::show("m_marker_7")
      shinyjs::hide("m_marker_8")
      shinyjs::hide("m_marker_9")
      shinyjs::show("m_marker6")
      shinyjs::show("m_marker7")
      shinyjs::hide("m_marker8")
      shinyjs::hide("m_marker9")
      shinyjs::hide("m_marker_box6")
      shinyjs::hide("m_marker10")
      shinyjs::show("m_marker11")
      shinyjs::hide("m_marker12")
    }
    else if (input$m_marker1 == 4) {
      shinyjs::hide("m_marker_6")
      shinyjs::hide("m_marker_7")
      shinyjs::show("m_marker_8")
      shinyjs::hide("m_marker_9")
      shinyjs::hide("m_marker6")
      shinyjs::hide("m_marker7")
      shinyjs::show("m_marker8")
      shinyjs::hide("m_marker9")
      shinyjs::hide("m_marker_box6")
      shinyjs::show("m_marker10")
      shinyjs::hide("m_marker11")
      shinyjs::show("m_marker12")
    }
    
    else if (input$m_marker1 == 5) {
      shinyjs::hide("m_marker_6")
      shinyjs::hide("m_marker_7")
      shinyjs::show("m_marker_8")
      shinyjs::show("m_marker_9")
      shinyjs::hide("m_marker6")
      shinyjs::hide("m_marker7")
      shinyjs::show("m_marker8")
      shinyjs::show("m_marker9")
      shinyjs::hide("m_marker_box6")
      shinyjs::show("m_marker10")
      shinyjs::hide("m_marker11")
      shinyjs::show("m_marker12")
    }
    
  })
  
  ########multiple hide celltype box########  
  shinyjs::hide("m_celltype_box3")
  shinyjs::hide("m_celltype_box4")
  shinyjs::hide("m_celltype_box5")
  shinyjs::hide("m_celltype_box7")
  shinyjs::hide("m_celltype_box8")
  shinyjs::hide("m_celltype_box9")
  shinyjs::hide("m_celltype_box10")
  
  observe({
    if (input$m_celltype1 == 1) {
      shinyjs::show("m_celltype_box2")
      shinyjs::show("m_celltype2")
      shinyjs::hide("m_celltype_box3")
      shinyjs::hide("m_celltype3")
      shinyjs::hide("m_celltype4")
      shinyjs::hide("m_celltype_box4")
      shinyjs::hide("m_celltype5")
      shinyjs::hide("m_celltype6")
      shinyjs::hide("m_celltype_box5")
      shinyjs::hide("m_celltype7")
    }
    else if (input$m_celltype1 == 2) {
      shinyjs::hide("m_celltype_box2")
      shinyjs::hide("m_celltype2")
      shinyjs::show("m_celltype_box3")
      shinyjs::show("m_celltype3")
      shinyjs::show("m_celltype4")
      shinyjs::hide("m_celltype_box4")
      shinyjs::hide("m_celltype5")
      shinyjs::hide("m_celltype6")
      shinyjs::hide("m_celltype_box5")
      shinyjs::hide("m_celltype7")
    }
    else if (input$m_celltype1 == 3) {
      shinyjs::hide("m_celltype_box2")
      shinyjs::hide("m_celltype2")
      shinyjs::hide("m_celltype_box3")
      shinyjs::hide("m_celltype3")
      shinyjs::hide("m_celltype4")
      shinyjs::show("m_celltype_box4")
      shinyjs::show("m_celltype5")
      shinyjs::show("m_celltype6")
      shinyjs::hide("m_celltype_box5")
      shinyjs::hide("m_celltype7")
    }
    else if (input$m_celltype1 == 4) {
      shinyjs::hide("m_celltype_box2")
      shinyjs::hide("m_celltype2")
      shinyjs::hide("m_celltype_box3")
      shinyjs::hide("m_celltype3")
      shinyjs::hide("m_celltype4")
      shinyjs::hide("m_celltype_box4")
      shinyjs::hide("m_celltype5")
      shinyjs::hide("m_celltype6")
      shinyjs::show("m_celltype_box5")
      shinyjs::show("m_celltype7")
    }
  })     
  
  ##################multiple Cluster-based plots####################### 
  shinyjs::hide("m_clusterbased2")
  shinyjs::hide("m_clusterbased_box2")
  shinyjs::hide("m_clusterbased_box3")     
  shinyjs::hide("m_clusterbased_box4")
  
  ##################multiple conditionbased####################### 
  shinyjs::hide("m_conditionbased_box3")
  shinyjs::hide("m_conditionbased_box4")     
  shinyjs::hide("m_conditionbased_box5")      
  
  #####################Tab1##############################
  ######################data Input##################
  datainput_multiple_sample_level<- eventReactive(input$multiple_sample_submit,{
    
    if (input$multiple_sample_format == "h5") {
      upload_multiple_sample_file <- input$multiple_sample_file
      upload_multiple_sample_file_names <- input$multiple_sample_file$name
      type=input$multiple_sample_format
      ext <- tools::file_ext(input$multiple_sample_file)
      req(input$multiple_sample_file)
      }
    else if (input$multiple_sample_format == "MFB") {
      
      file1 <- input$multiple_sample_file_mfb[['datapath']]
      filesdir = dirname(file1)
      
      file.rename(file1, paste0(filesdir,'/',input$multiple_sample_file_mfb$name))
      upload_multiple_sample_file <- filesdir
      upload_multiple_sample_file_names <- input$multiple_sample_file_mfb$name
    }
    else if (input$multiple_sample_format == "seurat_object") {
      upload_multiple_sample_file <- input$multiple_sample_file_rds
      upload_multiple_sample_file_names <- input$multiple_sample_file_rds$name
      type=input$multiple_sample_format
      ext <- tools::file_ext(input$multiple_sample_file_rds)
      req(input$multiple_sample_file_rds)
      }
    else if (input$multiple_sample_format == "matrix_count") {
      upload_multiple_sample_file <- input$multiple_sample_file_txt
      upload_multiple_sample_file_names <- input$multiple_sample_file_txt$name
      type=input$multiple_sample_format
      ext <- tools::file_ext(input$multiple_sample_file_txt)
      req(input$multiple_sample_file_txt)
      }
    
    
    
    source("scripts/multiple_file_upload.R")
    datainput_multiple_sample(index_multiple_sample_file = upload_multiple_sample_file$datapath, index_multiple_sample_file_names = upload_multiple_sample_file_names, index_multiple_sample_file1 = filesdir, index_multiple_sample_format=input$multiple_sample_format, index_multiple_sample_name = input$multiple_sample_name, index_multiple_sample_cell = input$multiple_sample_cell, index_multiple_sample_genes = input$multiple_sample_genes)
    
  })
  
  
  output$m_qc_before_filtering <- renderPlot({
    datainput_multiple_sample_level()[1]
  })
  ################m_QCplot############################
  observeEvent(input$download_m_qc_before_filtering, {
    showModal(modalDialog(
      title = strong("Download QC plot"),
      numericInput("m_qc_before_filtering_plot_height", label = h5("Figure height (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_qc_before_filtering_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_qc_before_filtering_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_qc_before_filtering_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_qc_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  
  output$m_qc_downloadoutput<- downloadHandler(
    filename = function(){
      paste("QC_before_filtering", input$multiple_sample_name, input$m_qc_before_filtering_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_sample_level()[[1]], width = input$m_qc_before_filtering_plot_width, height = input$m_qc_before_filtering_plot_height, dpi = input$m_qc_before_filtering_plot_dpi, units = "in")
    }
  )
  
  
  output$multiple_cell_table<- renderDataTable(DT::datatable((datainput_multiple_sample_level()[[2]]),
                                                             options = list(
                                                               scrollX = TRUE,
                                                               pageLength = 10,
                                                               bFilter=0
                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_multiple_cell_table <- downloadHandler(
    filename = function() { 
      paste("Number of cells", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_sample_level()[[2]], file)
    }
  )
  
  
  output$m_ff_before_filtering <- renderPlot({
    datainput_multiple_sample_level()[3]
  })
  
  #######################m_ff_plot#########################
  observeEvent(input$download_m_ff_before_filtering, {
    showModal(modalDialog(
      title = strong("Download feature-feature plot"),
      numericInput("m_ff_before_filtering_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_ff_before_filtering_plot_width", label = h5("Figure width (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_ff_before_filtering_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_ff_before_filtering_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_ff_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  
  output$m_ff_downloadoutput<- downloadHandler(
    filename = function(){
      paste("feature_feature_relationships_plot", input$multiple_sample_name, input$m_ff_before_filtering_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_sample_level()[[3]], width = input$m_ff_before_filtering_plot_width, height = input$m_ff_before_filtering_plot_height, dpi = input$m_ff_before_filtering_plot_dpi, units = "in")
    }
  )
  
  
  output$text_level<- renderText({
    paste(datainput_multiple_sample_level()[[4]])
  })
  
  
  
  output$m_so_before_filtering<- downloadHandler(
    filename = function(){
      paste("seuart_object_before_qc.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_sample_level()[[5]], file= file)
    }
  )
  
  ###############link to next tab###########################      
  observeEvent(input$link_m_qc_filtering, {
    newvalue <- "Sample Groups and QC Filtering"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  }) 
  
  ##########################Tab1.2###############################          
  ####################groupnames##############################  
  output$group1_samples1 <- renderUI ({
    
    samples <- req(datainput_multiple_sample_level()[[6]])
    #samples <- samples[!samples == input$group2_samples]
    
    shinyWidgets::pickerInput(
      inputId = "group1_samples",
      label = "Select group1 sample(s)",
      choices = sort(samples),
      multiple = T,
      selected = sort(samples)[1],
      #options = list(`actions-box` = TRUE)
    )
    
  })
  
  output$group2_samples2 <- renderUI ({
    samples <- req(datainput_multiple_sample_level()[[6]])
    remaining_samples <- samples[!samples %in% input$group1_samples]
    
    shinyWidgets::pickerInput(
      inputId = "group2_samples",
      label = "Select group2 sample(s)",
      choices = sort(remaining_samples),
      multiple = T,
      selected = if (length(remaining_samples) > 0) sort(remaining_samples)[1] else NULL,
      #options = list(`actions-box` = TRUE)
    )
    
  })
  
  output$group3_samples3 <- renderUI ({
    samples <- req(datainput_multiple_sample_level()[[6]])
    remaining_samples <- samples[!samples %in% c(input$group1_samples, input$group2_samples)]
    
    shinyWidgets::pickerInput(
      inputId = "group3_samples",
      label = "Select group3 sample(s)",
      choices = sort(remaining_samples),
      multiple = T,
      selected = if (length(remaining_samples) > 0) sort(remaining_samples)[1] else NULL,
      #options = list(`actions-box` = TRUE)
    )
    
  })
  
  output$group4_samples4 <- renderUI ({
    samples <- req(datainput_multiple_sample_level()[[6]])
    remaining_samples <- samples[!samples %in% c(input$group1_samples, input$group2_samples, input$group3_samples)]
    
    shinyWidgets::pickerInput(
      inputId = "group4_samples",
      label = "Select group4 sample(s)",
      choices = sort(remaining_samples),
      multiple = T,
      selected = if (length(remaining_samples) > 0) sort(remaining_samples)[1] else NULL,
      #options = list(`actions-box` = TRUE)
      )
    
  })
  
  output$group5_samples5 <- renderUI ({
    samples <- req(datainput_multiple_sample_level()[[6]])
    remaining_samples <- samples[!samples %in% c(input$group1_samples, input$group2_samples, input$group3_samples, input$group4_samples)]
    
    shinyWidgets::pickerInput(
      inputId = "group5_samples",
      label = "Select group5 sample(s)",
      choices = sort(remaining_samples),
      multiple = T,
      selected = if (length(remaining_samples) > 0) sort(remaining_samples)[1] else NULL,
      #options = list(`actions-box` = TRUE)
    )
    
  })
  
  
  output$group6_samples6 <- renderUI ({
    samples <- req(datainput_multiple_sample_level()[[6]])
    remaining_samples <- samples[!samples %in% c(input$group1_samples, input$group2_samples, input$group3_samples, input$group4_samples, input$group5_samples)]
    
    shinyWidgets::pickerInput(
      inputId = "group6_samples",
      label = "Select group6 sample(s)",
      choices = sort(remaining_samples),
      multiple = T,
      selected = if (length(remaining_samples) > 0) sort(remaining_samples)[1] else NULL,
      #options = list(`actions-box` = TRUE)
    )
    
  })
  
  
  
  
  ##############multiple QC after filtering###################   
  datainput_multiple_qc_filter_level <- eventReactive(input$multiple_sample_qc_filtering,{
    source("scripts/multiple_qc_filter.R")
    datainput_multiple_qc_filter(index_multiple_qc_input = datainput_multiple_sample_level()[[5]], index_multiple_qc_input1 = datainput_multiple_sample_level()[[7]], index_multiple_group_count = input$multiple_group_count, index_group1_name = input$group1_name, index_group1_samples = input$group1_samples, index_group2_name = input$group2_name, index_group2_samples = input$group2_samples, index_group3_name = input$group3_name, index_group3_samples = input$group3_samples, index_group4_name = input$group4_name, index_group4_samples = input$group4_samples, index_group5_name = input$group5_name, index_group5_samples = input$group5_samples, index_group6_name = input$group6_name, index_group6_samples = input$group6_samples, index_multiple_sample_min_count = input$multiple_sample_min_count, index_multiple_sample_max_count=input$multiple_sample_max_count, index_multiple_sample_max_mito_perc=input$multiple_sample_max_mito_perc)
  })
  
  
  
  output$m_qc_after_filtering<- renderPlot({
    datainput_multiple_qc_filter_level()[1]
  })
  
  
  #################m_QC_after_filtering_plot############
  observeEvent(input$download_m_qc_after_filtering, {
    showModal(modalDialog(
      title = strong("Download samples QC plot"),
      numericInput("m_qc_after_filtering_plot_height", label = h5("Figure height (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_qc_after_filtering_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_qc_after_filtering_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_qc_after_filtering_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_qc_after_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_qc_after_downloadoutput<- downloadHandler(
    filename = function(){
      paste("QC_after_filtering_sample_based", input$m_qc_after_filtering_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_qc_filter_level()[[1]], width = input$m_qc_after_filtering_plot_width, height = input$m_qc_after_filtering_plot_height, dpi = input$m_qc_after_filtering_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_qc_after_filtering2<- renderPlot({
    datainput_multiple_qc_filter_level()[2]
  })
  
  observeEvent(input$download_m_qc_after_filtering2, {
    showModal(modalDialog(
      title = strong("Download group QC plot"),
      numericInput("m_qc_after_filtering2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_qc_after_filtering2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_qc_after_filtering2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_qc_after_filtering2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_qc_after_downloadoutput2", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_qc_after_downloadoutput2<- downloadHandler(
    filename = function(){
      paste("QC_after_filtering_group_based", input$m_qc_after_filtering2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_qc_filter_level()[[2]], width = input$m_qc_after_filtering2_plot_width, height = input$m_qc_after_filtering2_plot_height, dpi = input$m_qc_after_filtering2_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_qc_after_filtering3<- renderPlot({
    datainput_multiple_qc_filter_level()[3]
  })
  
  
  observeEvent(input$download_m_qc_after_filtering3, {
    showModal(modalDialog(
      title = strong("Download samples bar plot"),
      numericInput("m_qc_after_filtering3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 10, width = "300px"),
      numericInput("m_qc_after_filtering3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_qc_after_filtering3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_qc_after_filtering3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_qc_after_downloadoutput3", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_qc_after_downloadoutput3<- downloadHandler(
    filename = function(){
      paste("Bar_plot_sample_based", input$m_qc_after_filtering3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_qc_filter_level()[[3]], width = input$m_qc_after_filtering3_plot_width, height = input$m_qc_after_filtering3_plot_height, dpi = input$m_qc_after_filtering3_plot_dpi, units = "in")
    }
  )
  
  
  output$m_qc_after_filtering4<- renderPlot({
    datainput_multiple_qc_filter_level()[4]
  })
  
  observeEvent(input$download_m_qc_after_filtering4, {
    showModal(modalDialog(
      title = strong("Download groups bar plot"),
      numericInput("m_qc_after_filtering4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_qc_after_filtering4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 6, width = "300px"),
      numericInput("m_qc_after_filtering4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_qc_after_filtering4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_qc_after_downloadoutput4", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_qc_after_downloadoutput4<- downloadHandler(
    filename = function(){
      paste("Bar_plot_group_based", input$m_qc_after_filtering4_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_qc_filter_level()[[4]], width = input$m_qc_after_filtering4_plot_width, height = input$m_qc_after_filtering4_plot_height, dpi = input$m_qc_after_filtering4_plot_dpi, units = "in")
    }
  )
  
  
  output$multiple_cell_table_after_qc<- renderDataTable(DT::datatable((datainput_multiple_qc_filter_level()[[5]]),
                                                                      options = list(
                                                                        scrollX = TRUE,
                                                                        pageLength = 10,
                                                                        bFilter=0
                                                                      ),rownames= FALSE, selection = "none"))
  
  output$download_multiple_cell_table_after_qc <- downloadHandler(
    filename = function() { 
      paste("Number of cells in samples after qc", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_qc_filter_level()[[5]], file)
    }
  )
  
  output$multiple_cell_table_after_qc2<- renderDataTable(DT::datatable((datainput_multiple_qc_filter_level()[[6]]),
                                                                       options = list(
                                                                         scrollX = TRUE,
                                                                         pageLength = 10,
                                                                         bFilter=0
                                                                       ),rownames= FALSE, selection = "none"))
  
  output$download_multiple_cell_table_after_qc2 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in groups after qc", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_qc_filter_level()[[6]], file)
    }
  )      
  
  ###################save seurat object after qc###################
  output$m_so_after_filtering<- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_qc.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_qc_filter_level()[[7]], file= file)
    }
  )
  
  ###############link to next tab###########################      
  observeEvent(input$link_m_normalization, {
    newvalue <- "Normalization and PCA analysis"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  })       
  
  
  
  ##########################Tab1.3###############################      
  ##############multiple Normalization & PCA###################      
  datainput_multiple_normalization_pca_level <- eventReactive(input$multiple_sample_normalization,{
    source("scripts/multiple_normalization_pca.R")
    datainput_multiple_normalization_pca(index_multiple_normalization_pca_input = datainput_multiple_qc_filter_level()[[7]], index_multiple_sample_normalization_method = input$multiple_sample_normalization_method, index_multiple_sample_scale_factor=input$multiple_sample_scale_factor, index_multiple_sample_var_genes = input$multiple_sample_var_genes, index_multiple_sample_normalization_variable_genes=input$multiple_sample_normalization_variable_genes, index_multiple_sample_pca_dim=input$multiple_sample_pca_dim)
  })
  
  
  output$m_pca_plot<-renderPlot({
    datainput_multiple_normalization_pca_level()[1]
  })
  
  observeEvent(input$download_m_pca_plot, {
    showModal(modalDialog(
      title = strong("Download PCA Plot"),
      numericInput("m_pca_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_pca_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_pca_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_pca_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_pca_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_pca_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_PCA_plot", input$m_pca_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_normalization_pca_level()[[1]], width = input$m_pca_plot_width, height = input$m_pca_plot_height, dpi = input$m_pca_plot_dpi, units = "in")
    }
  )
  
  output$m_elbow_plot<-renderPlot({
    datainput_multiple_normalization_pca_level()[2]
  })
  
  observeEvent(input$download_m_elbow_plot, {
    showModal(modalDialog(
      title = strong("Download Variable Features Plot"),
      numericInput("m_elbow_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_elbow_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_elbow_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_elbow_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_elbow_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_elbow_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_Elbow", input$m_elbow_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_normalization_pca_level()[[2]], width = input$m_elbow_plot_width, height = input$m_elbow_plot_height, dpi = input$m_elbow_plot_dpi, units = "in")
    }
  )
  
  output$m_pca2_plot<-renderPlot({
    datainput_multiple_normalization_pca_level()[3]
  })
  
  observeEvent(input$download_m_pca2_plot, {
    showModal(modalDialog(
      title = strong("Download PCA Plot"),
      numericInput("m_pca2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_pca2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_pca2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_pca2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_pca2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_pca2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_PCA_plot_sample_based", input$m_pca2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_normalization_pca_level()[[3]], width = input$m_pca2_plot_width, height = input$m_pca2_plot_height, dpi = input$m_pca2_plot_dpi, units = "in")
    }
  )
  
  
  output$m_pca3_plot<-renderPlot({
    datainput_multiple_normalization_pca_level()[4]
  })
  
  
  observeEvent(input$download_m_pca3_plot, {
    showModal(modalDialog(
      title = strong("Download PCA Plot"),
      numericInput("m_pca3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_pca3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_pca3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_pca3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_pca3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_pca3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_PCA_plot_group_based", input$m_pca3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_normalization_pca_level()[[4]], width = input$m_pca3_plot_width, height = input$m_pca3_plot_height, dpi = input$m_pca3_plot_dpi, units = "in")
    }
  )
  
  
  
  ###################save seurat object after normalization###################
  output$m_normalization<- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_normalization.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_normalization_pca_level()[[5]], file= file)
    }
  )
  
  #####################################link to next tab###########################      
  observeEvent(input$link_m_clustering, {
    newvalue <- "Clustering"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  })       
  
  
  
  #####################################################Tab1.4####################      
  ########################################multiple Clustering###################      
  datainput_multiple_clustering_level <- eventReactive(input$multiple_sample_clustering,{
    source("scripts/multiple_clustering.R")
    datainput_multiple_clustering(index_multiple_clustering_input = datainput_multiple_normalization_pca_level()[[5]], index_multiple_sample_normalization_method = input$multiple_sample_normalization_method, index_m_clustering1 = input$m_clustering1, index_m_clustering2 = input$m_clustering2, index_m_clustering3 = input$m_clustering3, index_m_clustering4 = input$m_clustering4, index_m_clustering5 = input$m_clustering5, index_m_clustering6 = input$m_clustering6, index_m_clustering7 = input$m_clustering7, index_m_clustering8 = input$m_clustering8, index_m_clustering9 = input$m_clustering9, index_m_clustering10 = input$m_clustering10, index_m_clustering11 = input$m_clustering11, index_m_clustering12 = input$m_clustering12, index_m_clustering13 = input$m_clustering13)
  })
  
  output$m_umap_tsne1_plot<-renderPlot({
    datainput_multiple_clustering_level()[1]
  })
  observeEvent(input$download_m_umap_tsne1_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_umap_tsne1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_umap_tsne1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_umap_tsne1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
    output$m_umap_tsne1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Cluster_plot", input$m_umap_tsne1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clustering_level()[[1]], width = input$m_umap_tsne1_plot_width, height = input$m_umap_tsne1_plot_height, dpi = input$m_umap_tsne1_plot_dpi, units = "in")
    }
  )
  
    
  output$m_umap_tsne_bar1_plot<-renderPlot({
    datainput_multiple_clustering_level()[2]
  }) 
  observeEvent(input$download_m_umap_tsne_bar1_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_umap_tsne_bar1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne_bar1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne_bar1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_umap_tsne_bar1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_umap_tsne_bar1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_umap_tsne_bar1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Cluster_based_bar_plot", input$m_umap_tsne_bar1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clustering_level()[[2]], width = input$m_umap_tsne_bar1_plot_width, height = input$m_umap_tsne_bar1_plot_height, dpi = input$m_umap_tsne_bar1_plot_dpi, units = "in")
    }
  )
  
  
  output$m_umap_tsne2_plot<-renderPlot({
    datainput_multiple_clustering_level()[3]
  })
  observeEvent(input$download_m_umap_tsne2_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_umap_tsne2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_umap_tsne2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_umap_tsne2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_umap_tsne2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Condition_based_plot", input$m_umap_tsne2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clustering_level()[[3]], width = input$m_umap_tsne2_plot_width, height = input$m_umap_tsne2_plot_height, dpi = input$m_umap_tsne2_plot_dpi, units = "in")
    }
  )
  
  
  output$m_umap_tsne_bar2_plot<-renderPlot({
    datainput_multiple_clustering_level()[4]
  })
  observeEvent(input$download_m_umap_tsne_bar2_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_umap_tsne_bar2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne_bar2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne_bar2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_umap_tsne_bar2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_umap_tsne_bar2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_umap_tsne_bar2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Condition_based_bar_plot", input$m_umap_tsne_bar2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clustering_level()[[4]], width = input$m_umap_tsne_bar2_plot_width, height = input$m_umap_tsne_bar2_plot_height, dpi = input$m_umap_tsne_bar2_plot_dpi, units = "in")
    }
  )
  
  
  output$m_umap_tsne3_plot<-renderPlot({
    datainput_multiple_clustering_level()[5]
  })
  observeEvent(input$download_m_umap_tsne3_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_umap_tsne3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_umap_tsne3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_umap_tsne3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_umap_tsne3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Sample_based_plot", input$m_umap_tsne3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clustering_level()[[5]], width = input$m_umap_tsne3_plot_width, height = input$m_umap_tsne3_plot_height, dpi = input$m_umap_tsne3_plot_dpi, units = "in")
    }
  )	  
  
  
  output$m_umap_tsne_bar3_plot<-renderPlot({
    datainput_multiple_clustering_level()[6]
  })
  observeEvent(input$download_m_umap_tsne_bar3_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_umap_tsne_bar3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne_bar3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_umap_tsne_bar3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_umap_tsne_bar3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_umap_tsne_bar3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_umap_tsne_bar3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Sample_based_bar_plot", input$m_umap_tsne_bar3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clustering_level()[[6]], width = input$m_umap_tsne_bar3_plot_width, height = input$m_umap_tsne_bar3_plot_height, dpi = input$m_umap_tsne_bar3_plot_dpi, units = "in")
    }
  )
  
  output$m_clustering_table1<- renderDataTable(DT::datatable((datainput_multiple_clustering_level()[[7]]),
                                                             options = list(
                                                               scrollX = TRUE,
                                                               pageLength = 10,
                                                               bFilter=0
                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_m_clustering_table1 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in clusters", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_clustering_level()[[7]], file)
    }
  ) 
  
  output$m_clustering_table2<- renderDataTable(DT::datatable((datainput_multiple_clustering_level()[[8]]),
                                                             options = list(
                                                               scrollX = TRUE,
                                                               pageLength = 10,
                                                               bFilter=0
                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_m_clustering_table2 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in clusters based on condition", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_clustering_level()[[8]], file)
    }
  ) 
  
  output$m_clustering_table3<- renderDataTable(DT::datatable((datainput_multiple_clustering_level()[[9]]),
                                                             options = list(
                                                               scrollX = TRUE,
                                                               pageLength = 10,
                                                               bFilter=0
                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_m_clustering_table3 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in clusters based on samples", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_clustering_level()[[9]], file)
    }
  ) 
  
  ###################save seurat object after clustering###################
  output$m_clustering<- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_clustering.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_clustering_level()[[10]], file= file)
    }
  )
  
  
  
  #####################################link to next tab###########################   
  observeEvent(input$link_m_doublet, {
    newvalue <- "Remove Doublets"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  }) 
  
  
  ###############################################Tab1.5##################################      
  ####################################multiple Doublet#####################################      
  datainput_multiple_doublet_level <- eventReactive(input$multiple_sample_doublet,{
    source("scripts/multiple_doublet.R")
    datainput_multiple_doublet(index_multiple_doublet_input = datainput_multiple_clustering_level()[[10]], index_m_doublet1 = input$m_doublet1, index_m_doublet3 = input$m_doublet3, index_m_clustering6 = input$m_clustering6, index_multiple_sample_normalization_method = input$multiple_sample_normalization_method, index_multiple_sample_pca_dim = input$multiple_sample_pca_dim)
  })
  
  output$m_doublet1_plot<-renderPlot({
    datainput_multiple_doublet_level()[1]
  })
  observeEvent(input$download_m_doublet1_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_", input$m_clustering6, input$m_doublet1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet_level()[[1]], width = input$m_doublet1_plot_width, height = input$m_doublet1_plot_height, dpi = input$m_doublet1_plot_dpi, units = "in")
    }
  )
  
  output$m_doublet2_plot<-renderPlot({
    datainput_multiple_doublet_level()[2]
  })
  observeEvent(input$download_m_doublet2_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_condition_based_", input$m_clustering6, input$m_doublet2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet_level()[[2]], width = input$m_doublet2_plot_width, height = input$m_doublet2_plot_height, dpi = input$m_doublet2_plot_dpi, units = "in")
    }
  )
  
  output$m_doublet3_plot<-renderPlot({
    datainput_multiple_doublet_level()[3]
  })
  observeEvent(input$download_m_doublet3_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_doublet3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_sample_based_", input$m_clustering6, input$m_doublet3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet_level()[[3]], width = input$m_doublet3_plot_width, height = input$m_doublet3_plot_height, dpi = input$m_doublet3_plot_dpi, units = "in")
    }
  )
  
  
  output$m_doublet4_plot<-renderPlot({
    datainput_multiple_doublet_level()[4]
  })
  
  observeEvent(input$download_m_doublet4_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet4_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet4_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_cluster_based_", input$m_clustering6, input$m_doublet4_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet_level()[[4]], width = input$m_doublet4_plot_width, height = input$m_doublet4_plot_height, dpi = input$m_doublet4_plot_dpi, units = "in")
    }
  )
  
  output$m_doublet1_table<- renderDataTable(DT::datatable((datainput_multiple_doublet_level()[[5]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_m_doublet1_table <- downloadHandler(
    filename = function() {
      paste("Number_of_singlets_and_doublet_counts", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_doublet_level()[[5]], file)
    }
  )
  
  
  #############submenu#########    
  datainput_multiple_doublet2_level <- eventReactive(input$multiple_sample_doublet2,{
    source("scripts/multiple_doublet1.R")
    datainput_multiple_doublet2(index_multiple_doublet1_input = datainput_multiple_doublet_level()[[6]], index_m_doublet2 = input$m_doublet2, index_m_doublet3 = input$m_doublet3, index_m_clustering6 = input$m_clustering6, index_multiple_sample_normalization_method = input$multiple_sample_normalization_method)
  })
  
  
  output$m_doublet5_plot<-renderPlot({
    datainput_multiple_doublet2_level()[1]
  })
  observeEvent(input$download_m_doublet5_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet5_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_doublet5_plot_width", label = h5("Figure width (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_doublet5_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet5_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet5_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet5_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_after_",input$m_doublet2,"_", input$m_clustering6, input$m_doublet5_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet2_level()[[1]], width = input$m_doublet5_plot_width, height = input$m_doublet5_plot_height, dpi = input$m_doublet5_plot_dpi, units = "in")
    }
  )
  
  
  output$m_doublet6_plot<-renderPlot({
    datainput_multiple_doublet2_level()[2]
  })
  
  observeEvent(input$download_m_doublet6_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet6_plot_height", label = h5("Figure height (upto 49 inces)"), value = 10, width = "300px"),
      numericInput("m_doublet6_plot_width", label = h5("Figure width (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_doublet6_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet6_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet6_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet6_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_condition_based_after_",input$m_doublet2,"_", input$m_clustering6, input$m_doublet6_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet2_level()[[2]], width = input$m_doublet6_plot_width, height = input$m_doublet6_plot_height, dpi = input$m_doublet6_plot_dpi, units = "in")
    }
  )
  
  output$m_doublet7_plot<-renderPlot({
    datainput_multiple_doublet2_level()[3]
  })
  
  observeEvent(input$download_m_doublet7_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet7_plot_height", label = h5("Figure height (upto 49 inces)"), value = 10, width = "300px"),
      numericInput("m_doublet7_plot_width", label = h5("Figure width (upto 49 inces)"), value = 14, width = "300px"),
      numericInput("m_doublet7_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet7_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet7_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet7_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_sample_based_after_",input$m_doublet2,"_", input$m_clustering6, input$m_doublet7_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet2_level()[[3]], width = input$m_doublet7_plot_width, height = input$m_doublet7_plot_height, dpi = input$m_doublet7_plot_dpi, units = "in")
    }
  )
  
  output$m_doublet8_plot<-renderPlot({
    datainput_multiple_doublet2_level()[4]
  })   
  observeEvent(input$download_m_doublet8_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet8_plot_height", label = h5("Figure height (upto 49 inces)"), value = 10, width = "300px"),
      numericInput("m_doublet8_plot_width", label = h5("Figure width (upto 49 inces)"), value = 14, width = "300px"),
      numericInput("m_doublet8_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet8_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet8_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet8_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_cluster_based_after_",input$m_doublet2,"_", input$m_clustering6, input$m_doublet8_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet2_level()[[4]], width = input$m_doublet8_plot_width, height = input$m_doublet8_plot_height, dpi = input$m_doublet8_plot_dpi, units = "in")
    }
  )
  
  output$m_doublet9_plot<-renderPlot({
    datainput_multiple_doublet2_level()[5]
  })  
  
  observeEvent(input$download_m_doublet9_plot, {
    showModal(modalDialog(
      title = strong("Download bar Plot"),
      numericInput("m_doublet9_plot_height", label = h5("Figure height (upto 49 inces)"), value = 6, width = "300px"),
      numericInput("m_doublet9_plot_width", label = h5("Figure width (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_doublet9_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet9_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet9_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_doublet9_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Bar_plot_cluster_based_after_",input$m_doublet2,"_", input$m_clustering6, input$m_doublet9_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet2_level()[[5]], width = input$m_doublet9_plot_width, height = input$m_doublet9_plot_height, dpi = input$m_doublet9_plot_dpi, units = "in")
    }
  )
  
  
  
  
  output$m_doublet2_table<- renderDataTable(DT::datatable((datainput_multiple_doublet2_level()[[6]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_m_doublet2_table <- downloadHandler(
    filename = function() {
      paste("Number_of_cell_after_", input$m_doublet2, '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_doublet2_level()[[6]], file)
    }
  )
  
  output$m_doublet3_table<- renderDataTable(DT::datatable((datainput_multiple_doublet2_level()[[7]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_m_doublet3_table <- downloadHandler(
    filename = function() {
      paste("Number_of_cell_after_", input$m_doublet2, '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_doublet2_level()[[7]], file)
    }
  )
  
  output$m_doublet4_table<- renderDataTable(DT::datatable((datainput_multiple_doublet2_level()[[8]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_m_doublet4_table <- downloadHandler(
    filename = function() {
      paste("Number_of_cell_after_", input$m_doublet2, '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_doublet2_level()[[8]], file)
    }
  )
  
  output$m_doublet5_table<- renderDataTable(DT::datatable((datainput_multiple_doublet2_level()[[9]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_m_doublet5_table <- downloadHandler(
    filename = function() {
      paste("Number_of_cell_after_", input$m_doublet2, '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_doublet2_level()[[9]], file)
    }
  )
  
  
  output$m_doublet10_plot<-renderPlot({
    datainput_multiple_doublet2_level()[11]
  })  
  
  observeEvent(input$download_m_doublet10_plot, {
    showModal(modalDialog(
      title = strong("Download Plot"),
      numericInput("m_doublet10_plot_height", label = h5("Figure height (upto 49 inces)"), value = 10, width = "300px"),
      numericInput("m_doublet10_plot_width", label = h5("Figure width (upto 49 inces)"), value = 14, width = "300px"),
      numericInput("m_doublet10_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_doublet10_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_doublet10_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_doublet10_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Doublet_plot_cluster_based_split_by_condition_and_sample_after_",input$m_doublet2,"_", input$m_clustering6, input$m_doublet10_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_doublet2_level()[[11]], width = input$m_doublet10_plot_width, height = input$m_doublet10_plot_height, dpi = input$m_doublet10_plot_dpi, units = "in")
    }
  )
  ###################save seurat object after doublet removal###################
  output$m_doublet<- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_doublet_removal.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_doublet2_level()[[10]], file= file)
    }
  )
  
  #####################################link to next tab###########################     
  observeEvent(input$link_m_marker, {
    newvalue <- "Markers Identification"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  })   
  
  
  ##########################Tab1.6###############################      
  ##############multiple Marker identification###################    
  observeEvent(input$multiple_sample_marker,{
    if(input$m_marker1 == 1){
      shinyjs::show("m_marker_box5")
      shinyjs::show("m_marker_box6")
      shinyjs::show("m_marker_box7")
    }
    else{
      shinyjs::show("m_marker_box5")
      shinyjs::hide("m_marker_box6")
      shinyjs::show("m_marker_box7")
    }
    
  })
  
  output$m_marker_6 <- renderUI ({
    clusters <- req(datainput_multiple_clustering_level()[[11]])
    
    shinyWidgets::pickerInput(
      inputId = "m_marker6",
      label = "Select one cluster for analsysis",
      choices = sort(clusters),
      multiple = F,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_marker_7 <- renderUI ({
    clusters <- req(datainput_multiple_clustering_level()[[11]])
    clusters <- clusters[!clusters == input$m_marker6]
    shinyWidgets::pickerInput(
      inputId = "m_marker7",
      label = "Identify markers distinguishing a cluster from other selected clusters",
      choices = sort(clusters),
      multiple = T,
      selected = sort(clusters)[1],
      options = list(`actions-box` = TRUE))
  })
  
  output$m_marker_8 <- renderUI ({
    clusters <- req(datainput_multiple_clustering_level()[[11]])
    
    shinyWidgets::pickerInput(
      inputId = "m_marker8",
      label = "Select one cluster to define markers",
      choices = sort(clusters),
      multiple = F,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_marker_9 <- renderUI ({
    clusters <- req(datainput_multiple_clustering_level()[[11]])
    clusters <- clusters[!clusters == input$m_marker8]
    shinyWidgets::pickerInput(
      inputId = "m_marker9",
      label = "Select the cluster to find the conserved markers between two clusters",
      choices = sort(clusters),
      multiple = T,
      selected = sort(clusters)[1],
      options = list(`actions-box` = TRUE))
  })
  
  
  datainput_multiple_marker_level <- eventReactive(input$multiple_sample_marker,{
    
    source("scripts/multiple_marker.R")
    datainput_multiple_marker(index_multiple_marker_input = datainput_multiple_doublet2_level()[[10]], index_m_marker1 = input$m_marker1, index_m_marker2 = input$m_marker2, index_m_marker3 = input$m_marker3, index_m_marker4 = input$m_marker4, index_m_marker5 = input$m_marker5, index_m_marker6 = input$m_marker6, index_m_marker7 = input$m_marker7, index_m_marker8 = input$m_marker8, index_m_marker9 = input$m_marker9, index_m_marker10 = input$m_marker10, index_multiple_sample_normalization_method = input$multiple_sample_normalization_method)
  })
  
  
  output$m_marker1_table<- renderDataTable(DT::datatable((datainput_multiple_marker_level()[[1]]),
                                                         options = list(
                                                           scrollX = TRUE,
                                                           pageLength = 10,
                                                           dom = "Blfrtip"
                                                           #bFilter=0
                                                         ),rownames= FALSE, selection = "none"))
  
  
  output$download_m_marker1_table <- downloadHandler(
    filename = function() {
      paste("Number_of_identified_markers_or_differentially_expressed_genes", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_marker_level()[[1]], file)
    }
  )
  
  
  output$m_marker1_plot<-renderPlot({
    datainput_multiple_marker_level()[3]
  })
  
  
  observeEvent(input$download_m_marker1_plot, {
    showModal(modalDialog(
      title = strong("Download Heatmap"),
      numericInput("m_marker1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_marker1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_marker1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_marker1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_marker1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_marker1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste( "Heatmap_with_Top5_expressed_genes", input$m_marker1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_marker_level()[[3]], width = input$m_marker1_plot_width, height = input$m_marker1_plot_height, dpi = input$m_marker1_plot_dpi, units = "in")
    }
  )
  
  
  ###################save seurat object after doublet removal###################
  output$m_marker<- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_marker_identification.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_marker_level()[[2]], file= file)
    }
  )
  
  #####################################link to next tab###########################     
  observeEvent(input$link_m_prediction, {
    newvalue <- "Cell Type Prediction"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  })  
  
  ##########################Tab1.7###############################      
  ##############multiple cell type################### 
  output$m_celltype7 <- renderUI({
    numberofclusters <- as.integer(length(levels(datainput_multiple_marker_level()[[2]])))
    lapply(1:numberofclusters, function(i) {
      column(3, textInput(paste("mcelltypenames", levels(datainput_multiple_marker_level()[[2]])[i], sep = ""),
                          paste("Cluster", levels(datainput_multiple_marker_level()[[2]])[i]), value = paste("Cluster", levels(datainput_multiple_marker_level()[[2]])[i])))
    })
  })
  
  observeEvent(input$multiple_sample_celltype,{
    if(input$m_celltype1 == 1){
      shinyjs::show("m_celltype_box7")
      shinyjs::show("m_celltype_box8")
      shinyjs::hide("m_celltype_box9")
      shinyjs::show("m_celltype_box10")
      shinyjs::show("m_celltype9")
      shinyjs::hide("m_celltype10")
    }
    else if (input$m_celltype1 == 2){
      shinyjs::show("m_celltype_box7")
      shinyjs::show("m_celltype_box8")
      shinyjs::show("m_celltype_box9")
      shinyjs::show("m_celltype_box10")
      shinyjs::hide("m_celltype9")
      shinyjs::show("m_celltype10")
    }
    else if (input$m_celltype1 == 3){
      shinyjs::show("m_celltype_box7")
      shinyjs::hide("m_celltype_box8")
      shinyjs::hide("m_celltype_box9")
      shinyjs::show("m_celltype_box10")
      shinyjs::hide("m_celltype9")
      shinyjs::hide("m_celltype10")
    }
    else if (input$m_celltype1 == 4){
      shinyjs::show("m_celltype_box7")
      shinyjs::hide("m_celltype_box8")
      shinyjs::hide("m_celltype_box9")
      shinyjs::show("m_celltype_box10")
      shinyjs::hide("m_celltype9")
      shinyjs::hide("m_celltype10")
    }
    
  })
  
  
  datainput_multiple_celltype_level <- eventReactive(input$multiple_sample_celltype,{
    source("scripts/multiple_celltype.R")
    datainput_multiple_celltype(index_multiple_celltype_input = datainput_multiple_marker_level()[[2]], index_cell_markers = datainput_multiple_marker_level()[[1]], index_m_celltype1 = input$m_celltype1, index_m_celltype2 = input$m_celltype2, index_m_celltype3 = input$m_celltype3, index_m_celltype4 = input$m_celltype4, index_m_celltype5 = input$m_celltype5, index_m_celltype6 = input$m_celltype6, index_m_celltype7 = c(input$mcelltypenames0,input$mcelltypenames1,input$mcelltypenames2,input$mcelltypenames3,input$mcelltypenames4,input$mcelltypenames5,input$mcelltypenames6,input$mcelltypenames7,input$mcelltypenames8,input$mcelltypenames9,input$mcelltypenames10,input$mcelltypenames11,input$mcelltypenames12,input$mcelltypenames13,input$mcelltypenames14,input$mcelltypenames15,input$mcelltypenames16,input$mcelltypenames17,input$mcelltypenames18,input$mcelltypenames19,input$mcelltypenames20,input$mcelltypenames21,input$mcelltypenames22,input$mcelltypenames23,input$mcelltypenames24,input$mcelltypenames25,input$mcelltypenames26,input$mcelltypenames27,input$mcelltypenames28,input$mcelltypenames29,input$mcelltypenames30,input$mcelltypenames31,input$mcelltypenames32,input$mcelltypenames33,input$mcelltypenames34,input$mcelltypenames35,input$mcelltypenames36,input$mcelltypenames37,input$mcelltypenames38,input$mcelltypenames39,input$mcelltypenames40,input$mcelltypenames41,input$mcelltypenames42,input$mcelltypenames43,input$mcelltypenames44,input$mcelltypenames45,input$mcelltypenames46,input$mcelltypenames47,input$mcelltypenames48,input$mcelltypenames49,input$mcelltypenames50,input$mcelltypenames51,input$mcelltypenames52,input$mcelltypenames53,input$mcelltypenames54,input$mcelltypenames55,input$mcelltypenames56,input$mcelltypenames57,input$mcelltypenames58,input$mcelltypenames59,input$mcelltypenames60,input$mcelltypenames61,input$mcelltypenames62,input$mcelltypenames63,input$mcelltypenames64,input$mcelltypenames65,input$mcelltypenames66,input$mcelltypenames67,input$mcelltypenames68,input$mcelltypenames69,input$mcelltypenames70,input$mcelltypenames71,input$mcelltypenames72,input$mcelltypenames73,input$mcelltypenames74,input$mcelltypenames75,input$mcelltypenames76,input$mcelltypenames77,input$mcelltypenames78,input$mcelltypenames79,input$mcelltypenames80,input$mcelltypenames81,input$mcelltypenames82,input$mcelltypenames83,input$mcelltypenames84,input$mcelltypenames85,input$mcelltypenames86,input$mcelltypenames87,input$mcelltypenames88,input$mcelltypenames89,input$mcelltypenames90,input$mcelltypenames91,input$mcelltypenames92,input$mcelltypenames93,input$mcelltypenames94,input$mcelltypenames95,input$mcelltypenames96,input$mcelltypenames97,input$mcelltypenames98,input$mcelltypenames99), index_m_celltype8 = input$m_celltype8, index_m_clustering6 = input$m_clustering6, index_multiple_sample_normalization_method = input$multiple_sample_normalization_method)
  })
  output$m_celltype1_plot<-renderPlot({
    datainput_multiple_celltype_level()[5]
  })
  
  observeEvent(input$download_m_celltype1_plot, {
    showModal(modalDialog(
      title = strong("Download Celltype"),
      numericInput("m_celltype1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_celltype1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 16, width = "300px"),
      numericInput("m_celltype1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_celltype1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_celltype1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_celltype1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Dimplot_with_celltype", input$m_celltype1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_celltype_level()[[5]], width = input$m_celltype1_plot_width, height = input$m_celltype1_plot_height, dpi = input$m_celltype1_plot_dpi, units = "in")
    }
  )
  
  output$m_celltype1_table<- renderDataTable(DT::datatable((datainput_multiple_celltype_level()[[6]]),
                                                           options = list(
                                                             scrollX = TRUE,
                                                             pageLength = 10,
                                                             dom = "Blfrtip"
                                                             #bFilter=0
                                                           ),rownames= TRUE, selection = "none"))
  
  output$download_m_celltype1_table <- downloadHandler(
    filename = function() { 
      paste("predicted_celltype_Scores", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_celltype_level()[[6]], file)
    }
  )
  
  output$m_celltype2_plot<-renderPlot({
    datainput_multiple_celltype_level()[7]
  })
  
  observeEvent(input$download_m_celltype2_plot, {
    showModal(modalDialog(
      title = strong("Download Heatmap"),
      numericInput("m_celltype2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_celltype2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_celltype2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_celltype2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_celltype2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_celltype2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("score_plots", input$m_celltype2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_celltype_level()[[7]], width = input$m_celltype2_plot_width, height = input$m_celltype2_plot_height, dpi = input$m_celltype2_plot_dpi, units = "in")
    }
  )
  
  
  output$m_celltype3_plot<-renderPlot({
    datainput_multiple_celltype_level()[8]
  })
  
  observeEvent(input$download_m_celltype3_plot, {
    showModal(modalDialog(
      title = strong("Download"),
      numericInput("m_celltype3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_celltype3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_celltype3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_celltype3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_celltype3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_celltype3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("score_plots", input$m_celltype3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_celltype_level()[[8]], width = input$m_celltype3_plot_width, height = input$m_celltype3_plot_height, dpi = input$m_celltype3_plot_dpi, units = "in")
    }
  )
  
 
  
  ###################save seurat object after doublet###################
  output$m_celltype<- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_celltypes.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_celltype_level()[[1]], file= file)
    }
  )
  #####################################link to next tab###########################     
  observeEvent(input$link_m_clusterbased, {
    newvalue <- "Cluster-based plots"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  })  
  
  
  ###################################Tab1.8###############################  
  ###########################Cluster-based plots####################### 
  observe({
    if(input$m_clusterbased1 == "gene_name_list"){
      shinyjs::show("m_clusterbased2")
    }
    else {
      shinyjs::hide("m_clusterbased2")
    }
  })  
  
  observeEvent(input$multiple_sample_clusterbased, {
    if(input$m_clusterbased4 ==  "seurat_clusters"){
      shinyjs::show("m_clusterbased_box3")
    }
    else {
      shinyjs::hide("m_clusterbased_box3")
    }
  })  
  
  observe({
    if(input$m_clusterbased3 ==  "RidgePlot"){
      shinyjs::hide("m_clusterbased5")
    }
    else {
      shinyjs::show("m_clusterbased5")
    }
  })
  
  observeEvent(input$multiple_sample_clusterbased,{
    shinyjs::show("m_clusterbased_box2")
    #shinyjs::show("m_clusterbased_box3")
    shinyjs::show("m_clusterbased_box4")
  })
  
  
  
  datainput_multiple_clusterbased_level <- eventReactive(input$multiple_sample_clusterbased,{
    source("scripts/multiple_clusterbased.R")
    datainput_multiple_clusterbased(index_multiple_clusterbased_input = datainput_multiple_celltype_level()[[1]], index_multiple_clusterbased_features = datainput_multiple_marker_level()[[1]], index_m_celltype_method = datainput_multiple_celltype_level()[[4]], index_m_clusterbased1 = input$m_clusterbased1, index_m_clusterbased2 = input$m_clusterbased2, index_m_clusterbased3 = input$m_clusterbased3, index_m_clusterbased4 = input$m_clusterbased4, index_m_clusterbased5 = input$m_clusterbased5)
  })  
  
  output$m_clusterbased1_plot<-renderPlot({
    datainput_multiple_clusterbased_level()[1]
  })
  
  observeEvent(input$download_m_clusterbased1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("m_clusterbased1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_clusterbased1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_clusterbased1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_clusterbased1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_clusterbased1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_clusterbased1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Plots_for_top_or_selected_markers",  input$m_clusterbased1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_clusterbased_level()[[1]], width = input$m_clusterbased1_plot_width, height = input$m_clusterbased1_plot_height, dpi = input$m_clusterbased1_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_clusterbased1_table<- renderDataTable(DT::datatable((datainput_multiple_clusterbased_level()[[3]]),
                                                               options = list(
                                                                 scrollX = TRUE,
                                                                 pageLength = 10,
                                                                 dom = "Blfrtip"
                                                                 #bFilter=0
                                                               ),rownames= FALSE, selection = "none"))
  
  output$download_m_clusterbased1_table <- downloadHandler(
    filename = function() { 
      paste("Top_or_selected_Cell_counts_proportion", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_clusterbased_level()[[3]], file)
    }
  )
  
  
  ###################save seurat object after doublet###################
  output$m_clusterbased <- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_plots.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_clusterbased_level()[[2]], file= file)
    }
  )
  
  #####################################link to next tab###########################     
  observeEvent(input$link_m_conditionbased, {
    newvalue <- "Condition based analysis"
    updateTabsetPanel(session, "multiple_tabsets", newvalue)
  })  
  
  
  ##########################Tab1.9###############################      
  ##############Condition-based analysis###################    
  observe({
    if(input$m_conditionbased7 == "VolcanoPlot"){
      shinyjs::hide("m_conditionbased8")
      shinyjs::hide("m_conditionbased9")
      #shinyjs::hide("m_conditionbased_box4")
    }
    else{
      shinyjs::show("m_conditionbased8")
      shinyjs::show("m_conditionbased9")
      #shinyjs::show("m_conditionbased_box4")
    }
  })
  
  observeEvent(input$multiple_sample_conditionbased,{
    shinyjs::show("m_conditionbased_box3")
    shinyjs::show("m_conditionbased_box4")
    shinyjs::show("m_conditionbased_box5")
  })
  
  
  
  output$m_conditionbased_1 <- renderUI ({
    clusters <- req(datainput_multiple_clustering_level()[[13]])
    
    shinyWidgets::pickerInput(
      inputId = "m_conditionbased1",
      label = "Select the Condition1",
      choices = sort(clusters),
      multiple = F,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_conditionbased_2 <- renderUI ({
    clusters <- req(datainput_multiple_clustering_level()[[13]])
    clusters <- clusters[!clusters == input$m_conditionbased1]
    shinyWidgets::pickerInput(
      inputId = "m_conditionbased2",
      label = "Select the Condition2",
      choices = sort(clusters),
      selected = sort(clusters)[1],
      multiple = F,
      options = list(`actions-box` = TRUE))
  })  
  
  
  
  datainput_multiple_conditionbased_level <- eventReactive(input$multiple_sample_conditionbased,{
    source("scripts/multiple_conditionbased.R")
    datainput_multiple_conditionbased(index_multiple_conditionbased_input = datainput_multiple_celltype_level()[[1]], index_multiple_sample_normalization_method = input$multiple_sample_normalization_method, index_m_conditionbased1 = input$m_conditionbased1, index_m_conditionbased2 = input$m_conditionbased2, index_m_conditionbased3 = input$m_conditionbased3, index_m_conditionbased4 = input$m_conditionbased4, index_m_conditionbased5 = input$m_conditionbased5, index_m_conditionbased6 = input$m_conditionbased6, index_m_conditionbased7 = input$m_conditionbased7, index_m_conditionbased8 = input$m_conditionbased8, index_m_conditionbased9 = input$m_conditionbased9)
  })  
  
  output$m_conditionbased1_plot<-renderPlot({
    datainput_multiple_conditionbased_level()[1]
  })
  
  observeEvent(input$download_m_conditionbased1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("m_conditionbased1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_conditionbased1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_conditionbased1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_conditionbased1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_conditionbased1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_conditionbased1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Plots_for_top_selected_markers",  input$m_conditionbased1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_multiple_conditionbased_level()[[1]], width = input$m_conditionbased1_plot_width, height = input$m_conditionbased1_plot_height, dpi = input$m_conditionbased1_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_conditionbased1_table<- renderDataTable(DT::datatable((datainput_multiple_conditionbased_level()[[2]]),
                                                                 options = list(
                                                                   scrollX = TRUE,
                                                                   pageLength = 10,
                                                                   dom = "Blfrtip"
                                                                   #bFilter=0
                                                                 ),rownames= FALSE, selection = "none"))
  
  output$download_m_conditionbased1_table <- downloadHandler(
    filename = function() { 
      paste("Top_or_selected_Cell_counts_proportion", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_multiple_conditionbased_level()[[2]], file)
    }
  )
  
  
  ###################save seurat object after doublet###################
  output$m_conditionbased <- downloadHandler(
    filename = function(){
      paste("multiple_sample_seuart_object_after_plots.RDS")
    },
    content = function(file){
      saveRDS(datainput_multiple_conditionbased_level()[[3]], file= file)
    }
  )
  
  
 
#################################################################Tab2#################################################################
#################################################Multiple samples subclustering######################################################      
  
  ##########multiple sidebar hide##########
  observeEvent(input[["subclustering_multiple_tabsets"]], {
    if(input[["subclustering_multiple_tabsets"]] == "Cell Stats"){
      showElement(selector = "#subclustering_multiple_sidebar")
      removeCssClass("subclustering_multiple_main_menu", "col-sm-12")
      addCssClass("subclustering_multiple_main_menu", "col-sm-8")
    }else{
      hideElement(selector = "#subclustering_multiple_sidebar")
      removeCssClass("subclustering_multiple_main_menu", "col-sm-8")
      addCssClass("subclustering_multiple_main_menu", "col-sm-12")
    }
  })
  
  ############subclustering multiple hide tabs########
  hideTab(inputId = "subclustering_multiple_tabsets", target = "Normalization and PCA analysis")
  hideTab(inputId = "subclustering_multiple_tabsets", target = "Clustering")
  hideTab(inputId = "subclustering_multiple_tabsets", target = "Markers Identification")
  hideTab(inputId = "subclustering_multiple_tabsets", target = "Cell Type Prediction")
  hideTab(inputId = "subclustering_multiple_tabsets", target = "Cluster-based plots")
  hideTab(inputId = "subclustering_multiple_tabsets", target = "Condition based analysis")
  
  
  ################hide menu bar######################
  #hideTab(inputId = "menu_tabs", target = "Multiple samples subclustering")
  
  ##########subclustering showtabbutton#############     
  observeEvent(input$link_m_subclustering_normalization, {
    showTab(inputId = "subclustering_multiple_tabsets", target = "Normalization and PCA analysis")
  })
  
  observeEvent(input$link_m_subclustering_clustering, {
    showTab(inputId = "subclustering_multiple_tabsets", target = "Clustering")
  })

  observeEvent(input$link_m_subclustering_marker, {
    showTab(inputId = "subclustering_multiple_tabsets", target = "Markers Identification")
  })
  
  observeEvent(input$link_m_subclustering_prediction, {
    showTab(inputId = "subclustering_multiple_tabsets", target = "Cell Type Prediction")
  })
  
  observeEvent(input$link_m_subclustering_clusterbased, {
    showTab(inputId = "subclustering_multiple_tabsets", target = "Cluster-based plots")
  })    
  
  observeEvent(input$link_m_subclustering_conditionbased, {
    showTab(inputId = "subclustering_multiple_tabsets", target = "Condition based analysis")
  })
  
  ################show menu bar######################
  observeEvent(input$multiple_sample_celltype, {
    #showTab(inputId = "menu_tabs", target = "Multiple samples subclustering")
    shinyjs::hide("m_subclustering0")
    shinyjs::show("m_subclustering1")
    shinyjs::show("subclustering_multiple_sample_submit")
    })   
  
  ########multiple samples subclustering hide stats########  
  shinyjs::hide("m_subclustering1")
  shinyjs::hide("m_subclustering_2")
  shinyjs::hide("m_subclustering_3")
  shinyjs::hide("m_subclustering3")
  shinyjs::hide("subclustering_multiple_sample_submit")
  shinyjs::hide("m_subclustering_box1")
  shinyjs::hide("m_subclustering_box2")
  shinyjs::hide("m_subclustering_box3")
  
  ##########multiple sidebar hide##########
  observeEvent(input[["subclustering_multiple_tabsets"]], {
    if(input[["subclustering_multiple_tabsets"]] == "Cell Stats"){
      showElement(selector = "#subclustering_multiple_sidebar")
      removeCssClass("subclustering_multiple_main_menu", "col-sm-12")
      addCssClass("subclustering_multiple_main_menu", "col-sm-8")
    }else{
      hideElement(selector = "#subclustering_multiple_sidebar")
      removeCssClass("subclustering_multiple_main_menu", "col-sm-8")
      addCssClass("subclustering_multiple_main_menu", "col-sm-12")
    }
  })
  
  
  ########multiple hide normalization########     
  shinyjs::hide("m_subclustering_pca_box1")
  shinyjs::hide("m_subclustering_elbow_box")
  shinyjs::hide("m_subclustering_pca_box2")
  shinyjs::hide("m_subclustering_pca_box3")
  shinyjs::hide("m_subclustering_pca_box4")
  
  observeEvent(input$subclustering_multiple_sample_normalization,{
    shinyjs::show("m_subclustering_pca_box1")
    shinyjs::show("m_subclustering_elbow_box")
    shinyjs::show("m_subclustering_pca_box2")
    shinyjs::show("m_subclustering_pca_box3")
    shinyjs::show("m_subclustering_pca_box4")
  })
  
  observe({
    if (input$subclustering_multiple_sample_normalization_method == "LogNormalize") {
      shinyjs::show("subclustering_multiple_sample_scale_factor")
      shinyjs::show("subclustering_multiple_sample_normalization_variable_genes")
      shinyjs::show("subclustering_multiple_sample_var_genes")
    }
    else if (input$subclustering_multiple_sample_normalization_method  == "SCTransform") {
      shinyjs::hide("subclustering_multiple_sample_scale_factor")
      shinyjs::hide("subclustering_multiple_sample_normalization_variable_genes")
      shinyjs::hide("subclustering_multiple_sample_var_genes")
    }
  })
  
  ########multiple hide clustering########     
  shinyjs::hide("m_subclustering_clustering_box1")
  shinyjs::hide("m_subclustering_clustering_box2")
  shinyjs::hide("m_subclustering_clustering_box3")
  shinyjs::hide("m_subclustering_clustering_box4")
  shinyjs::hide("m_subclustering_clustering_box5")
  shinyjs::hide("m_subclustering_clustering_box6") 
  shinyjs::hide("m_subclustering_clustering_box7") 
  shinyjs::hide("m_subclustering_clustering_box8")   
  shinyjs::hide("m_subclustering_clustering_box9") 
  shinyjs::hide("m_subclustering_clustering_box10") 
  shinyjs::hide("m_subclustering_clustering_box11")
  
  observeEvent(input$subclustering_multiple_sample_clustering,{
    shinyjs::show("m_subclustering_clustering_box1")
    shinyjs::show("m_subclustering_clustering_box2")
    shinyjs::show("m_subclustering_clustering_box3")
    shinyjs::show("m_subclustering_clustering_box4")
    shinyjs::show("m_subclustering_clustering_box5")
    shinyjs::show("m_subclustering_clustering_box6") 
    shinyjs::show("m_subclustering_clustering_box7") 
    shinyjs::show("m_subclustering_clustering_box8") 
    shinyjs::show("m_subclustering_clustering_box9") 
    shinyjs::show("m_subclustering_clustering_box10") 
    shinyjs::show("m_subclustering_clustering_box11") 
  })
  
  
  observe({
    if (input$m_subclustering_clustering6 == "umap") {
      shinyjs::show("m_subclustering_umap_box")
      shinyjs::hide("m_subclustering_tsne_box")
    }
    else if (input$m_subclustering_clustering6  == "tsne") {
      shinyjs::hide("m_subclustering_umap_box")
      shinyjs::show("m_subclustering_tsne_box")
    }
  })
  
  
  ########multiple hide doublet boxes########
  shinyjs::hide("m_subclustering_doublet_box2")
  shinyjs::hide("m_subclustering_doublet_box3")
  shinyjs::hide("m_subclustering_doublet_box4")
  shinyjs::hide("m_subclustering_doublet_box5")
  shinyjs::hide("m_subclustering_doublet_box6")
  shinyjs::hide("m_subclustering_doublet_box7")
  shinyjs::hide("m_subclustering_doublet_box8")
  shinyjs::hide("m_subclustering_doublet_box9")
  shinyjs::hide("m_subclustering_doublet_box10")
  shinyjs::hide("m_subclustering_doublet_box11")
  shinyjs::hide("m_subclustering_doublet_box12")
  shinyjs::hide("m_subclustering_doublet_box13")
  shinyjs::hide("m_subclustering_doublet_box14")
  shinyjs::hide("m_subclustering_doublet_box15")
  shinyjs::hide("m_subclustering_doublet_box16")
  shinyjs::hide("m_subclustering_doublet_box17")
  shinyjs::hide("m_subclustering_doublet_box18")
  
  observeEvent(input$subclustering_multiple_sample_doublet,{
    shinyjs::show("m_subclustering_doublet_box2")
    shinyjs::show("m_subclustering_doublet_box3")
    shinyjs::show("m_subclustering_doublet_box4")
    shinyjs::show("m_subclustering_doublet_box5")
    shinyjs::show("m_subclustering_doublet_box6")
    shinyjs::show("m_subclustering_doublet_box7")
  })
  observeEvent(input$subclustering_multiple_sample_doublet2,{
    
    shinyjs::show("m_subclustering_doublet_box8")
    shinyjs::show("m_subclustering_doublet_box9")
    shinyjs::show("m_subclustering_doublet_box10")
    shinyjs::show("m_subclustering_doublet_box11")
    shinyjs::show("m_subclustering_doublet_box12")
    shinyjs::show("m_subclustering_doublet_box13")
    shinyjs::show("m_subclustering_doublet_box14")
    shinyjs::show("m_subclustering_doublet_box15")
    shinyjs::show("m_subclustering_doublet_box16")
    shinyjs::show("m_subclustering_doublet_box17")
    shinyjs::show("m_subclustering_doublet_box18")
    
  })
  
  
  ########multiple hide markers box########  
  shinyjs::hide("m_subclustering_marker_box5")
  shinyjs::hide("m_subclustering_marker_box6")
  shinyjs::hide("m_subclustering_marker_box7")
  shinyjs::hide("m_subclustering_marker10")
  shinyjs::hide("m_subclustering_marker11")
  shinyjs::hide("m_subclustering_marker12")
  
  
  observe({
    if (input$m_subclustering_marker1 == 1) {
      shinyjs::hide("m_subclustering_marker_6")
      shinyjs::hide("m_subclustering_marker_7")
      shinyjs::hide("m_subclustering_marker_8")
      shinyjs::hide("m_subclustering_marker_9")
      shinyjs::hide("m_subclustering_marker6")
      shinyjs::hide("m_subclustering_marker7")
      shinyjs::hide("m_subclustering_marker8")
      shinyjs::hide("m_subclustering_marker9")
      shinyjs::hide("m_subclustering_marker10")
      shinyjs::show("m_subclustering_marker11")
      shinyjs::hide("m_subclustering_marker12")
    }
    else if (input$m_subclustering_marker1 == 2) {
      shinyjs::show("m_subclustering_marker_6")
      shinyjs::hide("m_subclustering_marker_7")
      shinyjs::hide("m_subclustering_marker_8")
      shinyjs::hide("m_subclustering_marker_9")
      shinyjs::show("m_subclustering_marker6")
      shinyjs::hide("m_subclustering_marker7")
      shinyjs::hide("m_subclustering_marker8")
      shinyjs::hide("m_subclustering_marker9")
      shinyjs::hide("m_subclustering_marker_box6")
      shinyjs::hide("m_subclustering_marker10")
      shinyjs::show("m_subclustering_marker11")
      shinyjs::hide("m_subclustering_marker12")
      
    }
    else if (input$m_subclustering_marker1 == 3) {
      shinyjs::show("m_subclustering_marker_6")
      shinyjs::show("m_subclustering_marker_7")
      shinyjs::hide("m_subclustering_marker_8")
      shinyjs::hide("m_subclustering_marker_9")
      shinyjs::show("m_subclustering_marker6")
      shinyjs::show("m_subclustering_marker7")
      shinyjs::hide("m_subclustering_marker8")
      shinyjs::hide("m_subclustering_marker9")
      shinyjs::hide("m_subclustering_marker_box6")
      shinyjs::hide("m_subclustering_marker10")
      shinyjs::show("m_subclustering_marker11")
      shinyjs::hide("m_subclustering_marker12")
    }
    else if (input$m_subclustering_marker1 == 4) {
      shinyjs::hide("m_subclustering_marker_6")
      shinyjs::hide("m_subclustering_marker_7")
      shinyjs::show("m_subclustering_marker_8")
      shinyjs::hide("m_subclustering_marker_9")
      shinyjs::hide("m_subclustering_marker6")
      shinyjs::hide("m_subclustering_marker7")
      shinyjs::show("m_subclustering_marker8")
      shinyjs::hide("m_subclustering_marker9")
      shinyjs::hide("m_subclustering_marker_box6")
      shinyjs::show("m_subclustering_marker10")
      shinyjs::hide("m_subclustering_marker11")
      shinyjs::show("m_subclustering_marker12")
    }
    
    else if (input$m_subclustering_marker1 == 5) {
      shinyjs::hide("m_subclustering_marker_6")
      shinyjs::hide("m_subclustering_marker_7")
      shinyjs::show("m_subclustering_marker_8")
      shinyjs::show("m_subclustering_marker_9")
      shinyjs::hide("m_subclustering_marker6")
      shinyjs::hide("m_subclustering_marker7")
      shinyjs::show("m_subclustering_marker8")
      shinyjs::show("m_subclustering_marker9")
      shinyjs::hide("m_subclustering_marker_box6")
      shinyjs::show("m_subclustering_marker10")
      shinyjs::hide("m_subclustering_marker11")
      shinyjs::show("m_subclustering_marker12")
    }
    
  })
  
  ########multiple hide celltype box########  
  shinyjs::hide("m_subclustering_celltype_box3")
  shinyjs::hide("m_subclustering_celltype_box4")
  shinyjs::hide("m_subclustering_celltype_box5")
  shinyjs::hide("m_subclustering_celltype_box7")
  shinyjs::hide("m_subclustering_celltype_box8")
  shinyjs::hide("m_subclustering_celltype_box9")
  shinyjs::hide("m_subclustering_celltype_box10")
  
  observe({
    if (input$m_subclustering_celltype1 == 1) {
      shinyjs::show("m_subclustering_celltype_box2")
      shinyjs::show("m_subclustering_celltype2")
      shinyjs::hide("m_subclustering_celltype_box3")
      shinyjs::hide("m_subclustering_celltype3")
      shinyjs::hide("m_subclustering_celltype4")
      shinyjs::hide("m_subclustering_celltype_box4")
      shinyjs::hide("m_subclustering_celltype5")
      shinyjs::hide("m_subclustering_celltype6")
      shinyjs::hide("m_subclustering_celltype_box5")
      shinyjs::hide("m_subclustering_celltype7")
    }
    else if (input$m_subclustering_celltype1 == 2) {
      shinyjs::hide("m_subclustering_celltype_box2")
      shinyjs::hide("m_subclustering_celltype2")
      shinyjs::show("m_subclustering_celltype_box3")
      shinyjs::show("m_subclustering_celltype3")
      shinyjs::show("m_subclustering_celltype4")
      shinyjs::hide("m_subclustering_celltype_box4")
      shinyjs::hide("m_subclustering_celltype5")
      shinyjs::hide("m_subclustering_celltype6")
      shinyjs::hide("m_subclustering_celltype_box5")
      shinyjs::hide("m_subclustering_celltype7")
    }
    else if (input$m_subclustering_celltype1 == 3) {
      shinyjs::hide("m_subclustering_celltype_box2")
      shinyjs::hide("m_subclustering_celltype2")
      shinyjs::hide("m_subclustering_celltype_box3")
      shinyjs::hide("m_subclustering_celltype3")
      shinyjs::hide("m_subclustering_celltype4")
      shinyjs::show("m_subclustering_celltype_box4")
      shinyjs::show("m_subclustering_celltype5")
      shinyjs::show("m_subclustering_celltype6")
      shinyjs::hide("m_subclustering_celltype_box5")
      shinyjs::hide("m_subclustering_celltype7")
    }
    else if (input$m_subclustering_celltype1 == 4) {
      shinyjs::hide("m_subclustering_celltype_box2")
      shinyjs::hide("m_subclustering_celltype2")
      shinyjs::hide("m_subclustering_celltype_box3")
      shinyjs::hide("m_subclustering_celltype3")
      shinyjs::hide("m_subclustering_celltype4")
      shinyjs::hide("m_subclustering_celltype_box4")
      shinyjs::hide("m_subclustering_celltype5")
      shinyjs::hide("m_subclustering_celltype6")
      shinyjs::show("m_subclustering_celltype_box5")
      shinyjs::show("m_subclustering_celltype7")
    }
  })     
  
  ##################multiple Cluster-based plots####################### 
  shinyjs::hide("m_subclustering_clusterbased2")
  shinyjs::hide("m_subclustering_clusterbased_box2")
  shinyjs::hide("m_subclustering_clusterbased_box3")     
  shinyjs::hide("m_subclustering_clusterbased_box4")
  
  ##################multiple conditionbased####################### 
  shinyjs::hide("m_subclustering_conditionbased_box3")
  shinyjs::hide("m_subclustering_conditionbased_box4")     
  shinyjs::hide("m_subclustering_conditionbased_box5")      
  
  ###################################################
  ######################data Input##################
  
  #######################TAB2.1################################
  ########multiple samples subclustering hide stats########   
  
  observe({
    if (input$m_subclustering1 == "seurat_clusters") {
      shinyjs::show("m_subclustering_2")
      shinyjs::hide("m_subclustering_3")
      shinyjs::hide("m_subclustering_4")
      shinyjs::hide("m_subclustering_5")
      shinyjs::hide("m_subclustering_6")
      shinyjs::show("m_subclustering2")
      shinyjs::hide("m_subclustering3")
      shinyjs::hide("m_subclustering4")
      shinyjs::hide("m_subclustering5")
      shinyjs::hide("m_subclustering6")
    }  
    else if (input$m_subclustering1 == "predicted") {
      shinyjs::hide("m_subclustering_2")
      shinyjs::show("m_subclustering_3")
      shinyjs::hide("m_subclustering_4")
      shinyjs::hide("m_subclustering_5")
      shinyjs::hide("m_subclustering_6")
      shinyjs::hide("m_subclustering2")
      shinyjs::show("m_subclustering3")
      shinyjs::hide("m_subclustering4")
      shinyjs::hide("m_subclustering5")
      shinyjs::hide("m_subclustering6")
    }  
   
  })
  
  observeEvent(input$subclustering_multiple_sample_submit,{
    shinyjs::show("m_subclustering_box1")
    shinyjs::show("m_subclustering_box2")
    shinyjs::show("m_subclustering_box3")
  })
  
  output$m_subclustering_2 <- renderUI ({
    clusters <- req(datainput_multiple_celltype_level()[[2]])
    shinyWidgets::pickerInput(
      inputId = "m_subclustering2",
      label = "Select one or multiple cluster(s) for analsysis",
      choices = sort(clusters),
      selected = sort(clusters)[1],
      multiple = T,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_subclustering_3 <- renderUI ({
    clusters <- req(datainput_multiple_celltype_level()[[3]])
    shinyWidgets::pickerInput(
      inputId = "m_subclustering3",
      label = "Select one or multiple cluster(s) for analsysis",
      choices = sort(clusters),
      multiple = T,
      options = list(`actions-box` = TRUE))
  })
  
  
  datainput_subclustering_multiple_sample_level<- eventReactive(input$subclustering_multiple_sample_submit,{
    source("scripts/subclustering_multiple_stats.R")
   #datainput_subclustering_multiple_sample(index_subclustering_multiple_sample_file = datainput_multiple_celltype_level()[[1]], index_m_subclustering1 = input$m_subclustering1, index_m_subclustering2 = input$m_subclustering2, index_m_subclustering3 = input$m_subclustering3, index_m_subclustering4 = input$m_subclustering4, index_m_subclustering5 = input$m_subclustering5, index_m_subclustering6 = input$m_subclustering6)
    datainput_subclustering_multiple_sample(index_subclustering_multiple_sample_file = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_celltype = datainput_multiple_celltype_level()[[4]], index_m_subclustering1 = input$m_subclustering1, index_m_subclustering2 = input$m_subclustering2, index_m_subclustering3 = input$m_subclustering3)
  })
  
  output$m_subclustering_qc <- renderPlot({
    datainput_subclustering_multiple_sample_level()[1]
  })
  
  observeEvent(input$download_m_subclustering_qc, {
    showModal(modalDialog(
      title = strong("Download QC plot"),
      numericInput("m_subclustering_qc_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_qc_plot_width", label = h5("Figure width (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_subclustering_qc_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_qc_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_qc_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  
  output$m_subclustering_qc_downloadoutput<- downloadHandler(
    filename = function(){
      paste("QC_for_the_selected_subclusters", input$m_subclustering_qc_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_sample_level()[[1]], width = input$m_subclustering_qc_plot_width, height = input$m_subclustering_qc_plot_height, dpi = input$m_subclustering_qc_plot_dpi, units = "in")
    }
  )
  
  output$subclustering_multiple_cell_table<- renderDataTable(DT::datatable((datainput_subclustering_multiple_sample_level()[[2]]),
                                                                           options = list(
                                                                             scrollX = TRUE,
                                                                             pageLength = 10,
                                                                             bFilter=0
                                                                           ),rownames= FALSE, selection = "none"))
  
  output$download_subclustering_multiple_cell_table <- downloadHandler(
    filename = function() { 
      paste("Number of cells", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_sample_level()[[2]], file)
    }
  )
  
  
  ###################save seurat object before qc###################
  output$m_subclustering_stats<- downloadHandler(
    filename = function(){
      paste("multiple_sample_subclustering_seuart_object.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_sample_level()[[3]], file= file)
    }
  )
  ###############link to next tab###########################      
  observeEvent(input$link_m_subclustering_normalization, {
    newvalue <- "Normalization and PCA analysis"
    updateTabsetPanel(session, "subclustering_multiple_tabsets", newvalue)
  })       
  
  
  
  
  ##########################Tab2.2###############################      
  ##############multiple Normalization & PCA###################      
  datainput_subclustering_multiple_normalization_pca_level <- eventReactive(input$subclustering_multiple_sample_normalization,{
    source("scripts/subclustering_multiple_normalization_pca.R")
    datainput_subclustering_multiple_normalization_pca(index_subclustering_multiple_normalization_pca_input = datainput_subclustering_multiple_sample_level()[[3]], index_subclustering_multiple_sample_normalization_method = input$subclustering_multiple_sample_normalization_method, index_subclustering_multiple_sample_scale_factor=input$subclustering_multiple_sample_scale_factor, index_subclustering_multiple_sample_var_genes = input$subclustering_multiple_sample_var_genes, index_subclustering_multiple_sample_normalization_variable_genes=input$subclustering_multiple_sample_normalization_variable_genes, index_subclustering_multiple_sample_pca_dim=input$subclustering_multiple_sample_pca_dim)
  })
  
  
  output$m_subclustering_pca_plot<-renderPlot({
    datainput_subclustering_multiple_normalization_pca_level()[1]
  })
  
  observeEvent(input$download_m_subclustering_pca_plot, {
    showModal(modalDialog(
      title = strong("Download PCA Plot"),
      numericInput("m_subclustering_pca_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_subclustering_pca_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_pca_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_pca_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_pca_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_subclustering_pca_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_PCA_plot", input$m_subclustering_pca_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_normalization_pca_level()[[1]], width = input$m_subclustering_pca_plot_width, height = input$m_subclustering_pca_plot_height, dpi = input$m_subclustering_pca_plot_dpi, units = "in")
    }
  )
  
  output$m_subclustering_elbow_plot<-renderPlot({
    datainput_subclustering_multiple_normalization_pca_level()[2]
  })
  
  observeEvent(input$download_m_subclustering_elbow_plot, {
    showModal(modalDialog(
      title = strong("Download Variable Features Plot"),
      numericInput("m_subclustering_elbow_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_subclustering_elbow_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_elbow_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_elbow_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_elbow_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_subclustering_elbow_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_Elbow", input$m_subclustering_elbow_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_normalization_pca_level()[[2]], width = input$m_subclustering_elbow_plot_width, height = input$m_subclustering_elbow_plot_height, dpi = input$m_subclustering_elbow_plot_dpi, units = "in")
    }
  )
  
  output$m_subclustering_pca2_plot<-renderPlot({
    datainput_subclustering_multiple_normalization_pca_level()[3]
  })
  
  observeEvent(input$download_m_subclustering_pca2_plot, {
    showModal(modalDialog(
      title = strong("Download PCA Plot"),
      numericInput("m_subclustering_pca2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_subclustering_pca2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_pca2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_pca2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_pca2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_subclustering_pca2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_PCA_plot_sample_based", input$m_subclustering_pca2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_normalization_pca_level()[[3]], width = input$m_subclustering_pca2_plot_width, height = input$m_subclustering_pca2_plot_height, dpi = input$m_subclustering_pca2_plot_dpi, units = "in")
    }
  )
  
  
  output$m_subclustering_pca3_plot<-renderPlot({
    datainput_subclustering_multiple_normalization_pca_level()[4]
  })
  
  
  observeEvent(input$download_m_subclustering_pca3_plot, {
    showModal(modalDialog(
      title = strong("Download PCA Plot"),
      numericInput("m_subclustering_pca3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 5, width = "300px"),
      numericInput("m_subclustering_pca3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_pca3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_pca3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_pca3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_subclustering_pca3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("After_normalization_PCA_plot_group_based", input$m_subclustering_pca3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_normalization_pca_level()[[4]], width = input$m_subclustering_pca3_plot_width, height = input$m_subclustering_pca3_plot_height, dpi = input$m_subclustering_pca3_plot_dpi, units = "in")
    }
  )
  
  
  
  ###################save seurat object after normalization###################
  output$m_subclustering_normalization<- downloadHandler(
    filename = function(){
      paste("subclustering_multiple_sample_seuart_object_after_normalization.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_normalization_pca_level()[[5]], file= file)
    }
  )
  
  #####################################link to next tab###########################      
  observeEvent(input$link_m_subclustering_clustering, {
    newvalue <- "Clustering"
    updateTabsetPanel(session, "subclustering_multiple_tabsets", newvalue)
  })       
  
  
  
  #####################################################Tab2.3####################      
  ########################################multiple Clustering###################      
  datainput_subclustering_multiple_clustering_level <- eventReactive(input$subclustering_multiple_sample_clustering,{
    source("scripts/subclustering_multiple_clustering.R")
    datainput_subclustering_multiple_clustering(index_subclustering_multiple_clustering_input = datainput_subclustering_multiple_normalization_pca_level()[[5]], index_subclustering_multiple_sample_normalization_method = input$subclustering_multiple_sample_normalization_method, index_m_subclustering_clustering1 = input$m_subclustering_clustering1, index_m_subclustering_clustering2 = input$m_subclustering_clustering2, index_m_subclustering_clustering3 = input$m_subclustering_clustering3, index_m_subclustering_clustering4 = input$m_subclustering_clustering4, index_m_subclustering_clustering5 = input$m_subclustering_clustering5, index_m_subclustering_clustering6 = input$m_subclustering_clustering6, index_m_subclustering_clustering7 = input$m_subclustering_clustering7, index_m_subclustering_clustering8 = input$m_subclustering_clustering8, index_m_subclustering_clustering9 = input$m_subclustering_clustering9, index_m_subclustering_clustering10 = input$m_subclustering_clustering10, index_m_subclustering_clustering11 = input$m_subclustering_clustering11, index_m_subclustering_clustering12 = input$m_subclustering_clustering12, index_m_subclustering_clustering13 = input$m_subclustering_clustering13)
  })
  
  output$m_subclustering_umap_tsne1_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[1]
  })
  observeEvent(input$download_m_subclustering_umap_tsne1_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_umap_tsne1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Cluster_plot", input$m_subclustering_umap_tsne1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[1]], width = input$m_subclustering_umap_tsne1_plot_width, height = input$m_subclustering_umap_tsne1_plot_height, dpi = input$m_subclustering_umap_tsne1_plot_dpi, units = "in")
    }
  )
  
  output$m_subclustering_umap_tsne_bar1_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[2]
  }) 
  observeEvent(input$download_m_subclustering_umap_tsne_bar1_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne_bar1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne_bar1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne_bar1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne_bar1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne_bar1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_umap_tsne_bar1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Cluster_bar_plot", input$m_subclustering_umap_tsne_bar1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[2]], width = input$m_subclustering_umap_tsne_bar1_plot_width, height = input$m_subclustering_umap_tsne_bar1_plot_height, dpi = input$m_subclustering_umap_tsne_bar1_plot_dpi, units = "in")
    }
  )
  
  
  
  
  
  output$m_subclustering_umap_tsne2_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[3]
  })
  observeEvent(input$download_m_subclustering_umap_tsne2_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_umap_tsne2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Condition_based_plot", input$m_subclustering_umap_tsne2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[3]], width = input$m_subclustering_umap_tsne2_plot_width, height = input$m_subclustering_umap_tsne2_plot_height, dpi = input$m_subclustering_umap_tsne2_plot_dpi, units = "in")
    }
  )
  
    
  
  output$m_subclustering_umap_tsne_bar2_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[4]
  })
  observeEvent(input$download_m_subclustering_umap_tsne_bar2_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne_bar2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne_bar2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne_bar2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne_bar2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne_bar2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_umap_tsne_bar2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Condition_based_bar_plot", input$m_subclustering_umap_tsne_bar2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[4]], width = input$m_subclustering_umap_tsne_bar2_plot_width, height = input$m_subclustering_umap_tsne_bar2_plot_height, dpi = input$m_subclustering_umap_tsne_bar2_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_subclustering_umap_tsne3_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[5]
  })
  observeEvent(input$download_m_subclustering_umap_tsne3_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_umap_tsne3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Sample_based_plot", input$m_subclustering_umap_tsne3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[5]], width = input$m_subclustering_umap_tsne3_plot_width, height = input$m_subclustering_umap_tsne3_plot_height, dpi = input$m_subclustering_umap_tsne3_plot_dpi, units = "in")
    }
  )	
  
  
  
  output$m_subclustering_umap_tsne_bar3_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[6]
  })
  observeEvent(input$download_m_subclustering_umap_tsne_bar3_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne_bar3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne_bar3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_umap_tsne_bar3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne_bar3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne_bar3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_umap_tsne_bar3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Sample_based_bar_plot", input$m_subclustering_umap_tsne_bar3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[6]], width = input$m_subclustering_umap_tsne_bar3_plot_width, height = input$m_subclustering_umap_tsne_bar3_plot_height, dpi = input$m_subclustering_umap_tsne_bar3_plot_dpi, units = "in")
    }
  )
  
  output$m_subclustering_clustering_table1<- renderDataTable(DT::datatable((datainput_subclustering_multiple_clustering_level()[[7]]),
                                                                           options = list(
                                                                             scrollX = TRUE,
                                                                             pageLength = 10,
                                                                             bFilter=0
                                                                           ),rownames= FALSE, selection = "none"))
  
  output$download_m_subclustering_clustering_table1 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in clusters", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_clustering_level()[[7]], file)
    }
  ) 
  
  output$m_subclustering_clustering_table2<- renderDataTable(DT::datatable((datainput_subclustering_multiple_clustering_level()[[8]]),
                                                                           options = list(
                                                                             scrollX = TRUE,
                                                                             pageLength = 10,
                                                                             bFilter=0
                                                                           ),rownames= FALSE, selection = "none"))
  
  output$download_m_subclustering_clustering_table2 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in clusters based on condition", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_clustering_level()[[8]], file)
    }
  ) 
  
  output$m_subclustering_clustering_table3<- renderDataTable(DT::datatable((datainput_subclustering_multiple_clustering_level()[[9]]),
                                                                           options = list(
                                                                             scrollX = TRUE,
                                                                             pageLength = 10,
                                                                             bFilter=0
                                                                           ),rownames= FALSE, selection = "none"))
  
  output$download_m_subclustering_clustering_table3 <- downloadHandler(
    filename = function() { 
      paste("Number of cells in clusters based on samples", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_clustering_level()[[9]], file)
    }
  ) 
  
  output$m_subclustering_umap_tsne4_plot<-renderPlot({
    datainput_subclustering_multiple_clustering_level()[15]
  })
  
  observeEvent(input$download_m_subclustering_umap_tsne4_plot, {
    showModal(modalDialog(
      title = strong("Download UMAP/ t-SNE Plot"),
      numericInput("m_subclustering_umap_tsne4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_subclustering_umap_tsne4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 15, width = "300px"),
      numericInput("m_subclustering_umap_tsne4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_umap_tsne4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_umap_tsne4_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_subclustering_umap_tsne4_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Sample_based_plot_split_by_clusters", input$m_subclustering_umap_tsne4_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clustering_level()[[15]], width = input$m_subclustering_umap_tsne4_plot_width, height = input$m_subclustering_umap_tsne4_plot_height, dpi = input$m_subclustering_umap_tsne4_plot_dpi, units = "in")
    }
  )	  
  
  ###################save seurat object after clustering###################
  output$m_subclustering_clustering<- downloadHandler(
    filename = function(){
      paste("subclustering_multiple_sample_seuart_object_after_clustering.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_clustering_level()[[10]], file= file)
    }
  )
  
  
  #####################################link to next tab###########################   
  observeEvent(input$link_m_subclustering_marker, {
    newvalue <- "Markers Identification"
    updateTabsetPanel(session, "subclustering_multiple_tabsets", newvalue)
  }) 
  
  
  ##########################Tab2.4###############################      
  ##############multiple Marker identification###################    
  observeEvent(input$subclustering_multiple_sample_marker,{
    if(input$m_subclustering_marker1 == 1){
      shinyjs::show("m_subclustering_marker_box5")
      shinyjs::show("m_subclustering_marker_box6")
      shinyjs::show("m_subclustering_marker_box7")
    }
    else{
      shinyjs::show("m_subclustering_marker_box5")
      shinyjs::hide("m_subclustering_marker_box6")
      shinyjs::show("m_subclustering_marker_box7")
    }
    
  })
  
  output$m_subclustering_marker_6 <- renderUI ({
    clusters <- req(datainput_subclustering_multiple_clustering_level()[[11]])
    
    shinyWidgets::pickerInput(
      inputId = "m_subclustering_marker6",
      label = "Select one cluster for analsysis",
      choices = sort(clusters),
      multiple = F,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_subclustering_marker_7 <- renderUI ({
    clusters <- req(datainput_subclustering_multiple_clustering_level()[[11]])
    clusters <- clusters[!clusters == input$m_subclustering_marker6]
    shinyWidgets::pickerInput(
      inputId = "m_subclustering_marker7",
      label = "Identify markers distinguishing a cluster from other selected clusters",
      choices = sort(clusters),
      multiple = T,
      selected = sort(clusters)[1],
      options = list(`actions-box` = TRUE))
  })
  
  output$m_subclustering_marker_8 <- renderUI ({
    clusters <- req(datainput_subclustering_multiple_clustering_level()[[11]])
    
    shinyWidgets::pickerInput(
      inputId = "m_subclustering_marker8",
      label = "Select one cluster to define markers",
      choices = sort(clusters),
      multiple = F,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_subclustering_marker_9 <- renderUI ({
    clusters <- req(datainput_subclustering_multiple_clustering_level()[[11]])
    clusters <- clusters[!clusters == input$m_subclustering_marker8]
    shinyWidgets::pickerInput(
      inputId = "m_subclustering_marker9",
      label = "Select the cluster to find the conserved markers between two clusters",
      choices = sort(clusters),
      multiple = T,
      selected = sort(clusters)[1],
      options = list(`actions-box` = TRUE))
  })
  
  
  datainput_subclustering_multiple_marker_level <- eventReactive(input$subclustering_multiple_sample_marker,{
    
    source("scripts/subclustering_multiple_marker.R")
    datainput_subclustering_multiple_marker(index_subclustering_multiple_marker_input = datainput_subclustering_multiple_clustering_level()[[10]], index_m_subclustering_marker1 = input$m_subclustering_marker1, index_m_subclustering_marker2 = input$m_subclustering_marker2, index_m_subclustering_marker3 = input$m_subclustering_marker3, index_m_subclustering_marker4 = input$m_subclustering_marker4, index_m_subclustering_marker5 = input$m_subclustering_marker5, index_m_subclustering_marker6 = input$m_subclustering_marker6, index_m_subclustering_marker7 = input$m_subclustering_marker7, index_m_subclustering_marker8 = input$m_subclustering_marker8, index_m_subclustering_marker9 = input$m_subclustering_marker9, index_m_subclustering_marker10 = input$m_subclustering_marker10, index_subclustering_multiple_sample_normalization_method = input$subclustering_multiple_sample_normalization_method)
  })
  
  
  output$m_subclustering_marker1_table<- renderDataTable(DT::datatable((datainput_subclustering_multiple_marker_level()[[1]]),
                                                                       options = list(
                                                                         scrollX = TRUE,
                                                                         pageLength = 10,
                                                                         dom = "Blfrtip"
                                                                         #bFilter=0
                                                                       ),rownames= FALSE, selection = "none"))
  
  
  output$download_m_subclustering_marker1_table <- downloadHandler(
    filename = function() {
      paste("Number_of_identified_markers_or_differentially_expressed_genes", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_marker_level()[[1]], file)
    }
  )
  
  
  output$m_subclustering_marker1_plot<-renderPlot({
    datainput_subclustering_multiple_marker_level()[3]
  })
  
  
  observeEvent(input$download_m_subclustering_marker1_plot, {
    showModal(modalDialog(
      title = strong("Download Heatmap"),
      numericInput("m_subclustering_marker1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_marker1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 12, width = "300px"),
      numericInput("m_subclustering_marker1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_marker1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_marker1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  
  output$m_subclustering_marker1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste( "Heatmap_with_Top5_expressed_genes", input$m_subclustering_marker1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_marker_level()[[3]], width = input$m_subclustering_marker1_plot_width, height = input$m_subclustering_marker1_plot_height, dpi = input$m_subclustering_marker1_plot_dpi, units = "in")
    }
  )
  
  
  ###################save seurat object after doublet removal###################
  output$m_subclustering_marker<- downloadHandler(
    filename = function(){
      paste("subclustering_multiple_sample_seuart_object_after_marker_identification.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_marker_level()[[2]], file= file)
    }
  )
  
  #####################################link to next tab###########################     
  observeEvent(input$link_m_subclustering_prediction, {
    newvalue <- "Cell Type Prediction"
    updateTabsetPanel(session, "subclustering_multiple_tabsets", newvalue)
  })  
  
  ##########################Tab2.5###############################      
  ##############multiple cell type################### 
  output$m_subclustering_celltype7 <- renderUI({
    numberofclusters <- as.integer(length(levels(datainput_subclustering_multiple_marker_level()[[2]])))
    lapply(1:numberofclusters, function(i) {
      column(3, textInput(paste("subclustering_mcelltypenames", levels(datainput_subclustering_multiple_marker_level()[[2]])[i], sep = ""),
                          paste("Cluster", levels(datainput_subclustering_multiple_marker_level()[[2]])[i]), value = paste("Cluster", levels(datainput_subclustering_multiple_marker_level()[[2]])[i])))
    })
  })
  
  observeEvent(input$subclustering_multiple_sample_celltype,{
    if(input$m_subclustering_celltype1 == 1){
      shinyjs::show("m_subclustering_celltype_box7")
      shinyjs::show("m_subclustering_celltype_box8")
      shinyjs::hide("m_subclustering_celltype_box9")
      shinyjs::show("m_subclustering_celltype_box10")
      shinyjs::show("m_subclustering_celltype9")
      shinyjs::hide("m_subclustering_celltype10")
    }
    else if (input$m_subclustering_celltype1 == 2){
      shinyjs::show("m_subclustering_celltype_box7")
      shinyjs::show("m_subclustering_celltype_box8")
      shinyjs::show("m_subclustering_celltype_box9")
      shinyjs::show("m_subclustering_celltype_box10")
      shinyjs::hide("m_subclustering_celltype9")
      shinyjs::show("m_subclustering_celltype10")
    }
    else if (input$m_subclustering_celltype1 == 3){
      shinyjs::show("m_subclustering_celltype_box7")
      shinyjs::hide("m_subclustering_celltype_box8")
      shinyjs::hide("m_subclustering_celltype_box9")
      shinyjs::show("m_subclustering_celltype_box10")
      shinyjs::hide("m_subclustering_celltype9")
      shinyjs::hide("m_subclustering_celltype10")
    }
    else if (input$m_subclustering_celltype1 == 4){
      shinyjs::show("m_subclustering_celltype_box7")
      shinyjs::hide("m_subclustering_celltype_box8")
      shinyjs::hide("m_subclustering_celltype_box9")
      shinyjs::show("m_subclustering_celltype_box10")
      shinyjs::hide("m_subclustering_celltype9")
      shinyjs::hide("m_subclustering_celltype10")
    }
    
  })
  
  
  datainput_subclustering_multiple_celltype_level <- eventReactive(input$subclustering_multiple_sample_celltype,{
    source("scripts/subclustering_multiple_celltype.R")
    datainput_subclustering_multiple_celltype(index_subclustering_multiple_celltype_input = datainput_subclustering_multiple_marker_level()[[2]], index_cell_markers = datainput_subclustering_multiple_marker_level()[[1]], index_m_subclustering_celltype1 = input$m_subclustering_celltype1, index_m_subclustering_celltype2 = input$m_subclustering_celltype2, index_m_subclustering_celltype3 = input$m_subclustering_celltype3, index_m_subclustering_celltype4 = input$m_subclustering_celltype4, index_m_subclustering_celltype5 = input$m_subclustering_celltype5, index_m_subclustering_celltype6 = input$m_subclustering_celltype6, 
                                              index_m_subclustering_celltype7 = c(input$subclustering_mcelltypenames0,input$subclustering_mcelltypenames1,input$subclustering_mcelltypenames2,input$subclustering_mcelltypenames3,input$subclustering_mcelltypenames4,input$subclustering_mcelltypenames5,input$subclustering_mcelltypenames6,input$subclustering_mcelltypenames7,input$subclustering_mcelltypenames8,input$subclustering_mcelltypenames9,input$subclustering_mcelltypenames10,input$subclustering_mcelltypenames11,input$subclustering_mcelltypenames12,input$subclustering_mcelltypenames13,input$subclustering_mcelltypenames14,input$subclustering_mcelltypenames15,input$subclustering_mcelltypenames16,input$subclustering_mcelltypenames17,input$subclustering_mcelltypenames18,input$subclustering_mcelltypenames19,input$subclustering_mcelltypenames20,input$subclustering_mcelltypenames21,input$subclustering_mcelltypenames22,input$subclustering_mcelltypenames23,input$subclustering_mcelltypenames24,input$subclustering_mcelltypenames25,input$subclustering_mcelltypenames26,input$subclustering_mcelltypenames27,input$subclustering_mcelltypenames28,input$subclustering_mcelltypenames29,input$subclustering_mcelltypenames30,input$subclustering_mcelltypenames31,input$subclustering_mcelltypenames32,input$subclustering_mcelltypenames33,input$subclustering_mcelltypenames34,input$subclustering_mcelltypenames35,input$subclustering_mcelltypenames36,input$subclustering_mcelltypenames37,input$subclustering_mcelltypenames38,input$subclustering_mcelltypenames39,input$subclustering_mcelltypenames40,input$subclustering_mcelltypenames41,input$subclustering_mcelltypenames42,input$subclustering_mcelltypenames43,input$subclustering_mcelltypenames44,input$subclustering_mcelltypenames45,input$subclustering_mcelltypenames46,input$subclustering_mcelltypenames47,input$subclustering_mcelltypenames48,input$subclustering_mcelltypenames49,input$subclustering_mcelltypenames50,input$subclustering_mcelltypenames51,input$subclustering_mcelltypenames52,input$subclustering_mcelltypenames53,input$subclustering_mcelltypenames54,input$subclustering_mcelltypenames55,input$subclustering_mcelltypenames56,input$subclustering_mcelltypenames57,input$subclustering_mcelltypenames58,input$subclustering_mcelltypenames59,input$subclustering_mcelltypenames60,input$subclustering_mcelltypenames61,input$subclustering_mcelltypenames62,input$subclustering_mcelltypenames63,input$subclustering_mcelltypenames64,input$subclustering_mcelltypenames65,input$subclustering_mcelltypenames66,input$subclustering_mcelltypenames67,input$subclustering_mcelltypenames68,input$subclustering_mcelltypenames69,input$subclustering_mcelltypenames70,input$subclustering_mcelltypenames71,input$subclustering_mcelltypenames72,input$subclustering_mcelltypenames73,input$subclustering_mcelltypenames74,input$subclustering_mcelltypenames75,input$subclustering_mcelltypenames76,input$subclustering_mcelltypenames77,input$subclustering_mcelltypenames78,input$subclustering_mcelltypenames79,input$subclustering_mcelltypenames80,input$subclustering_mcelltypenames81,input$subclustering_mcelltypenames82,input$subclustering_mcelltypenames83,input$subclustering_mcelltypenames84,input$subclustering_mcelltypenames85,input$subclustering_mcelltypenames86,input$subclustering_mcelltypenames87,input$subclustering_mcelltypenames88,input$subclustering_mcelltypenames89,input$subclustering_mcelltypenames90,input$subclustering_mcelltypenames91,input$subclustering_mcelltypenames92,input$subclustering_mcelltypenames93,input$subclustering_mcelltypenames94,input$subclustering_mcelltypenames95,input$subclustering_mcelltypenames96,input$subclustering_mcelltypenames97,input$subclustering_mcelltypenames98,input$subclustering_mcelltypenames99), 
                                              index_m_subclustering_celltype8 = input$m_subclustering_celltype8, index_m_subclustering_clustering6 = input$m_subclustering_clustering6, index_subclustering_multiple_sample_normalization_method = input$subclustering_multiple_sample_normalization_method)
  })
  output$m_subclustering_celltype1_plot<-renderPlot({
    datainput_subclustering_multiple_celltype_level()[5]
  })
  
  observeEvent(input$download_m_subclustering_celltype1_plot, {
    showModal(modalDialog(
      title = strong("Download Celltype"),
      numericInput("m_subclustering_celltype1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_celltype1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 16, width = "300px"),
      numericInput("m_subclustering_celltype1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_celltype1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_celltype1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_celltype1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Dimplot_with_celltype", input$m_subclustering_celltype1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_celltype_level()[[5]], width = input$m_subclustering_celltype1_plot_width, height = input$m_subclustering_celltype1_plot_height, dpi = input$m_subclustering_celltype1_plot_dpi, units = "in")
    }
  )
  
  output$m_subclustering_celltype1_table<- renderDataTable(DT::datatable((datainput_subclustering_multiple_celltype_level()[[6]]),
                                                                         options = list(
                                                                           scrollX = TRUE,
                                                                           pageLength = 10,
                                                                           dom = "Blfrtip"
                                                                           #bFilter=0
                                                                         ),rownames= FALSE, selection = "none"))
  
  output$download_m_subclustering_celltype1_table <- downloadHandler(
    filename = function() { 
      paste("predicted_celltype_Scores", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_celltype_level()[[6]], file)
    }
  )
  
  output$m_subclustering_celltype2_plot<-renderPlot({
    datainput_subclustering_multiple_celltype_level()[7]
  })
  
  observeEvent(input$download_m_subclustering_celltype2_plot, {
    showModal(modalDialog(
      title = strong("Download Heatmap"),
      numericInput("m_subclustering_celltype2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_celltype2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_celltype2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_celltype2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_celltype2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_celltype2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("score_plots", input$m_subclustering_celltype2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_celltype_level()[[7]], width = input$m_subclustering_celltype2_plot_width, height = input$m_subclustering_celltype2_plot_height, dpi = input$m_subclustering_celltype2_plot_dpi, units = "in")
    }
  )
  
  
  output$m_subclustering_celltype3_plot<-renderPlot({
    datainput_subclustering_multiple_celltype_level()[8]
  })
  
  observeEvent(input$download_m_subclustering_celltype3_plot, {
    showModal(modalDialog(
      title = strong("Download Heatmap"),
      numericInput("m_subclustering_celltype3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_celltype3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("m_subclustering_celltype3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_celltype3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_celltype3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_celltype3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("score_plots", input$m_subclustering_celltype3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_celltype_level()[[8]], width = input$m_subclustering_celltype3_plot_width, height = input$m_subclustering_celltype3_plot_height, dpi = input$m_subclustering_celltype3_plot_dpi, units = "in")
    }
  )
  
  
  
  ###################save seurat object after doublet###################
  output$m_subclustering_celltype<- downloadHandler(
    filename = function(){
      paste("subclustering_multiple_sample_seuart_object_after_celltypes.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_celltype_level()[[1]], file= file)
    }
  )
  #####################################link to next tab###########################     
  observeEvent(input$link_m_subclustering_clusterbased, {
    newvalue <- "Cluster-based plots"
    updateTabsetPanel(session, "subclustering_multiple_tabsets", newvalue)
  })  
  
  
  ###################################Tab2.6###############################  
  ###########################Cluster-based plots####################### 
  observe({
    if(input$m_subclustering_clusterbased1 == "gene_name_list"){
      shinyjs::show("m_subclustering_clusterbased2")
    }
    else {
      shinyjs::hide("m_subclustering_clusterbased2")
    }
  })  
  
  observeEvent(input$subclustering_multiple_sample_clusterbased, {
    if(input$m_subclustering_clusterbased4 ==  "seurat_clusters"){
      shinyjs::show("m_subclustering_clusterbased_box3")
    }
    else {
      shinyjs::hide("m_subclustering_clusterbased_box3")
    }
  })  
  
  observe({
    if(input$m_subclustering_clusterbased3 ==  "RidgePlot"){
      shinyjs::hide("m_subclustering_clusterbased5")
    }
    else {
      shinyjs::show("m_subclustering_clusterbased5")
    }
  })
  
  observeEvent(input$subclustering_multiple_sample_clusterbased,{
    shinyjs::show("m_subclustering_clusterbased_box2")
    #shinyjs::show("m_subclustering_clusterbased_box3")
    shinyjs::show("m_subclustering_clusterbased_box4")
  })
  
  
  
  datainput_subclustering_multiple_clusterbased_level <- eventReactive(input$subclustering_multiple_sample_clusterbased,{
    source("scripts/subclustering_multiple_clusterbased.R")
    datainput_subclustering_multiple_clusterbased(index_subclustering_multiple_clusterbased_input = datainput_subclustering_multiple_celltype_level()[[1]], index_subclustering_multiple_clusterbased_features = datainput_subclustering_multiple_marker_level()[[1]], index_m_subclustering_celltype_method = datainput_subclustering_multiple_celltype_level()[[4]], index_m_subclustering_clusterbased1 = input$m_subclustering_clusterbased1, index_m_subclustering_clusterbased2 = input$m_subclustering_clusterbased2, index_m_subclustering_clusterbased3 = input$m_subclustering_clusterbased3, index_m_subclustering_clusterbased4 = input$m_subclustering_clusterbased4, index_m_subclustering_clusterbased5 = input$m_subclustering_clusterbased5)
  })  
  
  output$m_subclustering_clusterbased1_plot<-renderPlot({
    datainput_subclustering_multiple_clusterbased_level()[1]
  })
  
  observeEvent(input$download_m_subclustering_clusterbased1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("m_subclustering_clusterbased1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_subclustering_clusterbased1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_subclustering_clusterbased1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_clusterbased1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_clusterbased1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_clusterbased1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Plots_for_top_or_selected_markers",  input$m_subclustering_clusterbased1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_clusterbased_level()[[1]], width = input$m_subclustering_clusterbased1_plot_width, height = input$m_subclustering_clusterbased1_plot_height, dpi = input$m_subclustering_clusterbased1_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_subclustering_clusterbased1_table<- renderDataTable(DT::datatable((datainput_subclustering_multiple_clusterbased_level()[[3]]),
                                                                             options = list(
                                                                               scrollX = TRUE,
                                                                               pageLength = 10,
                                                                               dom = "Blfrtip"
                                                                               #bFilter=0
                                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_m_subclustering_clusterbased1_table <- downloadHandler(
    filename = function() { 
      paste("Top_or_selected_Cell_counts_proportion", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_clusterbased_level()[[3]], file)
    }
  )
  
  
  ###################save seurat object after doublet###################
  output$m_subclustering_clusterbased <- downloadHandler(
    filename = function(){
      paste("subclustering_multiple_sample_seuart_object_after_plots.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_clusterbased_level()[[2]], file= file)
    }
  )
  
  #####################################link to next tab###########################     
  observeEvent(input$link_m_subclustering_conditionbased, {
    newvalue <- "Condition based analysis"
    updateTabsetPanel(session, "subclustering_multiple_tabsets", newvalue)
  })  
  
  
  ##########################Tab2.7###############################      
  ##############Condition-based analysis###################    
  observe({
    if(input$m_subclustering_conditionbased7 == "VolcanoPlot"){
      shinyjs::hide("m_subclustering_conditionbased8")
      shinyjs::hide("m_subclustering_conditionbased9")
      #shinyjs::hide("m_subclustering_conditionbased_box4")
    }
    else{
      shinyjs::show("m_subclustering_conditionbased8")
      shinyjs::show("m_subclustering_conditionbased9")
      #shinyjs::show("m_subclustering_conditionbased_box4")
    }
  })
  
  observeEvent(input$subclustering_multiple_sample_conditionbased,{
    shinyjs::show("m_subclustering_conditionbased_box3")
    shinyjs::show("m_subclustering_conditionbased_box4")
    shinyjs::show("m_subclustering_conditionbased_box5")
  })
  
  
  
  output$m_subclustering_conditionbased_1 <- renderUI ({
    clusters <- req(datainput_subclustering_multiple_clustering_level()[[13]])
    
    shinyWidgets::pickerInput(
      inputId = "m_subclustering_conditionbased1",
      label = "Select the Condition1",
      choices = sort(clusters),
      multiple = F,
      options = list(`actions-box` = TRUE))
  })
  
  output$m_subclustering_conditionbased_2 <- renderUI ({
    clusters <- req(datainput_subclustering_multiple_clustering_level()[[13]])
    clusters <- clusters[!clusters == input$m_subclustering_conditionbased1]
    shinyWidgets::pickerInput(
      inputId = "m_subclustering_conditionbased2",
      label = "Select the Condition2",
      choices = sort(clusters),
      selected = sort(clusters)[1],
      multiple = F,
      options = list(`actions-box` = TRUE))
  })  
  
  
  
  datainput_subclustering_multiple_conditionbased_level <- eventReactive(input$subclustering_multiple_sample_conditionbased,{
    source("scripts/subclustering_multiple_conditionbased.R")
    datainput_subclustering_multiple_conditionbased(index_subclustering_multiple_conditionbased_input = datainput_subclustering_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_normalization_method = input$subclustering_multiple_sample_normalization_method, index_m_subclustering_conditionbased1 = input$m_subclustering_conditionbased1, index_m_subclustering_conditionbased2 = input$m_subclustering_conditionbased2, index_m_subclustering_conditionbased3 = input$m_subclustering_conditionbased3, index_m_subclustering_conditionbased4 = input$m_subclustering_conditionbased4, index_m_subclustering_conditionbased5 = input$m_subclustering_conditionbased5, index_m_subclustering_conditionbased6 = input$m_subclustering_conditionbased6, index_m_subclustering_conditionbased7 = input$m_subclustering_conditionbased7, index_m_subclustering_conditionbased8 = input$m_subclustering_conditionbased8, index_m_subclustering_conditionbased9 = input$m_subclustering_conditionbased9)
  })  
  
  output$m_subclustering_conditionbased1_plot<-renderPlot({
    datainput_subclustering_multiple_conditionbased_level()[1]
  })
  
  observeEvent(input$download_m_subclustering_conditionbased1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("m_subclustering_conditionbased1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_subclustering_conditionbased1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 20, width = "300px"),
      numericInput("m_subclustering_conditionbased1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("m_subclustering_conditionbased1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("m_subclustering_conditionbased1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$m_subclustering_conditionbased1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Plots_for_top_selected_markers",  input$m_subclustering_conditionbased1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_subclustering_multiple_conditionbased_level()[[1]], width = input$m_subclustering_conditionbased1_plot_width, height = input$m_subclustering_conditionbased1_plot_height, dpi = input$m_subclustering_conditionbased1_plot_dpi, units = "in")
    }
  )
  
  
  
  output$m_subclustering_conditionbased1_table<- renderDataTable(DT::datatable((datainput_subclustering_multiple_conditionbased_level()[[2]]),
                                                                               options = list(
                                                                                 scrollX = TRUE,
                                                                                 pageLength = 10,
                                                                                 dom = "Blfrtip"
                                                                                 #bFilter=0
                                                                               ),rownames= FALSE, selection = "none"))
  
  output$download_m_subclustering_conditionbased1_table <- downloadHandler(
    filename = function() { 
      paste("Top_or_selected_Cell_counts_proportion", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_subclustering_multiple_conditionbased_level()[[2]], file)
    }
  )
  
  
  ###################save seurat object after doublet###################
  output$m_subclustering_conditionbased <- downloadHandler(
    filename = function(){
      paste("subclustering_multiple_sample_seuart_object_after_plots.RDS")
    },
    content = function(file){
      saveRDS(datainput_subclustering_multiple_conditionbased_level()[[3]], file= file)
    }
  )
  
  
##############################################Menu3#####################################################################   
#####################################Cell Cluster Correlation Network##########################################################
###################hidebox#################
  shinyjs::hide("s_cccn_box1")
  shinyjs::hide("s_cccn_box2")
  
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_cccn_box0")
      shinyjs::show("s_cccn_box1")
    }
    else{
      shinyjs::show("s_cccn_box0")
      shinyjs::hide("s_cccn_box1") 
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_cccn_box0")
      shinyjs::show("s_cccn_box1")
    }
    else{
      shinyjs::show("s_cccn_box0")
      shinyjs::hide("s_cccn_box1") 
    }
  })
  
  observeEvent(input$single_multiple_sample_cccn,{
    shinyjs::show("s_cccn_box2")
  })
  
  
  datainput_single_multiple_sample_cccn_level <- eventReactive(input$single_multiple_sample_cccn,{
    source("scripts/cccn.R")
    datainput_single_multiple_sample_cccn(index_multiple_sample_cccn_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_cccn_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_cccn_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_cccn_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_multiple_sample_normalization_method_cccn = input$multiple_sample_normalization_method, index_subclustering_multiple_sample_normalization_method_cccn = input$subclustering_multiple_sample_normalization_method, index_s_cccn1 = input$s_cccn1, index_s_cccn2 = input$s_cccn2, index_s_cccn3 = input$s_cccn3)
  })  
  
  output$s_cccn1_plot<-renderPlot({
    datainput_single_multiple_sample_cccn_level()[1]
  })
  
  # observeEvent(input$download_s_cccn1_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cccn1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cccn1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cccn1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cccn1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cccn1_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cccn1_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Cluster_based_correlation_matrix_plot",  input$s_cccn1_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cccn_level()[[1]], width = input$s_cccn1_plot_width, height = input$s_cccn1_plot_height, dpi = input$s_cccn1_plot_dpi, units = "in")
  #   }
  # )
  
  
  output$s_cccn2_plot<-renderPlot({
    datainput_single_multiple_sample_cccn_level()[2]
  })
  observeEvent(input$download_s_cccn2_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_cccn2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cccn2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cccn2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_cccn2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_cccn2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_cccn2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Cluster_based_correlation_network_plot",  input$s_cccn2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_cccn_level()[[2]], width = input$s_cccn2_plot_width, height = input$s_cccn2_plot_height, dpi = input$s_cccn2_plot_dpi, units = "in")
    }
  )
  
  output$s_cccn1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_cccn_level()[[3]]),
                                                       options = list(
                                                         scrollX = TRUE,
                                                         pageLength = 10,
                                                         bFilter=0
                                                       ),rownames= TRUE, selection = "none"))
  
  output$download_s_cccn1_table <- downloadHandler(
    filename = function() { 
      paste("Cluster_based_correlation_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_cccn_level()[[3]], file)
    }
  )
  
  
  
  
 
######################################################Menu4#####################################################################   
####################################################GO terms####################################################################  
###################hidebox#################
  shinyjs::hide("s_go_box1")
  shinyjs::hide("s_go_box2")
  shinyjs::hide("s_go_box3")

 
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_go_box0")
      shinyjs::show("s_go_box1")
      shinyjs::show("s_go_box2")
    }
    else{
      shinyjs::show("s_go_box0")
      shinyjs::hide("s_go_box1") 
      shinyjs::hide("s_go_box2") 
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_go_box0")
      shinyjs::show("s_go_box1")
      shinyjs::show("s_go_box2")
    }
    else{
      shinyjs::show("s_go_box0")
      shinyjs::hide("s_go_box1") 
      shinyjs::hide("s_go_box2") 
    }
  })
  
  observe({
      if(input$s_go1 == "gene_name_list"){
        shinyjs::show("s_go14")
        shinyjs::hide("s_go2")
        shinyjs::hide("s_go3")
        shinyjs::hide("s_go_3")
        shinyjs::hide("s_go4")
      }
      else {
        shinyjs::hide("s_go14")
        shinyjs::show("s_go2")
        shinyjs::show("s_go3")
        shinyjs::show("s_go_3")
        shinyjs::show("s_go4")
      }
    })
    
    
  observeEvent(input$single_multiple_sample_go,{
    shinyjs::show("s_go_box3")
    })
  
  observe({
   
    if (input$s_go1 == "multiple_sample" & input$s_go2 == "seurat_clusters"){
      output$s_go_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_go3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_go1 == "multiple_sample" & input$s_go2 == "predicted"){
      output$s_go_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_go3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }  
    else if (input$s_go1 == "multiple_sample_subclustering" & input$s_go2 == "seurat_clusters"){
      output$s_go_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_go3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_go1 == "multiple_sample_subclustering" & input$s_go2 == "predicted"){
      output$s_go_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_go3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
  })
  
  datainput_single_multiple_sample_go_level <- eventReactive(input$single_multiple_sample_go,{
    source("scripts/go.R")
    datainput_single_multiple_sample_go(index_multiple_sample_go_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_go_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_go_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_go_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_multiple_sample_go_input3 = datainput_multiple_marker_level()[[1]], index_subclustering_multiple_sample_go_input3 = datainput_subclustering_multiple_marker_level()[[1]], index_s_go1 = input$s_go1, index_s_go2 = input$s_go2, index_s_go3 = input$s_go3, index_s_go4 = input$s_go4, index_s_go5 = input$s_go5, index_s_go6 = input$s_go6, index_s_go7 = input$s_go7, index_s_go8 = input$s_go8, index_s_go9 = input$s_go9, index_s_go10 = input$s_go10, index_s_go11 = input$s_go11,index_s_go12= input$s_go12, index_s_go13 = input$s_go13, index_s_go14 = input$s_go14)
  })  
  
  output$s_go1_plot<-renderPlot({
    datainput_single_multiple_sample_go_level()[1]
  })

  observeEvent(input$download_s_go1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_go1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_go1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_go1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_go1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_go1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })

  output$s_go1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Go_terms_", input$s_go12, input$s_go1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_go_level()[[1]], width = input$s_go1_plot_width, height = input$s_go1_plot_height, dpi = input$s_go1_plot_dpi, units = "in")
    }
  )
  
  output$s_go1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_go_level()[[2]]),
                                                             options = list(
                                                               scrollX = TRUE,
                                                               pageLength = 10,
                                                               dom = "Blfrtip"
                                                               #bFilter=0
                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_s_go1_table <- downloadHandler(
    filename = function() { 
      paste("Go_terms_summary_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_go_level()[[2]], file)
    }
  )
  
#############################################################Menu5#####################################################################   
####################################################Pathway analsis####################################################################  
###################hidebox#################
  shinyjs::hide("s_pathway_box1")
  shinyjs::hide("s_pathway_box2")
  shinyjs::hide("s_pathway_box3")
  
  
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_pathway_box0")
      shinyjs::show("s_pathway_box1")
      shinyjs::show("s_pathway_box2")
    }
    else{
      shinyjs::show("s_pathway_box0")
      shinyjs::hide("s_pathway_box1") 
      shinyjs::hide("s_pathway_box2")
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_pathway_box0")
      shinyjs::show("s_pathway_box1")
      shinyjs::show("s_pathway_box2")
    }
    else{
      shinyjs::show("s_pathway_box0")
      shinyjs::hide("s_pathway_box1") 
      shinyjs::hide("s_pathway_box2")
    }
  })
  
  observe({
    if(input$s_pathway1 == "gene_name_list"){
      shinyjs::show("s_pathway14")
      shinyjs::hide("s_pathway2")
      shinyjs::hide("s_pathway3")
      shinyjs::hide("s_pathway_3")
      shinyjs::hide("s_pathway4")
    }
    else {
      shinyjs::hide("s_pathway14")
      shinyjs::show("s_pathway2")
      shinyjs::show("s_pathway3")
      shinyjs::show("s_pathway_3")
      shinyjs::show("s_pathway4")
    }
  })
  
  observeEvent(input$single_multiple_sample_pathway,{
    shinyjs::show("s_pathway_box3")
  })
  
  observe({
   
    if (input$s_pathway1 == "multiple_sample" & input$s_pathway2 == "seurat_clusters"){
      output$s_pathway_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_pathway3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_pathway1 == "multiple_sample" & input$s_pathway2 == "predicted"){
      output$s_pathway_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_pathway3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }  
    else if (input$s_pathway1 == "multiple_sample_subclustering" & input$s_pathway2 == "seurat_clusters"){
      output$s_pathway_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_pathway3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_pathway1 == "multiple_sample_subclustering" & input$s_pathway2 == "predicted"){
      output$s_pathway_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_pathway3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
  })
  
  datainput_single_multiple_sample_pathway_level <- eventReactive(input$single_multiple_sample_pathway,{
    source("scripts/pathway.R")
    datainput_single_multiple_sample_pathway(index_multiple_sample_pathway_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_pathway_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_pathway_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_pathway_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_multiple_sample_pathway_input3 = datainput_multiple_marker_level()[[1]], index_subclustering_multiple_sample_pathway_input3 = datainput_subclustering_multiple_marker_level()[[1]], index_s_pathway1 = input$s_pathway1, index_s_pathway2 = input$s_pathway2, index_s_pathway3 = input$s_pathway3, index_s_pathway4 = input$s_pathway4, index_s_pathway5 = input$s_pathway5, index_s_pathway6 = input$s_pathway6, index_s_pathway7 = input$s_pathway7, index_s_pathway8 = input$s_pathway8, index_s_pathway9 = input$s_pathway9, index_s_pathway10 = input$s_pathway10, index_s_pathway11 = input$s_pathway11,index_s_pathway12= input$s_pathway12, index_s_pathway13 = input$s_pathway13, index_s_pathway14 = input$s_pathway14)
  })  
  
  output$s_pathway1_plot<-renderPlot({
    datainput_single_multiple_sample_pathway_level()[1]
  })
  
  observeEvent(input$download_s_pathway1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_pathway1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_pathway1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_pathway1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_pathway1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_pathway1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  output$s_pathway1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("pathway_", input$s_pathway6, "_", input$s_pathway12, input$s_pathway1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_pathway_level()[[1]], width = input$s_pathway1_plot_width, height = input$s_pathway1_plot_height, dpi = input$s_pathway1_plot_dpi, units = "in")
    }
  )
  
  output$s_pathway1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_pathway_level()[[2]]),
                                                     options = list(
                                                       scrollX = TRUE,
                                                       pageLength = 10,
                                                       dom = "Blfrtip"
                                                       #bFilter=0
                                                     ),rownames= FALSE, selection = "none"))
  
  output$download_s_pathway1_table <- downloadHandler(
    filename = function() { 
      paste("Pathway_summary_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_pathway_level()[[2]], file)
    }
  )
  
######################################################Menu6#####################################################################   
####################################################GSEA terms####################################################################  
###################hidebox#################
  shinyjs::hide("s_gsea_box1")
  shinyjs::hide("s_gsea_box2")
  shinyjs::hide("s_gsea_box3")
  
 
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_gsea_box0")
      shinyjs::show("s_gsea_box1")
      shinyjs::show("s_gsea_box2")
    }
    else{
      shinyjs::show("s_gsea_box0")
      shinyjs::hide("s_gsea_box1")
      shinyjs::hide("s_gsea_box2")
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_gsea_box0")
      shinyjs::show("s_gsea_box1")
      shinyjs::show("s_gsea_box2")
    }
    else{
      shinyjs::show("s_gsea_box0")
      shinyjs::hide("s_gsea_box1")
      shinyjs::hide("s_gsea_box2")
    }
  })
  
  
  
  observeEvent(input$single_multiple_sample_gsea,{
    shinyjs::show("s_gsea_box3")
  })
  
  observe({
    
    if (input$s_gsea1 == "multiple_sample" & input$s_gsea2 == "seurat_clusters"){
      output$s_gsea_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_gsea3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_gsea1 == "multiple_sample" & input$s_gsea2 == "predicted"){
      output$s_gsea_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_gsea3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }  
    else if (input$s_gsea1 == "multiple_sample_subclustering" & input$s_gsea2 == "seurat_clusters"){
      output$s_gsea_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_gsea3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_gsea1 == "multiple_sample_subclustering" & input$s_gsea2 == "predicted"){
      output$s_gsea_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_gsea3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = T,
          options = list(`actions-box` = TRUE))
      })
    }
  })
  
  datainput_single_multiple_sample_gsea_level <- eventReactive(input$single_multiple_sample_gsea,{
    source("scripts/gsea.R")
    datainput_single_multiple_sample_gsea(index_multiple_sample_gsea_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_gsea_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_gsea_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_gsea_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_multiple_sample_gsea_input3 = datainput_multiple_marker_level()[[1]], index_subclustering_multiple_sample_gsea_input3 = datainput_subclustering_multiple_marker_level()[[1]], index_s_gsea1 = input$s_gsea1, index_s_gsea2 = input$s_gsea2, index_s_gsea3 = input$s_gsea3, index_s_gsea4 = input$s_gsea4, index_s_gsea5 = input$s_gsea5, index_s_gsea6 = input$s_gsea6, index_s_gsea7 = input$s_gsea7, index_s_gsea8 = input$s_gsea8, index_s_gsea9 = input$s_gsea9, index_s_gsea10 = input$s_gsea10, index_s_gsea11 = input$s_gsea11,index_s_gsea12= input$s_gsea12)
  })  
  
  output$s_gsea1_plot<-renderPlot({
    datainput_single_multiple_sample_gsea_level()[1]
  })
  
  observeEvent(input$download_s_gsea1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_gsea1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 10, width = "300px"),
      numericInput("s_gsea1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 15, width = "300px"),
      numericInput("s_gsea1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_gsea1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_gsea1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  output$s_gsea1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("gsea_plot",  input$s_gsea1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_gsea_level()[[1]], width = input$s_gsea1_plot_width, height = input$s_gsea1_plot_height, dpi = input$s_gsea1_plot_dpi, units = "in")
    }
  )
  
  output$s_gsea1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_gsea_level()[[2]]),
                                                     options = list(
                                                       scrollX = TRUE,
                                                       pageLength = 10,
                                                       dom = "Blfrtip"
                                                       #bFilter=0
                                                     ),rownames= FALSE, selection = "none"))
  
  output$download_s_gsea1_table <- downloadHandler(
    filename = function() { 
      paste("GSEA_summary_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_gsea_level()[[2]], file)
    }
  )
 
  
##############################################Menu7#####################################################################   
#####################################Cell-Cell Communication##########################################################
  ###################hidebox#################
  shinyjs::hide("s_cellchat_box1")
  shinyjs::hide("s_cellchat_box2")
  shinyjs::hide("s_cellchat_box3")
  shinyjs::hide("s_cellchat_box4")
  shinyjs::hide("s_cellchat_box5")
  shinyjs::hide("s_cellchat_box6")
  shinyjs::hide("s_cellchat_box7")
  shinyjs::hide("s_cellchat_box8")
  shinyjs::hide("s_cellchat_box9")
  shinyjs::hide("s_cellchat_box10")
  shinyjs::hide("s_cellchat_box11")
  
  
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_cellchat_box0")
      shinyjs::show("s_cellchat_box1")
    }
    else{
      shinyjs::show("s_cellchat_box0")
      shinyjs::hide("s_cellchat_box1") 
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_cellchat_box0")
      shinyjs::show("s_cellchat_box1")
    }
    else{
      shinyjs::show("s_cellchat_box0")
      shinyjs::hide("s_cellchat_box1") 
    }
  })
  
  observeEvent(input$single_multiple_sample_cellchat1,{
    shinyjs::show("s_cellchat_box2")
    shinyjs::show("s_cellchat_box3")
    shinyjs::show("s_cellchat_box4")
    shinyjs::show("s_cellchat_box5")
    shinyjs::show("s_cellchat_box6")
  })
  
  observeEvent(input$single_multiple_sample_cellchat2,{
    shinyjs::show("s_cellchat_box7")
    shinyjs::show("s_cellchat_box8")
    shinyjs::show("s_cellchat_box9")
    shinyjs::show("s_cellchat_box10")
    shinyjs::show("s_cellchat_box11")
  })
  
  datainput_single_multiple_sample_cellchat1_level <- eventReactive(input$single_multiple_sample_cellchat1,{
    source("scripts/cellchat1.R")
    datainput_single_multiple_sample_cellchat1(index_multiple_sample_cellchat1_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_cellchat1_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_cellchat1_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_cellchat1_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_s_cellchat1 = input$s_cellchat1, index_s_cellchat2 = input$s_cellchat2, index_s_cellchat3 = input$s_cellchat3, index_s_cellchat4 = input$s_cellchat4, index_s_cellchat5 = input$s_cellchat5, index_s_cellchat6 = input$s_cellchat6, index_s_cellchat7 = input$s_cellchat7, index_s_cellchat8 = input$s_cellchat8, index_s_cellchat9 = input$s_cellchat9, index_s_cellchat10 = input$s_cellchat10)
  })  
  
  output$s_cellchat1_plot<-renderPlot({
    #grid.newpage() 
    datainput_single_multiple_sample_cellchat1_level()[1]
  })
  # observeEvent(input$download_s_cellchat1_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat1_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat1_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Number_of_interactions", input$s_cellchat1_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat1_level()[[1]], width = input$s_cellchat1_plot_width, height = input$s_cellchat1_plot_height, dpi = input$s_cellchat1_plot_dpi, units = "in")
  #   }
  # )
  
  
  
  output$s_cellchat2_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat1_level()[2]
  })
  # observeEvent(input$download_s_cellchat2_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat2_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat2_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Interaction_weights_or_strength", input$s_cellchat2_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat1_level()[[2]], width = input$s_cellchat2_plot_width, height = input$s_cellchat2_plot_height, dpi = input$s_cellchat2_plot_dpi, units = "in")
  #   }
  # )
  
  
  output$s_cellchat3_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat1_level()[3]
  })
  # observeEvent(input$download_s_cellchat3_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat3_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat3_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Interactions_heatmap", input$s_cellchat3_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat1_level()[[3]], width = input$s_cellchat3_plot_width, height = input$s_cellchat3_plot_height, dpi = input$s_cellchat3_plot_dpi, units = "in")
  #   }
  # )
  
  
  output$s_cellchat4_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat1_level()[4]
  })
  # observeEvent(input$download_s_cellchat4_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat4_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat4_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Signaling_patterns", input$s_cellchat4_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat1_level()[[4]], width = input$s_cellchat4_plot_width, height = input$s_cellchat4_plot_height, dpi = input$s_cellchat4_plot_dpi, units = "in")
  #   }
  # )
  
  output$s_cellchat12_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat1_level()[5]
  })
  # observeEvent(input$download_s_cellchat12_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat12_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat12_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat12_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat12_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat12_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat12_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("communication_patterns", input$s_cellchat12_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat1_level()[[5]], width = input$s_cellchat12_plot_width, height = input$s_cellchat12_plot_height, dpi = input$s_cellchat12_plot_dpi, units = "in")
  #   }
  # )
  
  
  output$s_cellchat1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_cellchat1_level()[[6]]),
                                                           options = list(
                                                             scrollX = TRUE,
                                                             pageLength = 10,
                                                             dom = "Blfrtip"
                                                             #bFilter=0
                                                           ),rownames= FALSE, selection = "none"))
  
  output$download_s_cellchat1_table <- downloadHandler(
    filename = function() { 
      paste("cellchat_summary_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_cellchat1_level()[[6]], file)
    }
  )
  
  
##############################################submenu#########################################
observe({
      output$s_cellchat_11 <- renderUI ({
        clusters <- req(datainput_single_multiple_sample_cellchat1_level()[[7]])
        shinyWidgets::pickerInput(
          inputId = "s_cellchat11",
          label = "Select one signaligng pathway for vizualization",
          choices = clusters,
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    })
  
  datainput_single_multiple_sample_cellchat2_level <- eventReactive(input$single_multiple_sample_cellchat2,{
    source("scripts/cellchat2.R")
    datainput_single_multiple_sample_cellchat2(index_single_sample_cellchat2_input = datainput_single_multiple_sample_cellchat1_level()[[8]],  index_s_cellchat11 = input$s_cellchat11, index_s_cellchat12 = input$s_cellchat12)
  })  
  
  output$s_cellchat5_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[1]
  })
  # observeEvent(input$download_s_cellchat5_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat5_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat5_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat5_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat5_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat5_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat5_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Number_of_interactions_circle_plot", input$s_cellchat5_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[1]], width = input$s_cellchat5_plot_width, height = input$s_cellchat5_plot_height, dpi = input$s_cellchat5_plot_dpi, units = "in")
  #   }
  # )
  
  output$s_cellchat6_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[2]
  })
  # observeEvent(input$download_s_cellchat6_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat6_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat6_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat6_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat6_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat6_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat6_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Number_of_interactions_chord_plot", input$s_cellchat6_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[2]], width = input$s_cellchat6_plot_width, height = input$s_cellchat6_plot_height, dpi = input$s_cellchat6_plot_dpi, units = "in")
  #   }
  # )
  
  
  output$s_cellchat7_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[3]
  })
  # observeEvent(input$download_s_cellchat7_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat7_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat7_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat7_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat7_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat7_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat7_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("interaction_heatmap", input$s_cellchat7_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[3]], width = input$s_cellchat7_plot_width, height = input$s_cellchat7_plot_height, dpi = input$s_cellchat7_plot_dpi, units = "in")
  #   }
  # )
  
  output$s_cellchat11_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[4]
  })
  # observeEvent(input$download_s_cellchat11_plot, {
  #   showModal(modalDialog(
  #     title = strong("Download plot"),
  #     numericInput("s_cellchat11_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat11_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
  #     numericInput("s_cellchat11_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
  #     selectInput("s_cellchat11_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
  #     downloadBttn("s_cellchat11_plot_downloadoutput", "Download"),
  #     size = "s",
  #     easyClose = TRUE,
  #     #footer = NULL
  #   ))
  # })
  # output$s_cellchat11_plot_downloadoutput<- downloadHandler(
  #   filename = function(){
  #     paste("Hierachy_plot", input$s_cellchat11_plot_type, sep="")
  #   },
  #   content = function(file){
  #     ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[4]], width = input$s_cellchat11_plot_width, height = input$s_cellchat11_plot_height, dpi = input$s_cellchat11_plot_dpi, units = "in")
  #   }
  # )
  
  output$s_cellchat8_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[5]
  })
  observeEvent(input$download_s_cellchat8_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_cellchat8_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cellchat8_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cellchat8_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_cellchat8_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_cellchat8_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_cellchat8_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Bubble_plot", input$s_cellchat8_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[5]], width = input$s_cellchat8_plot_width, height = input$s_cellchat8_plot_height, dpi = input$s_cellchat8_plot_dpi, units = "in")
    }
  )
  
  output$s_cellchat9_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[6]
  })
  observeEvent(input$download_s_cellchat9_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_cellchat9_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cellchat9_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cellchat9_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_cellchat9_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_cellchat9_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_cellchat9_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Network_Analysis_contribution_Bar_plot", input$s_cellchat9_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[6]], width = input$s_cellchat9_plot_width, height = input$s_cellchat9_plot_height, dpi = input$s_cellchat9_plot_dpi, units = "in")
    }
  )
  
  output$s_cellchat10_plot<-renderPlot({
    datainput_single_multiple_sample_cellchat2_level()[7]
  })
  observeEvent(input$download_s_cellchat10_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_cellchat10_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cellchat10_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_cellchat10_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_cellchat10_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_cellchat10_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_cellchat10_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Interaction_table_", input$s_cellchat10, input$s_cellchat10_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_cellchat2_level()[[7]], width = input$s_cellchat10_plot_width, height = input$s_cellchat10_plot_height, dpi = input$s_cellchat10_plot_dpi, units = "in")
    }
  )
  
  
  
  output$s_cellchat2_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_cellchat2_level()[[8]]),
                                                           options = list(
                                                             scrollX = TRUE,
                                                             pageLength = 10,
                                                             dom = "Blfrtip"
                                                             #bFilter=0
                                                           ),rownames= FALSE, selection = "none"))
  
  output$download_s_cellchat7_table <- downloadHandler(
    filename = function() { 
      paste("cellchat_summary_table_", input$s_cellchat10, '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_cellchat2_level()[[8]], file)
    }
  )  
  
  ######################################################Menu8#####################################################################   
  #############################################Trajectory and Pseudotime analysis#################################################
  
  hideTab(inputId = "Coexpression_tabsets", target = "Transcription factor regulatory network analysis")
  observeEvent(input$link_s_tfrn, {
    showTab(inputId = "Coexpression_tabsets", target = "Transcription factor regulatory network analysis")
  })
  
  #####################################link to next tab###########################     
  observeEvent(input$link_s_tfrn, {
    newvalue <- "Transcription factor regulatory network analysis"
    updateTabsetPanel(session, "Coexpression_tabsets", newvalue)
  })  
  
  
  ###################hidebox#################
  shinyjs::hide("s_trajectory_box1")
  shinyjs::hide("s_trajectory_box2")
  shinyjs::hide("s_trajectory_box3")
  shinyjs::hide("s_trajectory_box4")
  shinyjs::hide("s_trajectory_box5")
  shinyjs::hide("s_trajectory_box6")
  shinyjs::hide("s_trajectory_box7")
  shinyjs::hide("s_trajectory_box8")
  shinyjs::hide("s_trajectory18")
  
  
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_trajectory_box0")
      shinyjs::show("s_trajectory_box1")
    }
    else{
      shinyjs::show("s_trajectory_box0")
      shinyjs::hide("s_trajectory_box1") 
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_trajectory_box0")
      shinyjs::show("s_trajectory_box1")
    }
    else{
      shinyjs::show("s_trajectory_box0")
      shinyjs::hide("s_trajectory_box1") 
    }
  })
  
  observeEvent(input$single_multiple_sample_trajectory1,{
    shinyjs::show("s_trajectory_box2")
    shinyjs::show("s_trajectory_box3")
    shinyjs::hide("s_trajectory_box4")
    shinyjs::hide("s_trajectory_box5")
    shinyjs::hide("s_trajectory_box6")
    shinyjs::hide("s_trajectory_box7")
    shinyjs::hide("s_trajectory_box8")
  })
  
  observeEvent(input$single_multiple_sample_trajectory2,{
    shinyjs::show("s_trajectory_box2")
    shinyjs::show("s_trajectory_box3")
    shinyjs::show("s_trajectory_box4")
    shinyjs::show("s_trajectory_box5")
    shinyjs::hide("s_trajectory_box6")
    shinyjs::hide("s_trajectory_box7")
    shinyjs::hide("s_trajectory_box8")
  })
  
  observeEvent(input$single_multiple_sample_trajectory3,{
    shinyjs::show("s_trajectory_box2")
    shinyjs::show("s_trajectory_box3")
    shinyjs::show("s_trajectory_box4")
    shinyjs::show("s_trajectory_box5")
    shinyjs::show("s_trajectory_box6")
    shinyjs::show("s_trajectory_box7")
    shinyjs::hide("s_trajectory_box8")
  })
  
  observeEvent(input$single_multiple_sample_trajectory4,{
    shinyjs::show("s_trajectory_box2")
    shinyjs::show("s_trajectory_box3")
    shinyjs::show("s_trajectory_box4")
    shinyjs::show("s_trajectory_box5")
    shinyjs::show("s_trajectory_box6")
    shinyjs::show("s_trajectory_box7")
    shinyjs::show("s_trajectory_box8")
  })  
  
  
  observe({
    
    if (input$s_trajectory1 == "multiple_sample" & input$s_trajectory2 == "seurat_clusters"){
      output$s_trajectory_10 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_trajectory10",
          label = "Select one cluster as root",
          choices = sort(clusters),
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_trajectory1 == "multiple_sample" & input$s_trajectory2 == "predicted"){
      output$s_trajectory_10 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_trajectory10",
          label = "Select one cluster as root",
          choices = sort(clusters),
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }  
    else if (input$s_trajectory1 == "multiple_sample_subclustering" & input$s_trajectory2 == "seurat_clusters"){
      output$s_trajectory_10 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_trajectory10",
          label = "Select one cluster as root",
          choices = sort(clusters),
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_trajectory1 == "multiple_sample_subclustering" & input$s_trajectory2 == "predicted"){
      output$s_trajectory_10 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_trajectory10",
          label = "Select one cluster as root",
          choices = sort(clusters),
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }
  })
  
  
  
  datainput_single_multiple_sample_trajectory1_level <- eventReactive(input$single_multiple_sample_trajectory1,{
    source("scripts/trajectory1.R")
    datainput_single_multiple_sample_trajectory1(index_multiple_sample_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_s_trajectory1 = input$s_trajectory1, index_s_trajectory2 = input$s_trajectory2, index_s_trajectory3 = input$s_trajectory3, index_s_trajectory4 = input$s_trajectory4, index_s_trajectory5 = input$s_trajectory5, index_s_trajectory6 = input$s_trajectory6, index_s_trajectory7 = input$s_trajectory7, index_s_trajectory8 = input$s_trajectory8)
  })  
  
  
  output$s_trajectory1_plot<-renderPlot({
    datainput_single_multiple_sample_trajectory1_level()[3]
  })
  
  observeEvent(input$download_s_trajectory1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_trajectory1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_trajectory1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_trajectory1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_trajectory1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Trajectory_Plot",  input$s_trajectory1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_trajectory1_level()[[3]], width = input$s_trajectory1_plot_width, height = input$s_trajectory1_plot_height, dpi = input$s_trajectory1_plot_dpi, units = "in")
    }
  )
  
  
  ########submenu2#########
  datainput_single_multiple_sample_trajectory2_level <- eventReactive(input$single_multiple_sample_trajectory2,{
    source("scripts/trajectory2.R")
    datainput_single_multiple_sample_trajectory2(index_trajectory2_input1 = datainput_single_multiple_sample_trajectory1_level()[[1]], index_trajectory2_input2 = datainput_single_multiple_sample_trajectory1_level()[[2]],  index_trajectory2_multiple_sample_input2 = datainput_multiple_celltype_level()[[4]], index_trajectory2_subclustering_multiple_sample_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_s_trajectory1 = input$s_trajectory1, index_s_trajectory2 = input$s_trajectory2, index_s_trajectory10 = input$s_trajectory10, index_s_trajectory11 = input$s_trajectory11, index_s_trajectory12 = input$s_trajectory12, index_s_trajectory13 = input$s_trajectory13, index_s_trajectory14 = input$s_trajectory14)
  })  
  
  output$s_trajectory2_plot<-renderPlot({
    datainput_single_multiple_sample_trajectory2_level()[3]
  })
  
  observeEvent(input$download_s_trajectory2_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_trajectory2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_trajectory2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_trajectory2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_trajectory2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("cells_in_Pseudotime",  input$s_trajectory2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_trajectory2_level()[[3]], width = input$s_trajectory2_plot_width, height = input$s_trajectory2_plot_height, dpi = input$s_trajectory2_plot_dpi, units = "in")
    }
  )
  
  output$s_trajectory3_plot<-renderPlot({
    datainput_single_multiple_sample_trajectory2_level()[4]
  })
  
  observeEvent(input$download_s_trajectory3_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_trajectory3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 16, width = "300px"),
      numericInput("s_trajectory3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_trajectory3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_trajectory3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_trajectory3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Cells_ordered_by_Monocle3_Pseudotime",  input$s_trajectory3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_trajectory2_level()[[4]], width = input$s_trajectory3_plot_width, height = input$s_trajectory3_plot_height, dpi = input$s_trajectory3_plot_dpi, units = "in")
    }
  )
  
  ########submenu3#########
  datainput_single_multiple_sample_trajectory3_level <- eventReactive(input$single_multiple_sample_trajectory3,{
    source("scripts/trajectory3.R")
    datainput_single_multiple_sample_trajectory3(index_trajectory3_input1 = datainput_single_multiple_sample_trajectory2_level()[[1]], index_trajectory3_input2 = datainput_single_multiple_sample_trajectory2_level()[[2]], index_s_trajectory15 = input$s_trajectory15, index_s_trajectory16 = input$s_trajectory16)
  })  
  
  output$s_trajectory4_plot<-renderPlot({
    datainput_single_multiple_sample_trajectory3_level()[3]
  })
  
  observeEvent(input$download_s_trajectory4_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_trajectory4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_trajectory4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_trajectory4_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  
  output$s_trajectory4_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("FeaturePlot_with_Pseudotime",  input$s_trajectory4_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_trajectory3_level()[[3]], width = input$s_trajectory4_plot_width, height = input$s_trajectory4_plot_height, dpi = input$s_trajectory4_plot_dpi, units = "in")
    }
  )
  
  output$s_trajectory1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_trajectory3_level()[[4]]),
                                                             options = list(
                                                               scrollX = TRUE,
                                                               pageLength = 10,
                                                               dom = "Blfrtip"
                                                               #bFilter=0
                                                             ),rownames= FALSE, selection = "none"))
  
  output$download_s_trajectory1_table <- downloadHandler(
    filename = function() { 
      paste("Genes_that_change_as_a_function_of_pseudotime", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_trajectory3_level()[[4]], file)
    }
  )
  
  ########submenu4#########
  observe({
    if(input$s_trajectory17 == "gene_name_list"){
      shinyjs::show("s_trajectory18")
    }
    else {
      shinyjs::hide("s_trajectory18")
    }
  }) 
  
  
  datainput_single_multiple_sample_trajectory4_level <- eventReactive(input$single_multiple_sample_trajectory4,{
    source("scripts/trajectory4.R")
    datainput_single_multiple_sample_trajectory4(index_trajectory4_input1 = datainput_single_multiple_sample_trajectory3_level()[[1]], index_trajectory4_input2 = datainput_single_multiple_sample_trajectory3_level()[[2]], index_trajectory4_input3 = datainput_single_multiple_sample_trajectory3_level()[[4]], index_s_trajectory17 = input$s_trajectory17, index_s_trajectory18 = input$s_trajectory18)
  })  
  
  output$s_trajectory5_plot<-renderPlot({
    datainput_single_multiple_sample_trajectory4_level()[3]
  })
  
  observeEvent(input$download_s_trajectory5_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_trajectory5_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory5_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_trajectory5_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_trajectory5_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_trajectory5_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_trajectory5_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("FeaturePlot_with_Pseudotime_for_selected_genes",  input$s_trajectory5_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_trajectory4_level()[[3]], width = input$s_trajectory5_plot_width, height = input$s_trajectory5_plot_height, dpi = input$s_trajectory5_plot_dpi, units = "in")
    }
  )
  
  
 
  
  ######################################################Menu9#####################################################################   
  ####################################################hdwgcna####################################################################  
  ###################hidebox#################
  shinyjs::hide("s_hdwgcna_box1")
  shinyjs::hide("s_hdwgcna_box2")
  shinyjs::hide("s_hdwgcna_box3")
  shinyjs::hide("s_hdwgcna_box4")
  shinyjs::hide("s_hdwgcna_box5")
  shinyjs::hide("s_hdwgcna_box6")
  shinyjs::hide("s_hdwgcna_box7")
  shinyjs::hide("s_hdwgcna_box8")
  shinyjs::hide("s_hdwgcna_box9")
  shinyjs::hide("s_hdwgcna_box10")
  
  observeEvent(input$multiple_sample_celltype, {
    if(input$m_marker1 == 1){
      shinyjs::hide("s_hdwgcna_box0")
      shinyjs::show("s_hdwgcna_box1")
      shinyjs::show("s_hdwgcna_box2")
      shinyjs::show("s_hdwgcna_box3")
      shinyjs::show("s_hdwgcna_box4")
      shinyjs::show("s_hdwgcna_box5")
    }
    else{
      shinyjs::show("s_hdwgcna_box0")
      shinyjs::hide("s_hdwgcna_box1") 
      shinyjs::hide("s_hdwgcna_box2")
      shinyjs::hide("s_hdwgcna_box3")
      shinyjs::hide("s_hdwgcna_box4")
      shinyjs::hide("s_hdwgcna_box5")
    }
  })
  observeEvent(input$subclustering_multiple_sample_celltype, {
    if(input$m_subclustering_marker1 == 1){
      shinyjs::hide("s_hdwgcna_box0")
      shinyjs::show("s_hdwgcna_box1")
      shinyjs::show("s_hdwgcna_box2")
      shinyjs::show("s_hdwgcna_box3")
      shinyjs::show("s_hdwgcna_box4")
      shinyjs::show("s_hdwgcna_box5")
    }
    else{
      shinyjs::show("s_hdwgcna_box0")
      shinyjs::hide("s_hdwgcna_box1") 
      shinyjs::hide("s_hdwgcna_box2")
      shinyjs::hide("s_hdwgcna_box3")
      shinyjs::hide("s_hdwgcna_box4")
      shinyjs::hide("s_hdwgcna_box5")
    }
  })
  
  observeEvent(input$single_multiple_sample_hdwgcna,{
    shinyjs::show("s_hdwgcna_box6")
    shinyjs::show("s_hdwgcna_box7")
    shinyjs::show("s_hdwgcna_box8")
    shinyjs::show("s_hdwgcna_box9")
    shinyjs::show("s_hdwgcna_box10")
  })
  
  observe({
    if (input$s_hdwgcna1 == "multiple_sample" & input$s_hdwgcna2 == "seurat_clusters"){
      output$s_hdwgcna_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_hdwgcna3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_hdwgcna1 == "multiple_sample" & input$s_hdwgcna2 == "predicted"){
      output$s_hdwgcna_3 <- renderUI ({
        clusters <- req(datainput_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_hdwgcna3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }  
    else if (input$s_hdwgcna1 == "multiple_sample_subclustering" & input$s_hdwgcna2 == "seurat_clusters"){
      output$s_hdwgcna_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[2]])
        shinyWidgets::pickerInput(
          inputId = "s_hdwgcna3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }
    else if (input$s_hdwgcna1 == "multiple_sample_subclustering" & input$s_hdwgcna2 == "predicted"){
      output$s_hdwgcna_3 <- renderUI ({
        clusters <- req(datainput_subclustering_multiple_celltype_level()[[3]])
        shinyWidgets::pickerInput(
          inputId = "s_hdwgcna3",
          label = "Select one or multiple cluster(s) for analsysis",
          choices = sort(clusters),
          selected = sort(clusters)[1],
          multiple = F,
          options = list(`actions-box` = TRUE))
      })
    }
  })
  
  observeEvent(input$single_multiple_sample_hdwgcna, {
    files_to_delete <- c(
      "www/combined_output.pdf",
      "www/PlotDendrogram.pdf",
      "www/ModuleUMAPPlot.pdf",
      "www/PlotModuleCorrelogram.pdf"
    )
    
    for (file in files_to_delete) {
      full_path <- file.path(getwd(), file)
      if (file.exists(full_path)) {
        unlink(full_path)
        cat("Deleted:", full_path, "\n")
      } else {
        cat("File not found, skipping:", full_path, "\n")
      }
    }
  })
  
  datainput_single_multiple_sample_hdwgcna_level <- eventReactive(input$single_multiple_sample_hdwgcna,{
    source("scripts/hdwgcna.R")
    datainput_single_multiple_sample_hdwgcna(index_multiple_sample_hdwgcna_input = datainput_multiple_celltype_level()[[1]], index_subclustering_multiple_sample_hdwgcna_input = datainput_subclustering_multiple_celltype_level()[[1]], index_multiple_sample_hdwgcna_input2 = datainput_multiple_celltype_level()[[4]], index_subclustering_multiple_sample_hdwgcna_input2 = datainput_subclustering_multiple_celltype_level()[[4]], index_multiple_sample_normalization_method_hdwgcna = input$multiple_sample_normalization_method, index_subclustering_multiple_sample_normalization_method_hdwgcna = input$subclustering_multiple_sample_normalization_method, index_s_hdwgcna1 = input$s_hdwgcna1, index_s_hdwgcna2 = input$s_hdwgcna2, index_s_hdwgcna3 = input$s_hdwgcna3, index_s_hdwgcna4 = input$s_hdwgcna4, index_s_hdwgcna5 = input$s_hdwgcna5, index_s_hdwgcna6 = input$s_hdwgcna6, index_s_hdwgcna7 = input$s_hdwgcna7, index_s_hdwgcna8 = input$s_hdwgcna8, index_s_hdwgcna9 = input$s_hdwgcna9, index_s_hdwgcna10 = input$s_hdwgcna10, index_s_hdwgcna11 = input$s_hdwgcna11, index_s_hdwgcna12 = input$s_hdwgcna12, index_s_hdwgcna13 = input$s_hdwgcna13, index_s_hdwgcna14 = input$s_hdwgcna14)
  })  
  
  output$s_hdwgcna1_plot<-renderPlot({
    datainput_single_multiple_sample_hdwgcna_level()[1]
  })
  observeEvent(input$download_s_hdwgcna1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_hdwgcna1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_hdwgcna1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_hdwgcna1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_hdwgcna1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("umap_plot", input$s_hdwgcna1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_hdwgcna_level()[[1]], width = input$s_hdwgcna1_plot_width, height = input$s_hdwgcna1_plot_height, dpi = input$s_hdwgcna1_plot_dpi, units = "in")
    }
  )
  
  output$s_hdwgcna2_plot<-renderPlot({
    datainput_single_multiple_sample_hdwgcna_level()[2]
  })
  observeEvent(input$download_s_hdwgcna2_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_hdwgcna2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_hdwgcna2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_hdwgcna2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_hdwgcna2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Soft_power_threshold_plots", input$s_hdwgcna2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_hdwgcna_level()[[2]], width = input$s_hdwgcna2_plot_width, height = input$s_hdwgcna2_plot_height, dpi = input$s_hdwgcna2_plot_dpi, units = "in")
    }
  )
  
  # output$text_level_test<- renderText({
  #   paste(datainput_single_multiple_sample_hdwgcna_level()[6])
  # })
  
  output$s_hdwgcna3_plot <- shiny::renderUI({
    # Construct the path to the PDF
    pdf_path <- file.path(datainput_single_multiple_sample_hdwgcna_level()[6],"/www/")
    
    # Ensure the file exists before rendering
    if (file.exists(pdf_path)) {
      tags$iframe(
        src = "PlotDendrogram.pdf",  # Use 'file:///' for local paths
        style = "width:100%; height:600px;",
        frameborder = 1
      )
    } else {
      tags$p("PDF file not found. Please ensure the file exists at the specified path.")
    }
  })
  
 
  output$s_hdwgcna4_plot<-renderPlot({
    datainput_single_multiple_sample_hdwgcna_level()[3]
  })
  observeEvent(input$download_s_hdwgcna4_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_hdwgcna4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_hdwgcna4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_hdwgcna4_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_hdwgcna4_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_ranked_by_eigengene_based_connectivity_kME", input$s_hdwgcna4_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_hdwgcna_level()[[3]], width = input$s_hdwgcna4_plot_width, height = input$s_hdwgcna4_plot_height, dpi = input$s_hdwgcna4_plot_dpi, units = "in")
    }
  )
  
  output$s_hdwgcna5_plot<-renderPlot({
    datainput_single_multiple_sample_hdwgcna_level()[4]
  })
  observeEvent(input$download_s_hdwgcna5_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_hdwgcna5_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna5_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna5_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_hdwgcna5_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_hdwgcna5_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_hdwgcna5_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_feature_plots", input$s_hdwgcna5_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_hdwgcna_level()[[4]], width = input$s_hdwgcna5_plot_width, height = input$s_hdwgcna5_plot_height, dpi = input$s_hdwgcna5_plot_dpi, units = "in")
    }
  )
  
  output$s_hdwgcna6_plot <- renderUI({
    # Construct the path to the PDF
    pdf_path <- file.path(datainput_single_multiple_sample_hdwgcna_level()[6],"/www/")
    
    # Ensure the file exists before rendering
    if (file.exists(pdf_path)) {
      tags$iframe(
        src = "PlotModuleCorrelogram.pdf",  # Use 'file:///' for local paths
        style = "width:100%; height:600px;",
        frameborder = 1
      )
    } else {
      tags$p("PDF file not found. Please ensure the file exists at the specified path.")
    }
  })
  
  output$s_hdwgcna7_plot<-renderPlot({
    datainput_single_multiple_sample_hdwgcna_level()[5]
  })
  observeEvent(input$download_s_hdwgcna7_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_hdwgcna7_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna7_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_hdwgcna7_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_hdwgcna7_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_hdwgcna7_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_hdwgcna7_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_with_Seurats_dot_plot", input$s_hdwgcna7_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_hdwgcna_level()[[5]], width = input$s_hdwgcna7_plot_width, height = input$s_hdwgcna7_plot_height, dpi = input$s_hdwgcna7_plot_dpi, units = "in")
    }
  )
  
  
    output$s_hdwgcna7_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_with_Seurats_dotPlot", input$s_hdwgcna7_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_hdwgcna_level()[[5]], width = input$s_hdwgcna7_plot_width, height = input$s_hdwgcna7_plot_height, dpi = input$s_hdwgcna7_plot_dpi, units = "in")
    }
  )
  
  output$s_hdwgcna8_plot <- renderUI({
     # Construct the path to the PDF
    pdf_path <- file.path(datainput_single_multiple_sample_hdwgcna_level()[6],"/www/")
    
    
    # Ensure the file exists before rendering
    if (file.exists(pdf_path)) {
      tags$iframe(
        src = "combined_output.pdf" ,  
        style = "width:100%; height:600px;",
        frameborder = 1
      )
      
    } else {
      tags$p("PDF file not found. Please ensure the file exists at the specified path.")
    }
    
  })
  
  
  output$s_hdwgcna9_plot <- renderUI({
    # Construct the path to the PDF
    pdf_path <- file.path(datainput_single_multiple_sample_hdwgcna_level()[6],"/www/")

    # Ensure the file exists before rendering
    if (file.exists(pdf_path)) {
      tags$iframe(
        src = "ModuleUMAPPlot.pdf" ,
        style = "width:100%; height:600px;",
        frameborder = 1
      )
    } else {
      tags$p("PDF file not found. Please ensure the file exists at the specified path.")
    }
  })
  
  
  output$s_hdwgcna1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_hdwgcna_level()[[7]]),
                                                       options = list(
                                                         scrollX = TRUE,
                                                         pageLength = 10,
                                                         dom = "Blfrtip"
                                                         #bFilter=0
                                                       ),rownames= FALSE, selection = "none"))
  
  output$download_s_hdwgcna1_table <- downloadHandler(
    filename = function() { 
      paste("Soft_power_threshold_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_hdwgcna_level()[[7]], file)
    }
  )
  
  output$s_hdwgcna2_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_hdwgcna_level()[[8]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            dom = "Blfrtip"
                                                            #bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_s_hdwgcna2_table <- downloadHandler(
    filename = function() { 
      paste("Module_assignment_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_hdwgcna_level()[[8]], file)
    }
  )
  
  output$s_hdwgcna3_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_hdwgcna_level()[[9]]),
                                                          options = list(
                                                            scrollX = TRUE,
                                                            pageLength = 10,
                                                            dom = "Blfrtip"
                                                            #bFilter=0
                                                          ),rownames= FALSE, selection = "none"))
  
  output$download_s_hdwgcna3_table <- downloadHandler(
    filename = function() { 
      paste("Top_N_hub_genes", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_hdwgcna_level()[[9]], file)
    }
  )

###################save object file###################
  output$s_hdwgcna <- downloadHandler(
    filename = function(){
      paste("Co_expression_network_analysis.RDS")
    },
    content = function(file){
      saveRDS(datainput_single_multiple_sample_hdwgcna_level()[10], file= file)
    }
  )

  ######################################################Menu9.2#####################################################################   
  ####################################################TFs####################################################################
  shinyjs::hide("s_tfrn_box5")
  shinyjs::hide("s_tfrn_box6")
  shinyjs::hide("s_tfrn_box7")
  
  observeEvent(input$single_multiple_sample_tfrn1, {
    shinyjs::show("s_tfrn_box5")
    shinyjs::show("s_tfrn_box6")
    shinyjs::hide("s_tfrn_box7")
  })
  observeEvent(input$single_multiple_sample_tfrn2, {
    shinyjs::show("s_tfrn_box7")
  })
  
  datainput_single_multiple_sample_tfrn1_level <- eventReactive(input$single_multiple_sample_tfrn1,{
    source("scripts/tfrn1.R")
    datainput_single_multiple_sample_tfrn1(index_multiple_sample_tfrn1_input = datainput_single_multiple_sample_hdwgcna_level()[[10]], index_s_tfrn1 = input$s_tfrn1, index_s_tfrn2 = input$s_tfrn2, index_s_tfrn3 = input$s_tfrn3, index_s_tfrn4 = input$s_tfrn4, index_s_tfrn5 = input$s_tfrn5, index_s_tfrn6 = input$s_tfrn6, index_s_tfrn7 = input$s_tfrn7, index_s_tfrn8 = input$s_tfrn8)
  })  
  
  output$s_tfrn1_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn1_level()[1]
  })
  observeEvent(input$download_s_tfrn1_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn1_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn1_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn1_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn1_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn1_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn1_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_regulatory_network_plot_Positive", input$s_tfrn1_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn1_level()[[1]], width = input$s_tfrn1_plot_width, height = input$s_tfrn1_plot_height, dpi = input$s_tfrn1_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn2_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn1_level()[2]
  })
  observeEvent(input$download_s_tfrn2_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn2_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn2_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn2_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn2_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn2_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn2_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_regulatory_network_plot_Negative", input$s_tfrn2_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn1_level()[[2]], width = input$s_tfrn2_plot_width, height = input$s_tfrn2_plot_height, dpi = input$s_tfrn2_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn3_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn1_level()[3]
  })
  observeEvent(input$download_s_tfrn3_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn3_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn3_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn3_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn3_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn3_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn3_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_regulatory_network_plot_both", input$s_tfrn3_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn1_level()[[3]], width = input$s_tfrn3_plot_width, height = input$s_tfrn3_plot_height, dpi = input$s_tfrn3_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn4_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn1_level()[4]
  })
  observeEvent(input$download_s_tfrn4_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn4_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn4_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn4_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn4_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn4_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn4_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Module_regulatory_network_plot_Module_UMAP", input$s_tfrn4_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn1_level()[[4]], width = input$s_tfrn4_plot_width, height = input$s_tfrn4_plot_height, dpi = input$s_tfrn4_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn1_table<- renderDataTable(DT::datatable((datainput_single_multiple_sample_tfrn1_level()[[5]]),
                                                       options = list(
                                                         scrollX = TRUE,
                                                         pageLength = 10,
                                                         dom = "Blfrtip"
                                                         #bFilter=0
                                                       ),rownames= FALSE, selection = "none"))
  
  output$download_s_tfrn1_table <- downloadHandler(
    filename = function() { 
      paste("TF_network_table", '.csv', sep='') },
    content = function(file){
      write.csv(datainput_single_multiple_sample_tfrn1_level()[[5]], file)
    }
  )   
  
################################################################################################################################  
######################################################submenu################################################################  
  observe({
    output$s_tfrn_11 <- renderUI ({
      clusters <- req(datainput_single_multiple_sample_tfrn1_level()[[6]])
      shinyWidgets::pickerInput(
        inputId = "s_tfrn11",
        label = "Select one TFs",
        choices = clusters,
        multiple = F,
        options = list(`actions-box` = TRUE))
    })
  })
  
  
  datainput_single_multiple_sample_tfrn2_level <- eventReactive(input$single_multiple_sample_tfrn2,{
    source("scripts/tfrn2.R")
    datainput_single_multiple_sample_tfrn2(index_multiple_sample_tfrn2_input = datainput_single_multiple_sample_tfrn1_level()[[7]], index_multiple_sample_tfrn2_input2 = datainput_single_multiple_sample_tfrn1_level()[[8]],index_multiple_sample_tfrn2_input3 = datainput_single_multiple_sample_tfrn1_level()[[9]],index_s_tfrn11 = input$s_tfrn11, index_s_tfrn12 = input$s_tfrn12, index_s_tfrn13 = input$s_tfrn13, index_s_tfrn14 = input$s_tfrn14)
  })  
  
  output$s_tfrn11_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn2_level()[1]
  })
  observeEvent(input$download_s_tfrn11_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn11_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn11_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn11_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn11_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn11_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn11_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Feature_plot_of_selected_TF", input$s_tfrn11_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn2_level()[[1]], width = input$s_tfrn11_plot_width, height = input$s_tfrn11_plot_height, dpi = input$s_tfrn11_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn12_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn2_level()[2]
  })
  observeEvent(input$download_s_tfrn12_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn12_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn12_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn12_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn12_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn12_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn12_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("Top_target_genes_within_TF_regulons", input$s_tfrn12_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn2_level()[[2]], width = input$s_tfrn12_plot_width, height = input$s_tfrn12_plot_height, dpi = input$s_tfrn12_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn13_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn2_level()[3]
  })
  observeEvent(input$download_s_tfrn13_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn13_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn13_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn13_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn13_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn13_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn13_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("TF_network_plot_positive", input$s_tfrn13_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn2_level()[[3]], width = input$s_tfrn13_plot_width, height = input$s_tfrn13_plot_height, dpi = input$s_tfrn13_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn14_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn2_level()[4]
  })
  observeEvent(input$download_s_tfrn14_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn14_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn14_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn14_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn14_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn14_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn14_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("TF_network_plot_negative", input$s_tfrn14_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn2_level()[[4]], width = input$s_tfrn14_plot_width, height = input$s_tfrn14_plot_height, dpi = input$s_tfrn14_plot_dpi, units = "in")
    }
  )
  
  output$s_tfrn15_plot<-renderPlot({
    datainput_single_multiple_sample_tfrn2_level()[5]
  })
  observeEvent(input$download_s_tfrn15_plot, {
    showModal(modalDialog(
      title = strong("Download plot"),
      numericInput("s_tfrn15_plot_height", label = h5("Figure height (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn15_plot_width", label = h5("Figure width (upto 49 inces)"), value = 8, width = "300px"),
      numericInput("s_tfrn15_plot_dpi", label = h5("Figure resolution (dpi:72 to 300)"), value = 300, width = "300px"),
      selectInput("s_tfrn15_plot_type", label = "Image format", choices = list("JPG" = ".jpg", "TIFF" =".tiff", "PDF" = ".pdf",  "SVG" = ".svg", "BMP" = ".bmp", "EPS" = ".eps", "PS" = ".ps"), selected = ".jpg"),
      downloadBttn("s_tfrn15_plot_downloadoutput", "Download"),
      size = "s",
      easyClose = TRUE,
      #footer = NULL
    ))
  })
  output$s_tfrn15_plot_downloadoutput<- downloadHandler(
    filename = function(){
      paste("TF_network_plot_both", input$s_tfrn15_plot_type, sep="")
    },
    content = function(file){
      ggsave(file,plot = datainput_single_multiple_sample_tfrn2_level()[[5]], width = input$s_tfrn15_plot_width, height = input$s_tfrn15_plot_height, dpi = input$s_tfrn15_plot_dpi, units = "in")
    }
  )
  
  ###################save object file###################
  output$s_tfrn <- downloadHandler(
    filename = function(){
      paste("TF_analysis.RDS")
    },
    content = function(file){
      saveRDS(datainput_single_multiple_sample_tfrn2_level()[6], file= file)
    }
  ) 
  
}
