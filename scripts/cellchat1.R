datainput_single_multiple_sample_cellchat1<- function(index_multiple_sample_cellchat1_input, index_subclustering_multiple_sample_cellchat1_input, index_multiple_sample_cellchat1_input2, index_subclustering_multiple_sample_cellchat1_input2, index_s_cellchat1, index_s_cellchat2, index_s_cellchat3, index_s_cellchat4, index_s_cellchat5, index_s_cellchat6, index_s_cellchat7, index_s_cellchat8, index_s_cellchat9, index_s_cellchat10){
  index_s_cellchat10 <- as.logical(index_s_cellchat10)
  
  if (index_s_cellchat1 == "multiple_sample" & index_s_cellchat2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_multiple_sample_cellchat1_input 
    current_clusters <- levels(single_multiple_sample_clustering)
    # Create new names for clusters (e.g., "Cluster0", "Cluster1", ...)
    new_cluster_names <- paste0("Cluster", current_clusters)
    # Rename clusters in the Seurat object
    names(new_cluster_names) <- current_clusters
    single_multiple_sample_clustering <- RenameIdents(single_multiple_sample_clustering, new_cluster_names)
    levels(single_multiple_sample_clustering)
    single_multiple_sample_clustering$new_cluster_names <- Idents(single_multiple_sample_clustering)
    data.input <- GetAssayData(single_multiple_sample_clustering, assay = "RNA", layer = "data") # normalized data matrix
    #Convert Seurat Object to CellChat Object
    cellchat <- createCellChat(object = single_multiple_sample_clustering, meta = single_multiple_sample_clustering@meta.data, group.by = "new_cluster_names") 
  }
  else if (index_s_cellchat1 == "multiple_sample_subclustering" & index_s_cellchat2 == "seurat_clusters"){
    single_multiple_sample_clustering <- index_subclustering_multiple_cellchat1_sample_input  
    current_clusters <- levels(single_multiple_sample_clustering)
    # Create new names for clusters (e.g., "Cluster0", "Cluster1", ...)
    new_cluster_names <- paste0("Cluster", current_clusters)
    # Rename clusters in the Seurat object
    names(new_cluster_names) <- current_clusters
    single_multiple_sample_clustering <- RenameIdents(single_multiple_sample_clustering, new_cluster_names)
    levels(single_multiple_sample_clustering)
    single_multiple_sample_clustering$new_cluster_names <- Idents(single_multiple_sample_clustering)
    data.input <- GetAssayData(single_multiple_sample_clustering, assay = "RNA", layer = "data") # normalized data matrix
    #Convert Seurat Object to CellChat Object
    cellchat <- createCellChat(object = single_multiple_sample_clustering, meta = single_multiple_sample_clustering@meta.data, group.by = "new_cluster_names") 
  }
  
  else if (index_s_cellchat1 == "multiple_sample" & index_s_cellchat2 == "predicted"){
    single_multiple_sample_clustering <- index_multiple_sample_cellchat1_input
    data.input <- GetAssayData(single_multiple_sample_clustering, assay = "RNA", layer = "data") # normalized data matrix
    #Convert Seurat Object to CellChat Object
    cellchat <- createCellChat(object = single_multiple_sample_clustering, meta = single_multiple_sample_clustering@meta.data, group.by = index_multiple_sample_cellchat1_input2) 
   }
  else if (index_s_cellchat1 == "multiple_sample_subclustering" & index_s_cellchat2 == "predicted"){
    single_multiple_sample_clustering <- index_subclustering_multiple_sample_cellchat1_input
    data.input <- GetAssayData(single_multiple_sample_clustering, assay = "RNA", layer = "data") # normalized data matrix
    #Convert Seurat Object to CellChat Object
    cellchat <- createCellChat(object = single_multiple_sample_clustering, meta = single_multiple_sample_clustering@meta.data, group.by = index_subclustering_multiple_sample_cellchat1_input2) 
    }
  
  #Set the Database for Cell-Cell Interaction
  if(index_s_cellchat3 == "PPI.human"){
  CellChatDB <- CellChatDB.human  
  cellchat@DB <- CellChatDB
  }
  else if(index_s_cellchat3 == "PPI.mouse"){
    CellChatDB <- CellChatDB.mouse 
    cellchat@DB <- CellChatDB
  }
  
  # Subset the expression data and signaling genes
  cellchat <- subsetData(cellchat)  # This subsets the expression data to include only the signaling genes
  # Identify overexpressed genes and interactions
  cellchat <- identifyOverExpressedGenes(cellchat, thresh.pc = index_s_cellchat4, thresh.fc = index_s_cellchat5, thresh.p = index_s_cellchat6)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # Project the data onto the PPI network
  if(index_s_cellchat3 == "PPI.human"){
  cellchat <- smoothData(cellchat, adj = PPI.human)  # for human, or PPI.mouse for mouse
  }
  else if(index_s_cellchat3 == "PPI.mouse"){
  cellchat <- smoothData(cellchat, adj = PPI.mouse)  
  }
  
  
  # Compute the communication probability
  cellchat <- computeCommunProb(cellchat, type = index_s_cellchat7)
  # Filter communication based on a probability threshold
  cellchat <- filterCommunication(cellchat, min.cells = index_s_cellchat8)
  
  
  # Infer the communication network
  cellchat <- computeCommunProbPathway(cellchat)
  
  
  # Aggregate the communication network
  cellchat <- aggregateNet(cellchat)
  
  
  #Visualize the Results
  #Plot Interaction Networks
  # Visualize the interaction network of all signaling pathways
  groupSize <- as.numeric(table(cellchat@idents))
  
  
  
  
  
  #Signaling Role Analysis
  # Perform a role analysis (sender, receiver, mediator, influencer)
  cellchat <- netAnalysis_computeCentrality(cellchat)
  # Visualize signaling roles using circle plot
  #par(mfrow = c(1,2), xpd=TRUE)
  #grid.newpage() 
  plots601 <- netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= index_s_cellchat10, title.name = "Number of interactions")
  plots602 <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= index_s_cellchat10, title.name ="Interaction weights/strength")
  #plots601+plots602
  #Heatmap of Interactions
  # Plot a heatmap of communication probabilities
  plots603 <- netVisual_heatmap(cellchat, measure = "weight")
  # Extract all communication interactions
  interaction_table <- subsetCommunication(cellchat)
  
  #Systems analysis of cell-cell communication network
  pathways_show_all <- as.data.frame(cellchat@netP$pathways)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups

  plots604 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 20)
  plots605 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 20)
  #index_s_cellchat9
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = index_s_cellchat9)
  cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = index_s_cellchat9)
  plots606 <-netAnalysis_river(cellchat,  pattern = "incoming")
  plots607 <-netAnalysis_river(cellchat,  pattern = "outgoing")
  
  return(list(plot1 = plots601, plot2 = plots602, plot3 = plots603, plot4 = plots604+plots605, plot5 = plots606+plots607, data1 =interaction_table, data3=pathways_show_all, data2=cellchat)) 
  
}
    