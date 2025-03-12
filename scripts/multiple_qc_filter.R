datainput_multiple_qc_filter <- function(index_multiple_qc_input, index_multiple_qc_input1, index_multiple_group_count, index_group1_name, index_group1_samples, index_group2_name, index_group2_samples, index_group3_name, index_group3_samples, index_group4_name, index_group4_samples, index_group5_name, index_group5_samples, index_group6_name, index_group6_samples, index_multiple_sample_min_count, index_multiple_sample_max_count, index_multiple_sample_max_mito_perc){
   multiple_list <- index_multiple_qc_input1

  if (index_multiple_group_count == 2)
  {
    group1_name <- index_group1_name
    group1 <- index_group1_samples
    multiple_group1 <- merge(multiple_list[[group1[1]]], y = multiple_list[group1[-1]], add.cell.ids = names(multiple_list)[group1],project="group1_name")
    multiple_group1@meta.data$condition <- group1_name
    
    
    group2_name <- index_group2_name
    group2 <- index_group2_samples
    multiple_group2 <- merge(multiple_list[[group2[1]]], y = multiple_list[group2[-1]], add.cell.ids = names(multiple_list)[group2],project="group2_name")
    multiple_group2@meta.data$condition <- group2_name
    
    #mering group1
    group1_vs_group2 <- merge(multiple_group1, y = multiple_group2, add.cell.ids = c(group1_name, group2_name), project = paste(group1_name, "_vs_", group2_name, sep=""))
    #before filtering
    groups_table1 <- table(group1_vs_group2$orig.ident) %>% as.data.frame
    colnames(groups_table1) <- c("Samples", "Cell counts before QC")
    groups_table2 <- table(group1_vs_group2$condition) %>% as.data.frame
    colnames(groups_table2) <- c("Groups", "Cell counts before QC")
    group1_vs_group2[["percent.mt"]] <- PercentageFeatureSet(group1_vs_group2, pattern = "^MT-")
    groups_merged <- subset(group1_vs_group2, subset = nFeature_RNA > index_multiple_sample_min_count & nFeature_RNA < index_multiple_sample_max_count & percent.mt < index_multiple_sample_max_mito_perc)
    groups_table3 <- table(groups_merged$orig.ident) %>% as.data.frame
    colnames(groups_table3) <- c("Samples", "Cell counts after QC")
    groups_table4 <- table(groups_merged$condition) %>% as.data.frame
    colnames(groups_table4) <- c("Groups", "Cell counts after QC")
    
    
    sample_count <- inner_join(groups_table1, groups_table3)
    sample_count_bar <- reshape2::melt(sample_count,id.vars = 1)
    colnames(sample_count_bar) <- c("Samples", "variable", "Cell_counts")
    
    group_count <- inner_join(groups_table2, groups_table4) 
    group_count_bar <- reshape2::melt(group_count,id.vars = 1)
    colnames(group_count_bar) <- c("Groups", "variable", "Cell_counts")  
  }
  
  else if (index_multiple_group_count == 1)
  {
    group1_name <- index_group1_name
    group1 <- index_group1_samples
    group1_vs_group2 <- merge(multiple_list[[group1[1]]], y = multiple_list[group1[-1]], add.cell.ids = names(multiple_list)[group1],project="group1_name")
    group1_vs_group2@meta.data$condition <- group1_name
           
    #before filtering
    groups_table1 <- table(group1_vs_group2$orig.ident) %>% as.data.frame
    colnames(groups_table1) <- c("Samples", "Cell counts before QC")
    groups_table2 <- table(group1_vs_group2$condition) %>% as.data.frame
    colnames(groups_table2) <- c("Groups", "Cell counts before QC")
    group1_vs_group2[["percent.mt"]] <- PercentageFeatureSet(group1_vs_group2, pattern = "^MT-")
    groups_merged <- subset(group1_vs_group2, subset = nFeature_RNA > index_multiple_sample_min_count & nFeature_RNA < index_multiple_sample_max_count & percent.mt < index_multiple_sample_max_mito_perc)
    groups_table3 <- table(groups_merged$orig.ident) %>% as.data.frame
    colnames(groups_table3) <- c("Samples", "Cell counts after QC")
    groups_table4 <- table(groups_merged$condition) %>% as.data.frame
    colnames(groups_table4) <- c("Groups", "Cell counts after QC")
    
    
    sample_count <- inner_join(groups_table1, groups_table3)
    sample_count_bar <- reshape2::melt(sample_count,id.vars = 1)
    colnames(sample_count_bar) <- c("Samples", "variable", "Cell_counts")
    
    group_count <- inner_join(groups_table2, groups_table4) 
    group_count_bar <- reshape2::melt(group_count,id.vars = 1)
    colnames(group_count_bar) <- c("Groups", "variable", "Cell_counts")  
  }
  
 else if (index_multiple_group_count == 3)
  {
    group1_name <- index_group1_name
    group1 <- index_group1_samples
    multiple_group1 <- merge(multiple_list[[group1[1]]], y = multiple_list[group1[-1]], add.cell.ids = names(multiple_list)[group1],project="group1_name")
    multiple_group1@meta.data$condition <- group1_name
    
    
    group2_name <- index_group2_name
    group2 <- index_group2_samples
    multiple_group2 <- merge(multiple_list[[group2[1]]], y = multiple_list[group2[-1]], add.cell.ids = names(multiple_list)[group2],project="group2_name")
    multiple_group2@meta.data$condition <- group2_name
    
    group3_name <- index_group3_name
    group3 <- index_group3_samples
    multiple_group3 <- merge(multiple_list[[group3[1]]], y = multiple_list[group3[-1]], add.cell.ids = names(multiple_list)[group3],project="group3_name")
    multiple_group3@meta.data$condition <- group3_name
    
    #mering group1
    group1_vs_group2_vs_group3 <- merge(multiple_group1, y = c(multiple_group2, multiple_group3), add.cell.ids = c(group1_name, group2_name, group3_name), project = paste(group1_name, "_vs_", group2_name, "_vs_", group3_name, sep=""))
    #before filtering
    groups_table1 <- table(group1_vs_group2_vs_group3$orig.ident) %>% as.data.frame
    colnames(groups_table1) <- c("Samples", "Cell counts before QC")
    groups_table2 <- table(group1_vs_group2_vs_group3$condition) %>% as.data.frame
    colnames(groups_table2) <- c("Groups", "Cell counts before QC")
    group1_vs_group2_vs_group3[["percent.mt"]] <- PercentageFeatureSet(group1_vs_group2_vs_group3, pattern = "^MT-")
    groups_merged <- subset(group1_vs_group2_vs_group3, subset = nFeature_RNA > index_multiple_sample_min_count & nFeature_RNA < index_multiple_sample_max_count & percent.mt < index_multiple_sample_max_mito_perc)
    groups_table3 <- table(groups_merged$orig.ident) %>% as.data.frame
    colnames(groups_table3) <- c("Samples", "Cell counts after QC")
    groups_table4 <- table(groups_merged$condition) %>% as.data.frame
    colnames(groups_table4) <- c("Groups", "Cell counts after QC")
    
    
    sample_count <- inner_join(groups_table1, groups_table3)
    sample_count_bar <- reshape2::melt(sample_count,id.vars = 1)
    colnames(sample_count_bar) <- c("Samples", "variable", "Cell_counts")
    
    group_count <- inner_join(groups_table2, groups_table4) 
    group_count_bar <- reshape2::melt(group_count,id.vars = 1)
    colnames(group_count_bar) <- c("Groups", "variable", "Cell_counts")  
  }
  else if (index_multiple_group_count == 4)
  {
    group1_name <- index_group1_name
    group1 <- index_group1_samples
    multiple_group1 <- merge(multiple_list[[group1[1]]], y = multiple_list[group1[-1]], add.cell.ids = names(multiple_list)[group1],project="group1_name")
    multiple_group1@meta.data$condition <- group1_name
    
    
    group2_name <- index_group2_name
    group2 <- index_group2_samples
    multiple_group2 <- merge(multiple_list[[group2[1]]], y = multiple_list[group2[-1]], add.cell.ids = names(multiple_list)[group2],project="group2_name")
    multiple_group2@meta.data$condition <- group2_name
    
    group3_name <- index_group3_name
    group3 <- index_group3_samples
    multiple_group3 <- merge(multiple_list[[group3[1]]], y = multiple_list[group3[-1]], add.cell.ids = names(multiple_list)[group3],project="group3_name")
    multiple_group3@meta.data$condition <- group3_name
    
    group4_name <- index_group4_name
    group4 <- index_group4_samples
    multiple_group4 <- merge(multiple_list[[group4[1]]], y = multiple_list[group4[-1]], add.cell.ids = names(multiple_list)[group4],project="group4_name")
    multiple_group4@meta.data$condition <- group4_name
    
    #mering group1
    group1_vs_group2_vs_group3_vs_group4 <- merge(multiple_group1, y = c(multiple_group2, multiple_group3, multiple_group4), add.cell.ids = c(group1_name, group2_name, group3_name, group4_name), project = paste(group1_name, "_vs_", group2_name, "_vs_", group3_name, "_vs_", group4_name, sep=""))
    #before filtering
    groups_table1 <- table(group1_vs_group2_vs_group3_vs_group4$orig.ident) %>% as.data.frame
    colnames(groups_table1) <- c("Samples", "Cell counts before QC")
    groups_table2 <- table(group1_vs_group2_vs_group3_vs_group4$condition) %>% as.data.frame
    colnames(groups_table2) <- c("Groups", "Cell counts before QC")
    group1_vs_group2_vs_group3_vs_group4[["percent.mt"]] <- PercentageFeatureSet(group1_vs_group2_vs_group3_vs_group4, pattern = "^MT-")
    groups_merged <- subset(group1_vs_group2_vs_group3_vs_group4, subset = nFeature_RNA > index_multiple_sample_min_count & nFeature_RNA < index_multiple_sample_max_count & percent.mt < index_multiple_sample_max_mito_perc)
    groups_table3 <- table(groups_merged$orig.ident) %>% as.data.frame
    colnames(groups_table3) <- c("Samples", "Cell counts after QC")
    groups_table4 <- table(groups_merged$condition) %>% as.data.frame
    colnames(groups_table4) <- c("Groups", "Cell counts after QC")
    
    
    sample_count <- inner_join(groups_table1, groups_table3)
    sample_count_bar <- reshape2::melt(sample_count,id.vars = 1)
    colnames(sample_count_bar) <- c("Samples", "variable", "Cell_counts")
    
    group_count <- inner_join(groups_table2, groups_table4) 
    group_count_bar <- reshape2::melt(group_count,id.vars = 1)
    colnames(group_count_bar) <- c("Groups", "variable", "Cell_counts")  
  }
   
   
   else if (index_multiple_group_count == 5)
   {
     group1_name <- index_group1_name
     group1 <- index_group1_samples
     multiple_group1 <- merge(multiple_list[[group1[1]]], y = multiple_list[group1[-1]], add.cell.ids = names(multiple_list)[group1],project="group1_name")
     multiple_group1@meta.data$condition <- group1_name
     
     
     group2_name <- index_group2_name
     group2 <- index_group2_samples
     multiple_group2 <- merge(multiple_list[[group2[1]]], y = multiple_list[group2[-1]], add.cell.ids = names(multiple_list)[group2],project="group2_name")
     multiple_group2@meta.data$condition <- group2_name
     
     group3_name <- index_group3_name
     group3 <- index_group3_samples
     multiple_group3 <- merge(multiple_list[[group3[1]]], y = multiple_list[group3[-1]], add.cell.ids = names(multiple_list)[group3],project="group3_name")
     multiple_group3@meta.data$condition <- group3_name
     
     group4_name <- index_group4_name
     group4 <- index_group4_samples
     multiple_group4 <- merge(multiple_list[[group4[1]]], y = multiple_list[group4[-1]], add.cell.ids = names(multiple_list)[group4],project="group4_name")
     multiple_group4@meta.data$condition <- group4_name
     
     group5_name <- index_group5_name
     group5 <- index_group5_samples
     multiple_group5 <- merge(multiple_list[[group5[1]]], y = multiple_list[group5[-1]], add.cell.ids = names(multiple_list)[group5],project="group5_name")
     multiple_group5@meta.data$condition <- group5_name
     
     #mering group1
     group1_vs_group2_vs_group3_vs_group4_vs_group5 <- merge(multiple_group1, y = c(multiple_group2, multiple_group3, multiple_group4, multiple_group5), add.cell.ids = c(group1_name, group2_name, group3_name, group4_name, group5_name), project = paste(group1_name, "_vs_", group2_name, "_vs_", group3_name, "_vs_", group4_name, "_vs_", group5_name, sep=""))
     #before filtering
     groups_table1 <- table(group1_vs_group2_vs_group3_vs_group4_vs_group5$orig.ident) %>% as.data.frame
     colnames(groups_table1) <- c("Samples", "Cell counts before QC")
     groups_table2 <- table(group1_vs_group2_vs_group3_vs_group4_vs_group5$condition) %>% as.data.frame
     colnames(groups_table2) <- c("Groups", "Cell counts before QC")
     group1_vs_group2_vs_group3_vs_group4_vs_group5[["percent.mt"]] <- PercentageFeatureSet(group1_vs_group2_vs_group3_vs_group4_vs_group5, pattern = "^MT-")
     groups_merged <- subset(group1_vs_group2_vs_group3_vs_group4_vs_group5, subset = nFeature_RNA > index_multiple_sample_min_count & nFeature_RNA < index_multiple_sample_max_count & percent.mt < index_multiple_sample_max_mito_perc)
     groups_table3 <- table(groups_merged$orig.ident) %>% as.data.frame
     colnames(groups_table3) <- c("Samples", "Cell counts after QC")
     groups_table4 <- table(groups_merged$condition) %>% as.data.frame
     colnames(groups_table4) <- c("Groups", "Cell counts after QC")
     
     
     sample_count <- inner_join(groups_table1, groups_table3)
     sample_count_bar <- reshape2::melt(sample_count,id.vars = 1)
     colnames(sample_count_bar) <- c("Samples", "variable", "Cell_counts")
     
     group_count <- inner_join(groups_table2, groups_table4) 
     group_count_bar <- reshape2::melt(group_count,id.vars = 1)
     colnames(group_count_bar) <- c("Groups", "variable", "Cell_counts")  
   }
   
   else if (index_multiple_group_count == 6)
   {
     group1_name <- index_group1_name
     group1 <- index_group1_samples
     multiple_group1 <- merge(multiple_list[[group1[1]]], y = multiple_list[group1[-1]], add.cell.ids = names(multiple_list)[group1],project="group1_name")
     multiple_group1@meta.data$condition <- group1_name
     
     
     group2_name <- index_group2_name
     group2 <- index_group2_samples
     multiple_group2 <- merge(multiple_list[[group2[1]]], y = multiple_list[group2[-1]], add.cell.ids = names(multiple_list)[group2],project="group2_name")
     multiple_group2@meta.data$condition <- group2_name
     
     group3_name <- index_group3_name
     group3 <- index_group3_samples
     multiple_group3 <- merge(multiple_list[[group3[1]]], y = multiple_list[group3[-1]], add.cell.ids = names(multiple_list)[group3],project="group3_name")
     multiple_group3@meta.data$condition <- group3_name
     
     group4_name <- index_group4_name
     group4 <- index_group4_samples
     multiple_group4 <- merge(multiple_list[[group4[1]]], y = multiple_list[group4[-1]], add.cell.ids = names(multiple_list)[group4],project="group4_name")
     multiple_group4@meta.data$condition <- group4_name
     
     group5_name <- index_group5_name
     group5 <- index_group5_samples
     multiple_group5 <- merge(multiple_list[[group5[1]]], y = multiple_list[group5[-1]], add.cell.ids = names(multiple_list)[group5],project="group5_name")
     multiple_group5@meta.data$condition <- group5_name
     
     group6_name <- index_group6_name
     group6 <- index_group6_samples
     multiple_group6 <- merge(multiple_list[[group6[1]]], y = multiple_list[group6[-1]], add.cell.ids = names(multiple_list)[group6],project="group6_name")
     multiple_group6@meta.data$condition <- group6_name
     
     #mering group1
     group1_vs_group2_vs_group3_vs_group4_vs_group5_vs_group6 <- merge(multiple_group1, y = c(multiple_group2, multiple_group3, multiple_group4, multiple_group5, multiple_group6), add.cell.ids = c(group1_name, group2_name, group3_name, group4_name, group5_name, group6_name), project = paste(group1_name, "_vs_", group2_name, "_vs_", group3_name, "_vs_", group4_name, "_vs_", group5_name, "_vs_", group6_name, sep=""))
     #before filtering
     groups_table1 <- table(group1_vs_group2_vs_group3_vs_group4_vs_group5_vs_group6$orig.ident) %>% as.data.frame
     colnames(groups_table1) <- c("Samples", "Cell counts before QC")
     groups_table2 <- table(group1_vs_group2_vs_group3_vs_group4_vs_group5_vs_group6$condition) %>% as.data.frame
     colnames(groups_table2) <- c("Groups", "Cell counts before QC")
     group1_vs_group2_vs_group3_vs_group4_vs_group5_vs_group6[["percent.mt"]] <- PercentageFeatureSet(group1_vs_group2_vs_group3_vs_group4_vs_group5_vs_group6, pattern = "^MT-")
     groups_merged <- subset(group1_vs_group2_vs_group3_vs_group4_vs_group5_vs_group6, subset = nFeature_RNA > index_multiple_sample_min_count & nFeature_RNA < index_multiple_sample_max_count & percent.mt < index_multiple_sample_max_mito_perc)
     groups_table3 <- table(groups_merged$orig.ident) %>% as.data.frame
     colnames(groups_table3) <- c("Samples", "Cell counts after QC")
     groups_table4 <- table(groups_merged$condition) %>% as.data.frame
     colnames(groups_table4) <- c("Groups", "Cell counts after QC")
     
     
     sample_count <- inner_join(groups_table1, groups_table3)
     sample_count_bar <- reshape2::melt(sample_count,id.vars = 1)
     colnames(sample_count_bar) <- c("Samples", "variable", "Cell_counts")
     
     group_count <- inner_join(groups_table2, groups_table4) 
     group_count_bar <- reshape2::melt(group_count,id.vars = 1)
     colnames(group_count_bar) <- c("Groups", "variable", "Cell_counts")  
   }
   
  plots6 <- VlnPlot(groups_merged, features = "nFeature_RNA", ncol = 1)
  plots7 <- VlnPlot(groups_merged, features = "nCount_RNA", ncol = 1)
  plots8 <- VlnPlot(groups_merged, features = "percent.mt", ncol = 1)
  plots9 <- VlnPlot(groups_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1, group.by = "condition")
  
  plots10 <- ggplot(sample_count_bar, aes(x=Samples, y=Cell_counts, fill = variable)) +
    geom_bar(stat="identity", position = position_dodge())+
    geom_text(aes(label=Cell_counts), vjust=1.6, position = position_dodge(0.9), color="white", size=3.5)+
    theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
    theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  plots11 <- ggplot(group_count_bar, aes(x=Groups, y=Cell_counts, fill = variable)) +
    geom_bar(stat="identity", position = position_dodge())+
    geom_text(aes(label=Cell_counts), vjust=1.6, position = position_dodge(0.9), color="white", size=3.5)+
    theme(panel.background = element_blank(), panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(), plot.margin=unit(c(1,1,1,1),"line")) +
    theme(axis.text.x=element_blank())+ guides(fill=guide_legend(title="Cell count"))+theme(axis.text.x = element_text(angle = 90, vjust = 1))
  
  return(list(plot1 = plots6+plots7+plots8, plot2 = plots9, plot3 = plots10, plot4 = plots11, data1 = groups_table3, data2 = groups_table4, data3 = groups_merged))
}