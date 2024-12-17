#doubletfinder
library(DoubletFinder)
multiple_sample_clustering.split <- SplitObject(multiple_sample_clustering, split.by = "orig.ident") 

# loop through samples to find doublets
for (i in 1:length(multiple_sample_clustering.split)) {
  # print the sample we are on
  print(paste0("orig.ident ",i))
  
  # Pre-process seurat object with standard seurat workflow
  multiple_sample_clustering.sample <- NormalizeData(multiple_sample_clustering.split[[i]])
  multiple_sample_clustering.sample <- FindVariableFeatures(multiple_sample_clustering.sample)
  multiple_sample_clustering.sample <- ScaleData(multiple_sample_clustering.sample)
  multiple_sample_clustering.sample <- RunPCA(multiple_sample_clustering.sample, nfeatures.print = 10)
  
  # Find significant PCs
  stdv <- multiple_sample_clustering.sample[["pca"]]@stdev
  sum.stdv <- sum(multiple_sample_clustering.sample[["pca"]]@stdev)
  percent.stdv <- (stdv / sum.stdv) * 100
  cumulative <- cumsum(percent.stdv)
  co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
  co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                       percent.stdv[2:length(percent.stdv)]) > 0.1), 
              decreasing = T)[1] + 1
  min.pc <- min(co1, co2)
  min.pc
  
  # finish pre-processing
  multiple_sample_clustering.sample <- RunUMAP(multiple_sample_clustering.sample, dims = 1:min.pc)
  multiple_sample_clustering.sample <- FindNeighbors(object = multiple_sample_clustering.sample, dims = 1:min.pc)              
  multiple_sample_clustering.sample <- FindClusters(object = multiple_sample_clustering.sample, resolution = 0.1)
  
  # pK identification (no ground-truth)
  sweep.list <- paramSweep(multiple_sample_clustering.sample, PCs = 1:min.pc, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  
  # # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
  bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
  optimal.pk <- bcmvn.max$pK
  optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
  
  ## Homotypic doublet proportion estimate
  annotations <- multiple_sample_clustering.sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp.poi <- round(optimal.pk * nrow(multiple_sample_clustering.sample@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
  
  # run DoubletFinder
  multiple_sample_clustering.sample <- doubletFinder(seu = multiple_sample_clustering.sample, 
                                                     PCs = 1:min.pc, 
                                                     pK = 0.09,
                                                     nExp = nExp.poi, reuse.pANN = FALSE, sct = FALSE)
  metadata <- multiple_sample_clustering.sample@meta.data
  colnames(metadata)[10] <- "doublet_finder"
  multiple_sample_clustering.sample@meta.data <- metadata 
  
  # subset and save
  multiple_sample_clustering.singlets <- subset(multiple_sample_clustering.sample, doublet_finder == "Singlet")
  multiple_sample_clustering.split[[i]] <- multiple_sample_clustering.singlets
  remove(multiple_sample_clustering.singlets)
}

# converge multiple_sample_clustering.split
multiple_sample_clustering.singlets <- merge(x = multiple_sample_clustering.split[[1]],
                                             y = multiple_sample_clustering.split[2:length(multiple_sample_clustering.split)], , add.cell.ids = names(multiple_sample_clustering.split),project="merged")
multiple_sample_clustering.singlets

# compare before and after doublet finder
multiple_sample_clustering
table(multiple_sample_clustering$orig.ident)

multiple_sample_clustering.singlets

table(multiple_sample_clustering.singlets$orig.ident)
