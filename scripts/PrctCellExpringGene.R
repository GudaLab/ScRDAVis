# updated 1/31/2020 to accommodate V3.1
# updated 2/4/2020 to output "NA" for genes not detected in certain subgroups

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper,object=object))
	prct1 = unlist(lapply(genes,calc_helper1,object=object))
    result = data.frame(Markers = genes, Cell_count = prct, Cell_proportion = prct1)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  #counts = object[['RNA']]@counts
  counts = object[['RNA']]$counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)
    #sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

calc_helper1 <- function(object,genes){
  #counts = object[['RNA']]@counts
  counts = object[['RNA']]$counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    #sum(counts[genes,]>0)
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

