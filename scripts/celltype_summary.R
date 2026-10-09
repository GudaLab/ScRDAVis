# Shared helpers for the "Cell Type Prediction" tabs (main analysis and sub-clustering).
#
# For any annotation column written by the prediction step (own labels = cell_type,
# ScType = sctype_classification, GPTCelltype = GPTCelltype, SingleR = singleR_labels)
# these helpers build three extra outputs that use the assigned cell type names:
#   celltype_counts  - number of cells per cell type name
#   celltype_markers - marker gene table with the cell type name added to every cluster
#   celltype_heatmap - heatmap of the top marker genes with cells grouped by cell type name
#
# Own labels, ScType and GPTCelltype name whole clusters, so the cluster -> name link is exact.
# SingleR names individual cells, so each cluster is shown with its most frequent SingleR label.

celltype_summary_clean_labels <- function(values) {
  values <- as.character(values)
  values[is.na(values) | !nzchar(trimws(values))] <- "Unassigned"
  values
}

# One row per cluster: the cell type name given to it (or the most frequent one).
celltype_cluster_label_map <- function(object, annotation_column) {
  cluster_ids <- Idents(object)
  cluster_levels <- levels(cluster_ids)
  cluster_ids <- as.character(cluster_ids)
  labels <- celltype_summary_clean_labels(object@meta.data[[annotation_column]])

  rows <- lapply(cluster_levels[cluster_levels %in% cluster_ids], function(cluster_id) {
    cluster_labels <- labels[cluster_ids == cluster_id]
    label_counts <- sort(table(cluster_labels), decreasing = TRUE)
    data.frame(
      cluster = cluster_id,
      cell_type = names(label_counts)[1],
      cells_in_cluster = length(cluster_labels),
      pct_of_cluster = round(100 * as.numeric(label_counts[1]) / length(cluster_labels), 2),
      n_labels = length(label_counts),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Number of cells per cell type name.
build_celltype_count_table <- function(object, annotation_column, label_map = NULL) {
  if (is.null(label_map)) {
    label_map <- celltype_cluster_label_map(object, annotation_column)
  }
  labels <- celltype_summary_clean_labels(object@meta.data[[annotation_column]])
  label_counts <- table(labels)

  # Cell types follow the cluster order; names that never lead a cluster come last.
  remaining <- setdiff(names(sort(label_counts, decreasing = TRUE)), label_map$cell_type)
  label_order <- c(unique(label_map$cell_type), remaining)
  label_order <- label_order[label_order %in% names(label_counts)]

  counts <- as.integer(label_counts[label_order])
  count_table <- data.frame(
    label_order,
    counts,
    round(100 * counts / sum(counts), 2),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  colnames(count_table) <- c("Cell type", "Cell counts", "Percentage of cells")

  # Cluster-level methods: list the Seurat cluster(s) that carry each name.
  if (all(label_map$n_labels == 1)) {
    clusters_per_label <- vapply(label_order, function(label) {
      paste(label_map$cluster[label_map$cell_type == label], collapse = ", ")
    }, character(1))
    count_table <- data.frame(
      count_table[, 1, drop = FALSE],
      `Seurat clusters` = unname(clusters_per_label),
      count_table[, -1, drop = FALSE],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  rownames(count_table) <- NULL
  count_table
}

celltype_summary_valid_markers <- function(markers) {
  is.data.frame(markers) && nrow(markers) > 0 && all(c("cluster", "gene") %in% colnames(markers))
}

# Marker table from "Identify markers in all clusters" with the cell type name of each cluster.
build_celltype_marker_table <- function(markers, label_map) {
  if (!celltype_summary_valid_markers(markers)) {
    return(NULL)
  }
  marker_table <- as.data.frame(markers, stringsAsFactors = FALSE)
  map_index <- match(as.character(marker_table$cluster), label_map$cluster)
  cluster_position <- match("cluster", colnames(marker_table))

  new_columns <- data.frame(cell_type = label_map$cell_type[map_index], stringsAsFactors = FALSE)
  if (any(label_map$n_labels > 1)) {
    # Per-cell methods (SingleR): report how much of the cluster carries the name shown.
    new_columns$cell_type_pct_of_cluster <- label_map$pct_of_cluster[map_index]
  }

  marker_table <- data.frame(
    marker_table[, seq_len(cluster_position), drop = FALSE],
    new_columns,
    marker_table[, -seq_len(cluster_position), drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  rownames(marker_table) <- NULL
  marker_table
}

# Heatmap of the top marker genes per cluster with the cells grouped by cell type name.
build_celltype_marker_heatmap <- function(object, markers, annotation_column, label_map, top_n = 5) {
  if (!celltype_summary_valid_markers(markers) || !"avg_log2FC" %in% colnames(markers)) {
    return(NULL)
  }
  group_values <- celltype_summary_clean_labels(object@meta.data[[annotation_column]])
  group_levels <- c(unique(label_map$cell_type), setdiff(unique(group_values), label_map$cell_type))
  group_levels <- group_levels[group_levels %in% group_values]

  top_markers <- as.data.frame(markers) %>% group_by(cluster) %>% slice_max(n = top_n, order_by = avg_log2FC)
  top_markers <- as.data.frame(top_markers)
  top_markers$cell_type <- label_map$cell_type[match(as.character(top_markers$cluster), label_map$cluster)]
  # Keep the genes of clusters that share a name next to each other.
  top_markers <- top_markers[order(match(top_markers$cell_type, group_levels), match(as.character(top_markers$cluster), label_map$cluster)), , drop = FALSE]
  heatmap_features <- unique(as.character(top_markers$gene))

  heatmap_object <- object
  heatmap_object@meta.data$celltype_heatmap_group <- factor(group_values, levels = group_levels)

  # Leave room above the colour bar for long cell type names (drawn vertically).
  top_margin <- max(30, 6 * max(nchar(group_levels)) + 10)
  DoHeatmap(heatmap_object, features = heatmap_features, angle = 90, size = 4, group.by = "celltype_heatmap_group") +
    theme(plot.margin = margin(t = top_margin, r = 20, b = 10, l = 10, unit = "pt"))
}

build_celltype_summary_outputs <- function(object, annotation_column, markers = NULL) {
  outputs <- list(celltype_counts = NULL, celltype_markers = NULL, celltype_heatmap = NULL)
  safe_build <- function(expr, what) {
    tryCatch(expr, error = function(e) {
      message("Cell type summary: ", what, " could not be created: ", conditionMessage(e))
      NULL
    })
  }

  label_map <- safe_build(celltype_cluster_label_map(object, annotation_column), "cluster to cell type map")
  if (is.null(label_map)) {
    return(outputs)
  }
  outputs["celltype_counts"] <- list(safe_build(build_celltype_count_table(object, annotation_column, label_map), "cell count table"))
  outputs["celltype_markers"] <- list(safe_build(build_celltype_marker_table(markers, label_map), "marker table"))
  outputs["celltype_heatmap"] <- list(safe_build(build_celltype_marker_heatmap(object, markers, annotation_column, label_map), "marker heatmap"))
  outputs
}

# Pads the result of a prediction method to the eight positions the app reads by index
# (object, clusters, labels, column name, dimplot, score table, score plot 1, score plot 2)
# and appends the three named cell type summary outputs.
finalize_celltype_result <- function(result, object, annotation_column, markers = NULL) {
  base_names <- c("data1", "data2", "data3", "text_summary", "plot1", "data2", "plot2", "plot3")
  while (length(result) < length(base_names)) {
    position <- length(result) + 1
    result[position] <- list(NULL)
    names(result)[position] <- base_names[position]
  }
  c(result, build_celltype_summary_outputs(object, annotation_column, markers))
}
