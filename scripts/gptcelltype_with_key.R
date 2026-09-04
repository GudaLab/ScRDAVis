gptcelltype_with_key <- function(input, api_key, tissuename = NULL, model = "gpt-4", topgenenumber = 10) {
  api_key <- trimws(as.character(api_key)[1])
  on.exit({
    api_key <- NULL
  }, add = TRUE)

  if (is.na(api_key) || !nzchar(api_key)) {
    stop("An OpenAI API key is required for GPTCelltype.", call. = FALSE)
  }

  if (is.list(input) && !is.data.frame(input)) {
    markers_by_cluster <- vapply(input, paste, collapse = ",", FUN.VALUE = character(1))
  } else {
    required_columns <- c("avg_log2FC", "gene", "cluster")
    missing_columns <- setdiff(required_columns, colnames(input))
    if (length(missing_columns) > 0) {
      stop(
        "GPTCelltype marker input is missing required column(s): ",
        paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }

    positive_markers <- input[input$avg_log2FC > 0, , drop = FALSE]
    markers_by_cluster <- tapply(
      positive_markers$gene,
      positive_markers$cluster,
      function(genes) paste0(genes[seq_len(min(length(genes), topgenenumber))], collapse = ",")
    )
  }

  if (length(markers_by_cluster) == 0) {
    stop("No positive marker genes were available for GPTCelltype.", call. = FALSE)
  }

  batch_count <- ceiling(length(markers_by_cluster) / 30)
  batch_ids <- if (batch_count > 1) {
    as.numeric(cut(seq_along(markers_by_cluster), batch_count))
  } else {
    rep(1, length(markers_by_cluster))
  }

  all_results <- lapply(seq_len(batch_count), function(batch_index) {
    marker_indices <- which(batch_ids == batch_index)
    result <- character()
    max_attempts <- 3L

    for (attempt in seq_len(max_attempts)) {
      response <- openai::create_chat_completion(
        model = model,
        messages = list(list(
          role = "user",
          content = paste0(
            "Identify cell types of ", tissuename,
            " cells using the following markers separately for each\n row. ",
            "Only provide the cell type name. Do not show numbers before the name.\n ",
            "Some can be a mixture of multiple cell types.\n",
            paste(markers_by_cluster[marker_indices], collapse = "\n")
          )
        )),
        openai_api_key = api_key
      )

      result <- trimws(strsplit(response$choices[, "message.content"], "\n")[[1]])
      result <- result[nzchar(result)]
      result <- sub("^(?:[-*]\\s*|[0-9]+[.)]\\s*)", "", result, perl = TRUE)

      if (length(result) == length(marker_indices)) {
        break
      }
    }

    if (length(result) != length(marker_indices)) {
      stop(
        "GPTCelltype model '", model, "' returned ", length(result),
        " labels for ", length(marker_indices),
        " clusters after ", max_attempts, " attempts. Please try the run again.",
        call. = FALSE
      )
    }

    names(result) <- names(markers_by_cluster)[marker_indices]
    result
  })

  gsub(",$", "", unlist(all_results))
}
