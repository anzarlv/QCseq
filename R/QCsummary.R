#' Return a clear summary report for the gene expression data.
#'
#' Returns a clear summary report for the gene expression data using the key
#' quality control metrics (nGenes, nTranscripts, percent_mit, and percent_ribo)
#' and the signature QC score that combines the metrics.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#' @param threshold Tpye numeric, sets the QC score threshold for what is an
#' acceptable score. Deafult = 0.5.
#'
#' @return Returns text (type character). It is a clear quality control summary
#' report specific to the quality control metrics extracted form the scRNAseq
#' gene expression matrix.
#'
#' @examples
#' # TODO
#'
#' @references
#' # TODO
#'
#' @export

QCsummary <- function(expr_matrix, threshold = 0.5) {
  df <- QCseq::QCscore(expr_matrix)

  # version 1... just testing out what looks nice
  # TODO... make this more systematic, find a paper with good thresholds?
  qc_features <- df[, sapply(df, is.numeric), drop = FALSE]

  qc_means <- colMeans(qc_features, na.rm = TRUE)
  qc_sds <- apply(qc_features, 2, sd, na.rm = TRUE)

  prop_low <- mean(df$qc_score_rank < threshold, na.rm = TRUE)
  pct_low <- round(prop_low * 100, 1)

  metric_names <- names(qc_means)
  metric_text <- paste(
    paste0(metric_names, " (mean = ", round(qc_means, 2),
           ", SD = ", round(qc_sds, 2), ")"),
    collapse = "; "
  )

  summary_text <- paste0(
    "Quality control summary: Across ", nrow(df), " cells, the mean QC rank score was ",
    round(mean(df$qc_score_rank, na.rm = TRUE), 3),
    " (SD = ", round(sd(df$qc_score_rank, na.rm = TRUE), 3), "). ",
    "Key QC metrics include ", metric_text, ". ",
    "Approximately ", pct_low, "% of cells had a QC rank score below ",
    threshold,
    ", suggesting they may represent low-quality or damaged cells that should be considered for filtering. ",
    "Overall, these results indicate that the dataset has a ",
    ifelse(mean(df$qc_score_rank, na.rm = TRUE) > 0.6, "generally good",
           ifelse(mean(df$qc_score_rank, na.rm = TRUE) > 0.4, "moderate", "poor")),
    " quality profile."
  )
  print(typeof(summary_text))
  return(summary_text)
}


# here we need justification for thresholds... find a paper?
