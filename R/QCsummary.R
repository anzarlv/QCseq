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
#' # Using sample_scRNAseq dataset available with the package
#' data("sample_scRNAseq") # Access sample data
#' qc_summary <- QCsummary(sample_scRNAseq) # Output QC summary paragraph
#' qc_summary # Display the output paragraph
#'
#' @references
#' R Core Team (2025). _R: A Language and Environment for Statistical
#' Computing_. R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/.
#'
#' @importFrom stats sd
#'
#' @export

QCsummary <- function(expr_matrix, threshold = 0.5) {
  # Compute numeric QC metrics using QCscore function
  df <- QCseq::QCscore(expr_matrix)
  qc_features <- df[, sapply(df, is.numeric), drop = FALSE]

  # Calculate mean and standard deviation for each numeric QC metric
  qc_means <- colMeans(qc_features, na.rm = TRUE)
  qc_sds <- apply(qc_features, 2, sd, na.rm = TRUE)

  # Determine proportion (%) of cells below the quality threshold
  prop_low <- mean(df$qc_score_rank < threshold, na.rm = TRUE)
  pct_low <- round(prop_low * 100, 1)

  # Get names for all QC metrics for readable output
  metric_names <- names(qc_means)

  # Summary text: each QC metric with mean and SD
  metric_text <- paste(
    paste0(metric_names, " (mean = ", round(qc_means, 2),
           ", SD = ", round(qc_sds, 2), ")"),
    collapse = "; "
  )

  # Final summary paragraph
  # I used my own judgement to determine that a QC score > 0.6 is "generally
  # good", a QC score between 0.4 and 0.5 is "moderate", and a QC score below
  # 0.4 is "poor".
  summary_text <- paste0(
    "Quality control summary: Across ", nrow(df), " cells, the mean QC rank
    score was ", round(mean(df$qc_score_rank, na.rm = TRUE), 3),
    " (SD = ", round(sd(df$qc_score_rank, na.rm = TRUE), 3), "). ",
    "Key QC metrics include ", metric_text, ". ",
    "Approximately ", pct_low, "% of cells had a QC rank score below ",
    threshold, ", suggesting they may represent low-quality or damaged cells
    that should be considered for filtering. ", "Overall, these results indicate
    that the dataset has a ", ifelse(mean(df$qc_score_rank, na.rm = TRUE) > 0.6,
                                     "generally good",
                                     ifelse(mean(df$qc_score_rank, na.rm = TRUE)
                                            > 0.4, "moderate", "poor")),
    " quality profile."
  )

  return(summary_text)
}
