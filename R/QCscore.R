#' Computing unified quality control score
#'
#' Computes a unified quality control score from 3 independent QC metrics
#' (nGenes, nTranscripts, and percent_mit) using median absolute deviation (MAD)
#' standardization and rank-based scores.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#'
#' @return Returns a dataframe where rows are Cells, and columns are:
#' number of genes (nGenes), number of transcripts (nTranscripts),
#' percentage of mitochondrial genes (percent_mit), percentage of ribosomal
#' genes (percent_ribo), and unified QC score (qc_score_rank).
#'
#' @examples
#' # TODO
#'
#' @references
#' # TODO
#' # ChatGPT
#' # https://pmc.ncbi.nlm.nih.gov/articles/PMC8789062/
#' # above is for inspiration on using MAD standardization
#'
#' @export


QCscore <- function(expr_matrix) {
  # Call QCcompute to extract QC metrics
  qc_df <- QCseq::QCcompute(expr_matrix)

  # Keep QC metrics of interest
  qc_sub <- qc_df[, c("nGenes", "nTranscipts", "percent_mit")]

  # Standardize QC metrics and svae in df_scaled dataframe
  # Positive for nGenes & nTranscripts
  # Negative for percent_mit
  df_scaled <- data.frame(
    nGenes =  QCseq::robust_scale(qc_sub$nGenes),
    nTranscipts = QCseq::robust_scale(qc_sub$nTranscipts),
    percent_mit = -(QCseq::robust_scale(qc_sub$percent_mit))
  )
  rownames(df_scaled) <- qc_df$Cell

  # Rank-based QC scores
  rank_mat <- apply(df_scaled, 2, function(col) {
    # Give the average rank if multiple values are tied
    rk <- rank(col, na.last = "keep", ties.method = "average")
    rk / max(rk, na.rm = TRUE)
  })
  qc_score_rank <- rowMeans(rank_mat, na.rm = TRUE)

  # Attach back
  qc_df$qc_score_rank <- qc_score_rank[qc_df$Cell]

  # Return final dataframe that includes both independent QC metrics + QC score
  return(qc_df)
}

# Helper Functions
robust_scale <- function(x) {
  m   <- median(x, na.rm = TRUE)
  mad <- mad(x, constant = 1.4826, na.rm = TRUE) # 1.4826 is the constant for
                                                 # the mad function
                                                 # (R documentation)
  if (!is.finite(mad) || mad == 0) {
    iq <- IQR(x, na.rm = TRUE)
    if (!is.finite(iq) || iq == 0) return(scale(x))
    return((x - m) / (iq / 1.349)) # used for scaling
  }
  (x - m) / mad
}
