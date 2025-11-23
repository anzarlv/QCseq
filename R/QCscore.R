#' Computing unified quality control score
#'
#' Computes a unified quality control score from 3 independent QC metrics
#' (n_genes, n_transcripts, and percent_mit) using median absolute deviation (MAD)
#' standardization and rank-based scores.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#'
#' @return Returns a dataframe where rows are Cells, and columns are:
#' number of genes (n_genes), number of transcripts (n_transcripts),
#' percentage of mitochondrial genes (percent_mit), percentage of ribosomal
#' genes (percent_ribo), and unified QC score (qc_score_rank).
#'
#' @examples
#' # Using sample_scRNAseq dataset available with the package
#' data("sample_scRNAseq") # Access sample data
#' qc_results_with_score <- QCscore(sample_scRNAseq) # QC metrics + Score
#' head(qc_results_with_score) # Display results head
#'
#' @references
#' Bacher, R., Chu, L.-F., Argus, C., Bolin, J. M., Knight, P., Thomson, J.,
#' Stewart, R., & Kendziorski, C. (2021). Enhancing biological signals and
#' detection rates in single-cell RNA-seq experiments with cDNA library
#' equalization. Nucleic Acids Research, 50(2), e12–e12.
#' https://doi.org/10.1093/nar/gkab1071.
#'
#' OpenAI. (2025). ChatGPT. ChatGPT 5; OpenAI. https://chatgpt.com/.
#'
#' R Core Team (2025). _R: A Language and Environment for Statistical
#' Computing_. R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/.
#'
#' @importFrom stats IQR median
#'
#' @export

QCscore <- function(expr_matrix) {
  # Call QCcompute to extract QC metrics
  qc_df <- QCseq::QCcompute(expr_matrix)

  # Keep QC metrics of interest
  qc_sub <- qc_df[, c("n_genes", "n_transcripts", "percent_mit")]

  # Standardize QC metrics and svae in df_scaled dataframe
  # Positive for n_genes & n_transcripts
  # Negative for percent_mit
  df_scaled <- data.frame(
    n_genes =  robust_scale(qc_sub$n_genes),
    n_transcripts = robust_scale(qc_sub$n_transcripts),
    percent_mit = -(robust_scale(qc_sub$percent_mit))
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

# Internal helper
# ChatGPT (OpenAI, 2025) helped with setting up the template for this helper
# Bacher et al.'s (2021) research team inspired me to use MAD standardization
#' @noRd
robust_scale <- function(x) {
  CONSTANT = 1.4826
  m   <- median(x, na.rm = TRUE)
  mad <- mad(x, constant = CONSTANT, na.rm = TRUE)
  # 1.4826 is the constant for
  # the mad function
  # (R documentation)
  if (!is.finite(mad) || mad == 0) {
    iq <- IQR(x, na.rm = TRUE)
    if (!is.finite(iq) || iq == 0) return(scale(x))
    return((x - m) / (iq / 1.349)) # used for scaling
  }
  (x - m) / mad
}
