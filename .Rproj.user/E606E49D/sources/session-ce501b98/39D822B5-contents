
QCscore <- function(expr_matrix) {
  qc_df <- QCseq::QCcompute(expr_matrix)

  # Step 1: keep QC metrics of interest
  qc_sub <- qc_df[, c("nGenes", "nTranscipts", "percent_mit")]

  # Step 2: robust standardization (positive for nGenes & nTranscripts, negative for mito%)
  df_scaled <- data.frame(
    nGenes =  QCseq::robust_scale(qc_sub$nGenes),
    nTranscipts = QCseq::robust_scale(qc_sub$nTranscipts),
    percent_mit = -(QCseq::robust_scale(qc_sub$percent_mit))
  )
  rownames(df_scaled) <- qc_df$Cell

  # Step 3: rank-based QC score
  rank_mat <- apply(df_scaled, 2, function(col) {
    rk <- rank(col, na.last = "keep", ties.method = "average")
    rk / max(rk, na.rm = TRUE)
  })
  qc_score_rank <- rowMeans(rank_mat, na.rm = TRUE)

  # Step 4: attach back
  qc_df$qc_score_rank <- qc_score_rank[qc_df$Cell]

  return(qc_df)
}


# Helper Functions
robust_scale <- function(x) {
  m   <- median(x, na.rm = TRUE)
  mad <- mad(x, constant = 1.4826, na.rm = TRUE)
  if (!is.finite(mad) || mad == 0) {
    iq <- IQR(x, na.rm = TRUE)
    if (!is.finite(iq) || iq == 0) return(scale(x))
    return((x - m) / (iq / 1.349))
  }
  (x - m) / mad
}
