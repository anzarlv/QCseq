#' PCA plot labeled by QC score
#'
#' The function performs a principal component analysis (PCA) on key scRNAseq
#' QC metrics including nGenes, nTranscripts, percent_mit, and percent_ribo.
#' PCA is then applied to numeric QC features to reduce the data to 2 main
#' dimensions (PC1 and PC2), which capture the greatest sources of variation
#' among the cells.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#'
#' @return Returns ggplot PCA scatter plot labeled by QC rank score.
#'
#' @examples
#' # Using sample_scRNAseq dataset available with the package
#' data("sample_scRNAseq") # Access sample data
#' QCpca(sample_scRNAseq) # Output PCA plot based on QC scores
#'
#' @references
#' R Core Team (2025). _R: A Language and Environment for Statistical
#' Computing_. R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/.
#'
#' Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4,
#' https://ggplot2.tidyverse.org.
#'
#' @import ggplot2
#' @importFrom stats prcomp
#'
#' @export

QCpca <- function(expr_matrix) {

  df <- QCseq::QCscore(expr_matrix)

  # Select numeric QC metrics (excluding qc_rank_score and non-numeric columns)
  qc_features <- df[, sapply(df, is.numeric)]
  qc_features <- qc_features[, setdiff(colnames(qc_features), "qc_score_rank"), drop = FALSE]

  # Run PCA (scale = TRUE for unit variance normalization)
  pca_res <- prcomp(qc_features, scale. = TRUE)

  # Prepare data frame for plotting
  pca_df <- as.data.frame(pca_res$x[, 1:2]) # We want PCA1 and PCA2
  pca_df$qc_score_rank <- df$qc_score_rank

  # Variance explained (for axis labels)
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

  # PCA plot colored by QC rank score
  # We use ggplot2 (Wickham, 2016) package here
  p <- ggplot2::ggplot(pca_df, aes(x = .data$PC1, y = .data$PC2, color = .data$qc_score_rank)) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal(base_size = 14) +
    labs(
      title = "PCA of QC Metrics (colored by QC Rank Score)",
      x = paste0("PC1 (", var_explained[1], "% variance)"),
      y = paste0("PC2 (", var_explained[2], "% variance)"),
      color = "QC Rank Score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  return(p)
}
