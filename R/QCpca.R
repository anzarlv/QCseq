# Purpose: Script for QCpca function of QCseq package
# Author: Anzar Alvi
# Date: November 27th, 2025
# Version: 0.1.0
# Bugs and Issues: NA

#' PCA plot labeled by QC score
#'
#' The function performs a principal component analysis (PCA) on key scRNAseq
#' QC metrics including n_genes, n_transcripts, percent_mit, and percent_ribo.
#' PCA is then applied to numeric QC features to reduce the data to 2 main
#' dimensions (PCx and PCy), specified by the users. PC1 and PC2 are default,
#' which capture the greatest sources of variation among the cells.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#' @param pc_x Integer specifying which principal component to use on the
#' x-axis (default is 1, i.e., PC1).
#' @param pc_y Integer specifying which principal component to use on the
#' y-axis (default is 2, i.e., PC2).
#'
#' @return Returns ggplot PCA scatter plot labeled by QC rank score.
#'
#' @examples
#' # Using sample_scRNAseq dataset available with the package
#' data("sample_scRNAseq") # Access sample data
#' QCpca(sample_scRNAseq) # Output PCA plot based on QC scores
#' QCpca(sample_scRNAseq, pc_x = 2, pc_y = 3) # Use PC2 vs PC3
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

QCpca <- function(expr_matrix, pc_x = 1, pc_y = 2) {

  df <- QCseq::QCscore(expr_matrix)

  # Select numeric QC metrics (excluding qc_rank_score and non-numeric columns)
  qc_features <- df[, sapply(df, is.numeric)]
  qc_features <- qc_features[, setdiff(colnames(qc_features), "qc_score_rank"), drop = FALSE]

  # Run PCA (scale = TRUE for unit variance normalization)
  pca_res <- prcomp(qc_features, scale. = TRUE)

  # Prepare data frame for plotting
  # We want 2 PC's (default PC1 and PC2, but user can change it)
  pca_df <- as.data.frame(pca_res$x[, c(pc_x, pc_y), drop = FALSE])
  colnames(pca_df) <- c("PCx", "PCy")
  pca_df$qc_score_rank <- df$qc_score_rank

  # Variance explained (for axis labels)
  var_explained <- round(summary(pca_res)$importance[2, c(pc_x, pc_y)] * 100, 1)

  # PCA plot colored by QC rank score
  # We use ggplot2 (Wickham, 2016) package here
  p <- ggplot2::ggplot(pca_df, aes(x = .data$PCx, y = .data$PCy, color = .data$qc_score_rank)) +
    geom_point(alpha = 0.8, size = 1.8) +
    scale_color_viridis_c(option = "plasma") +
    theme_minimal(base_size = 14) +
    labs(
      title = "PCA of QC Metrics (colored by QC Rank Score)",
      x = paste0("PC", pc_x, " (", var_explained[1], "% variance)"),
      y = paste0("PC", pc_y, " (", var_explained[2], "% variance)"),
      color = "QC Rank Score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
  return(p)
}

# [END]
