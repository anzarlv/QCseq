#' Returns a histogram of quality control scores across all samples.
#'
#' Returns a publication-level histogram where x-axis is quality control score
#' (from QCscore function) and y-axis is frequency.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#'
#' @return Returns ggplot histogram plot.
#'
#' @examples
#' # TODO
#'
#' @references
#' # TODO
#'
#' @import ggplot2
#' @export

QChistogram <- function(expr_matrix) {
  main_df <- QCseq::QCscore(expr_matrix)
  print(main_df)
  # generate histogram of QC scores
  p <- ggplot2::ggplot(main_df, aes(x = qc_score_rank)) +
    geom_histogram(
      bins = 50,
      fill = "#4C72B0",
      color = "black",
      alpha = 0.9
    ) +
    labs(
      title = "Distribution of QC Rank Scores Across Cells",
      x = "QC Rank Score",
      y = "Number of Cells"
    ) +
    theme_classic(base_size = 12) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text = element_text(size = 11, face = "plain"),
      axis.line = element_line(colour = "black"),
      panel.grid = element_blank(), # remove gridlines
      axis.ticks = element_line(colour = "black")
    )
  return(p)
}
