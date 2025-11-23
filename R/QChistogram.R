# Purpose: Script for QChistogram function of QCseq package
# Author: Anzar Alvi
# Date: November 27th, 2025
# Version: 0.1.0
# Bugs and Issues: NA

#' Returns a histogram of quality control scores across all samples.
#'
#' Returns a publication-level histogram where x-axis is quality control score
#' (from QCscore function) and y-axis is frequency. Users can customize colour.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#'
#' @param bar_colour Hex code for colour to allow users to customize the colour
#' of the bars in the histogram, default is #4C72B0 (light blue)
#'
#' @return Returns ggplot histogram plot.
#'
#' @examples
#' # Using sample_scRNAseq dataset available with the package
#' data("sample_scRNAseq") # Access sample data
#' QChistogram(sample_scRNAseq) # Output histogram for QC scores
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
#' @export

QChistogram <- function(expr_matrix, bar_colour = "#4C72B0") {
  main_df <- QCseq::QCscore(expr_matrix)

  # Generate histogram of QC scores
  # We use ggplot2 (Wickham, 2016) package here
  p <- ggplot2::ggplot(main_df, aes(x = .data$qc_score_rank)) +
    geom_histogram(
      bins = 50,
      fill = bar_colour,
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

# [END]
