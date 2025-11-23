# Purpose: Script for QCcompute function of QCseq package
# Author: Anzar Alvi
# Date: November 27th, 2025
# Version: 0.1.0
# Bugs and Issues: NA

#' Computing independent quality control metrics
#'
#' Computes the following quality control metrics from a scRNA-seq gene
#' expression matrix: 1) number of detected genes per cell, 2) total transcripts
#' per cell, 3) percentage of mitochondrial genes, and 4) percentage of
#' ribosomal genes.
#'
#' @param expr_matrix A Seurat object, RNA assay, gene expression matrix
#' where columns are cells and rows are genes.
#'
#' @return Returns a dataframe where rows are Cells, and columns are:
#' number of genes (n_genes), number of transcripts (n_transcripts),
#' percentage of mitochondrial genes (percent_mit), and percentage of ribosomal
#' genes (percent_ribo).
#'
#' @examples
#' # Using sample_scRNAseq dataset available with the package
#' data("sample_scRNAseq") # Access sample data
#' qc_results <- QCcompute(sample_scRNAseq) # Compute QC metrics for each sample
#' head(qc_results) # Display results head
#'
#' @references
#' Hao, Y., Stuart, T. A., Kowalski, M. H., Choudhary, S., Hoffman, P., Hartman,
#' A., Srivastava, A., Molla, G., Shaista Madad, Fernandez-Granda, C., & Rahul
#' atija. (2023). Dictionary learning for integrative, multimodal and scalable
#' single-cell analysis. Nature Biotechnology, 42, 293–304.
#' https://doi.org/10.1038/s41587-023-01767-y
#'
#' Lu, J., Sheng, Y., Qian, W., Pan, M., Zhao, X., & Ge, Q. (2023).
#' scRNA‐seq data analysis method to improve analysis performance.
#' Iet Nanobiotechnology, 17(3), 246–256. https://doi.org/10.1049/nbt2.12115.
#'
#' OpenAI. (2025). ChatGPT. ChatGPT 5; OpenAI. https://chatgpt.com/.
#'
#' R Core Team (2025). _R: A Language and Environment for Statistical
#' Computing_. R Foundation for Statistical Computing, Vienna, Austria.
#' https://www.R-project.org/.
#'
#' @export
#' @import Seurat

QCcompute <- function(expr_matrix) {
  # Defensive coding below
  # Input must be a Seurat object with an RNA assay and a counts layer
  if (!inherits(expr_matrix, "Seurat")) {
    stop("`expr_matrix` must be a Seurat object (e.g., like `sample_scRNAseq`).")
  }

  if (!"RNA" %in% Seurat::Assays(expr_matrix)) {
    stop("Seurat object must contain an 'RNA' assay.")
  }

  # This assumes the RNA assay has a 'counts' layer; if not, this will error.
  main_df <- Seurat::GetAssayData(expr_matrix, assay = "RNA", layer = "counts")

  if (!is.matrix(main_df) && !"dgCMatrix" %in% class(main_df)) {
    stop("Extracted RNA counts must be a matrix-like object ",
         "with genes in rows and cells in columns.")
  }

  if (anyNA(main_df)) {
    stop("Expression matrix must not contain NA values; use 0 for undetected genes.")
  }

  # Switch to data frame format
  main_df <- data.frame(main_df)

  # Lu et al.'s (2023) reserach teams helped decide what QC metrics to
  # extract
  # Count number of total detected genes per cell (only non-zero counts)
  n_genes <- colSums(main_df > 0)

  # Count total number of transcripts per cell (unique molecular identifiers)
  n_transcripts <- colSums(main_df)

  # Calculate mitochondrial genes percentage
  percent_mit <- get_gene_percentage(main_df, "^MT-", n_transcripts)

  # Calculate ribosomal genes percentage
  # RPL: large subunit of ribosome
  # RPS: small subunit of ribosome
  percent_ribo <- get_gene_percentage(main_df, "^RPL|^RPS", n_transcripts)

  # Create final dataframe that will be returned
  final_qc_df <- data.frame(
    Cell = colnames(main_df),
    n_genes = n_genes,
    n_transcripts = n_transcripts,
    percent_mit = percent_mit,
    percent_ribo = percent_ribo,
    stringsAsFactors = FALSE
  )

  # Set rownames to null
  rownames(final_qc_df) <- NULL

  # Return final output
  return(final_qc_df)
}

# Internal helper function
# ChatGPT (OpenAI, 2025) helped with setting up this helper
#' @noRd
get_gene_percentage <- function(expr_matrix, pattern, total_counts) {
  genes <- grep(pattern, rownames(expr_matrix), value = TRUE)
  if (length(genes) > 0) {
    gene_counts <- colSums(expr_matrix[genes, , drop = FALSE])
    percent <- gene_counts / total_counts * 100
  } else {
    percent <- rep(NA, ncol(expr_matrix))
  }
  return(percent)
}

# [END]

