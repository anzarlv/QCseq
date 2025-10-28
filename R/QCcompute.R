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
#' number of genes (nGenes), number of transcripts (nTranscripts),
#' percentage of mitochondrial genes (percent_mit), and percentage of ribosomal
#' genes (percent_ribo).
#'
#' @examples
#' # TODO
#'
#' @references
#' # TODO
#' # ChatGPT
#' # https://ietresearch.onlinelibrary.wiley.com/doi/full/10.1049/nbt2.12115
#' # above is for inspiration on what QC metrics to extract
#'
#' @export
#' @import Seurat
#' @import Matrix

QCcompute <- function(expr_matrix) {
  # Properly extract the expression matrix
  main_df <- Seurat::GetAssayData(expr_matrix, assay = "RNA", slot = "counts")

  # Count number of total detected genes per cell (only non-zero counts)
  nGenes <- Matrix::colSums(main_df > 0)

  # Count total number of transcripts per cell (unique molecular identifiers)
  nTranscripts <- Matrix::colSums(main_df)

  # Calculate mitochondrial genes percentage
  percent_mit <- QCseq::get_gene_percentage(main_df, "^MT-", nTranscripts)

  # Calculate ribosomal genes percentage
  # RPL: large subunit of ribosome
  # RPS: small subunit of ribosome
  percent_ribo <- QCseq::get_gene_percentage(main_df, "^RPL|^RPS", nTranscripts)

  # Create final dataframe that will be returned
  final_qc_df <- data.frame(
    Cell = colnames(main_df),
    nGenes = nGenes,
    nTranscipts = nTranscripts,
    percent_mit = percent_mit,
    percent_ribo = percent_ribo,
    stringsAsFactors = FALSE
  )

  # Set rownames to null
  rownames(final_qc_df) <- NULL

  # Return final output
  return(final_qc_df)

}

# Helper Functions, ChatGPT helped with setting up this helper function
get_gene_percentage <- function(expr_matrix, pattern, total_counts) {
  genes <- grep(pattern, rownames(expr_matrix), value = TRUE)
  if (length(genes) > 0) {
    gene_counts <- Matrix::colSums(expr_matrix[genes, , drop = FALSE])
    percent <- gene_counts / total_counts * 100
  } else {
    percent <- rep(NA, ncol(expr_matrix))
  }
  return(percent)
}
