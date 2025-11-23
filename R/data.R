# Purpose: Describe the details of the sample data provided in the package
# (sample_scRNAseq)
# Author: Anzar Alvi
# Date: November 27th, 2025
# Version: 0.1.0
# Bugs and Issues: NA

#' Single-cell RNA sequencing (scRNAseq) gene counts
#'
#' A dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from
#' 10X Genomics. There are 2,700 single cells that were sequenced on the
#' Illumina NextSeq 500.
#'
#' @source 10x Genomics
#'
#' @format
#' A \code{Seurat} object containing 13,714 genes and 2,700 cells with:
#' \itemize{
#'   \item One active assay named \code{"RNA"}.
#'   \item Within the \code{"RNA"} assay, a single \code{"counts"} layer storing a
#'         sparse matrix of raw UMI counts with genes in rows and cells in columns
#'         (13,714 x 2,700).
#' }
#'
#' Each row of the counts matrix corresponds to a gene, and each column
#' corresponds to a single cell.
#'
#' @examples
#' \dontrun{
#' sample_scRNAseq
#' }
"sample_scRNAseq"

# [END]

