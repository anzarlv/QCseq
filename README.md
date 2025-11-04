# QCseq

“Quality-Control Sequencing (QCseq)”: An R package providing
sample-level QC scores to assess data quality.

## Description

`QCseq` is an R package developed to compute and visualize quality
control (QC) metrics for single-cell RNA sequencing (scRNAseq) gene
expression data, providing sample-level QC scores to assess data
quality. As of right now, researchers must manually examine multiple
independent QC metrics such as the number of detected genes, total
transcriptions, and mitochondrial RNA percentage to assess sample
quality, a time-consuming and subjective task. QCseq addresses this by
integrating multiple QC metrics into a single, robust QC score ranging
from 0 to 100. This unified score provides a more objective and
reproducible measure of data quality, allowing researchers to quickly
identify low-quality samples prior to any downstream analyses
(e.g. differential gene expression analysis). Unlike some existing
packages which only output metrics independently \[INSERT PACKAGES HERE
AND CITE\], QCseq systematically provides a unified score to reflect
overall sample integrity. The `QCseq` package was developed using
`R version 4.5.1 (2025-06-13)`,
`Platform: aarch64-apple-darwin20 (64-bit)` and
`Running under: macOS Sequoia 15.6.1`.

## Installation

To install the latest version of the package:

    install.packages("devtools")
    library("devtools")
    devtools::install_github("anzarlv/QCseq", build_vignettes = TRUE)
    library("QCseq")

To run the Shiny app:

    # UNDER CONSTRUCTION

## Overview

    ls("package:QCseq")
    data(package = "QCseq") 
    browseVignettes("QCseq")

`QCseq` contains 5 functions.

1.  ***QCcompute*** for calculating independent quality control metrics
    (number of genes, number of transcripts, percentage of mitochondrial
    genes, and percentage of ribosomal genes) from single-cell RNA
    sequencing (scRNAseq) data.

2.  ***QCscore*** for calculating a unified quality control score from
    independent QC metrics using median absolute deviation (MAD)
    standardization and rank-based scores.

3.  ***QCsummary*** for a clear summary report for the gene expression
    data using the key quality control metrics extracted and the
    signature QC score that combines the metrics.

4.  ***QChistogram*** for a publication-level histogram where x-axis is
    quality control score (from QCscore function) and y-axis is
    frequency.

5.  ***QCpca*** for a plot of principal components analysis (PCA)
    applied to numeric QC features to reduce the data to 2 main
    dimensions (PC1 and PC2), which capture the greatest sources of
    variation among the cells.

The package also contains one scRNA sequencing dataset, called
sample\_scRNAseq. Refer to package vignettes for more details. An
overview of the package is illustrated below.

![](./inst/extdata/QCseq_diagram.png)

## Contributions

The author of the package is Anzar Alvi.

The author wrote the *QCcompute* function, which calculates four
independent QC metrics (number of genes, number of transcripts,
percentage of mitochondrial genes, and percentage of ribosomal genes).
The *QCcompute* function makes use of the GetAssayData function from the
`Seurat` (Source) package to properly extract the expression matrix from
the input Seurat object, and the colSums function from the `Matrix`
(Source) package to calculate, the total number of genes, total number
of transcripts, and gene percentages. ChatGPT (Source) was used to set
up the template for the helper function *get\_gene\_percentage* which is
called in the *QCcompute* function. The function was completed by me
after having the template.

The author wrote the *QCscore* function, which computes a unified QC
score using median absolute deviation (MAD) standardization (Source) and
rank-based scoring from 3 independent metrics: 1) number of genes, 2)
number of transcripts, and 3) percentage of mitochondrial genes. The
*QCscore* function makes use of base R packages including the `stats`
package to compute the median and MAD for different QC metrics (called
in the robust\_scale helper function). ChatGPT (Source) was used to set
up the template for the helper function *robust\_scale* which is called
in the *QCscore* function. The function was completed by me after having
the template.

The author wrote the *QCsummary* function, which returns a clear summary
report for the gene expression data using the key quality control
metrics (nGenes, nTranscripts, percent\_mit, and percent\_ribo) and the
signature QC score that combines the metrics. The *QCsummary* function
makes use of the `stats` package to derive insights (such as sd) from
the metrics. No artificial intelligence tools were used for the
development of this function.

The *QChistogram* function returns a publication-level histogram of QC
signature scores for all inputted samples in the assay matrix. The
*QChistogram* function makes use of the ggplot function from the
`ggplot2` package to plot the histogram. No artificial intelligence
tools were used for the development of this function.

The *QCpca* function performs a principal component analysis (PCA) on
key scRNAseq QC metrics. The PCA reduces the data to 2 main dimensions
(PC1 and PC2), which capture the greatest sources of variation among the
cells. The function returns a `ggplot2` gplot PCA scatter plot labeled
by QC rank score. Furthermore, the *QCpca* function makes use of the
prcomp function from the `stats` package.

## References

\[TODO\]

## Acknowledgements

This package was developed as part of an assessment for 2025 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `QCseq` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/TestingPackage/issues).
