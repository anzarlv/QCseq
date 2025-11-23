# QCseq

“Quality-Control Sequencing (QCseq)”: An R package providing
sample-level QC scores to assess data quality.

## Description

`QCseq` is an R package developed to compute and visualize quality
control (QC) metrics for single-cell RNA sequencing (scRNAseq) gene
expression data, providing sample-level QC scores to assess data
quality.

For background, scRNAseq gene expression data lets you quantify the
level of expression of different genes for individual cells. This is
very useful for researchers as the pattern of gene expression can be
used as a “fingerprint” of the state of a sample during the time of
measurement. Details about which genes are relevant to the fingerprint
(this is referred to as genes that are “differentially expressed”) can
give researchers insight into a particular condition of interest or even
biomarkers (Andrews & Hemberg, 2018).

As of right now, researchers must manually examine multiple independent
QC metrics such as the number of detected genes, total transcriptions,
and mitochondrial RNA percentage to assess sample quality, a
time-consuming and subjective task. QCseq addresses this by integrating
multiple QC metrics into a single, robust QC score ranging from 0 to
100. This unified score provides a more objective and reproducible
measure of data quality, allowing researchers to quickly identify
low-quality samples prior to any downstream analyses (e.g. differential
gene expression analysis). Unlike some existing packages which only
output metrics independently such as RNAseqQC (Ziebell, 2024), Seurat
(Hao et al., 2023), or scater (McCarthy et al., 2017), QCseq
systematically provides a unified score to reflect overall sample
integrity. The `QCseq` package was developed using
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
Lu et al.’s (2023) research team helped the author determine the best QC
metrics to extract. The *QCcompute* function makes use of the
GetAssayData function from the `Seurat` (Hao et al., 2023) package to
properly extract the expression matrix from the input Seurat object, and
the colSums function from base R to calculate, the total number of
genes, total number of transcripts, and gene percentages. ChatGPT
(OpenAI, 2025) was used to set up the template for the helper function
*get\_gene\_percentage* which is called in the *QCcompute* function. The
function was completed by me after having the template.

The author wrote the *QCscore* function, which computes a unified QC
score using median absolute deviation (MAD) standardization (Source) and
rank-based scoring from 3 independent metrics: 1) number of genes, 2)
number of transcripts, and 3) percentage of mitochondrial genes. The
idea to use this standardization and rank-based scoring method came from
Bacher et al.’s (2021) research team. The *QCscore* function makes use
of base R packages including the `stats` (R Core Team, 2025) package to
compute the median and MAD for different QC metrics (called in the
robust\_scale helper function). ChatGPT (OpenAI, 2025) was used to set
up the template for the helper function *robust\_scale* which is called
in the *QCscore* function. The function was completed by me after having
the template.

The author wrote the *QCsummary* function, which returns a clear summary
report for the gene expression data using the key quality control
metrics (n\_genes, n\_transcripts, percent\_mit, and percent\_ribo) and
the signature QC score that combines the metrics. The *QCsummary*
function makes use of the `stats` (R Core Team, 2025) package to derive
insights (such as sd) from the metrics. No artificial intelligence tools
were used for the development of this function.

The *QChistogram* function returns a publication-level histogram of QC
signature scores for all inputted samples in the assay matrix. The
*QChistogram* function makes use of the ggplot function from the
`ggplot2` (Wickham, 2016) package to plot the histogram. No artificial
intelligence tools were used for the development of this function.

The *QCpca* function performs a principal component analysis (PCA) on
key scRNAseq QC metrics. The PCA reduces the data to 2 main dimensions
(PC1 and PC2), which capture the greatest sources of variation among the
cells. The function returns a `ggplot2` (Wickham, 2016) gplot PCA
scatter plot labeled by QC rank score. Furthermore, the *QCpca* function
makes use of the prcomp function from the `stats` package.

## References

Andrews, T. S., & Hemberg, M. (2018). Identifying cell populations with
scRNASeq. Molecular Aspects of Medicine, 59, 114–122.
<https://doi.org/10.1016/j.mam.2017.07.002>.

Bacher, R., Chu, L.-F., Argus, C., Bolin, J. M., Knight, P., Thomson,
J., Stewart, R., & Kendziorski, C. (2021). Enhancing biological signals
and detection rates in single-cell RNA-seq experiments with cDNA library
equalization. *Nucleic Acids Research*, *50*(2), e12–e12.
<https://doi.org/10.1093/nar/gkab1071>.

Hao, Y., Stuart, T. A., Kowalski, M. H., Choudhary, S., Hoffman, P.,
Hartman, A., Srivastava, A., Molla, G., Shaista Madad, Fernandez-Granda,
C., & Rahul atija. (2023). Dictionary learning for integrative,
multimodal and scalable single-cell analysis. Nature Biotechnology, 42,
293–304. <https://doi.org/10.1038/s41587-023-01767-y>

Lu, J., Sheng, Y., Qian, W., Pan, M., Zhao, X., & Ge, Q. (2023).
scRNA‐seq data analysis method to improve analysis performance. *Iet
Nanobiotechnology*, *17*(3), 246–256.
<https://doi.org/10.1049/nbt2.12115>.

McCarthy DJ, Campbell KR, Lun ATL, Willis QF (2017). “Scater:
pre-processing, quality control, normalisation and visualisation of
single-cell RNA-seq data in R.” \_Bioinformatics\_, \*33\*, 1179-1186.
<doi:10.1093/bioinformatics/btw777>
<https://doi.org/10.1093/bioinformatics/btw777>.

OpenAI. (2025). *ChatGPT*. ChatGPT 5; OpenAI. <https://chatgpt.com/>.

R Core Team (2025). \_R: A Language and Environment for Statistical
Computing\_. R Foundation for Statistical Computing, Vienna, Austria.
<https://www.R-project.org/>.

Wickham H (2016). *ggplot2: Elegant Graphics for Data Analysis*.
Springer-Verlag New York. ISBN 978-3-319-24277-4,
[https://ggplot2.tidyverse.org](https://ggplot2.tidyverse.org/).

Ziebell F (2024). \_RNAseqQC: Quality Control for RNA-Seq Data\_.
<doi:10.32614/CRAN.package.RNAseqQC>
<https://doi.org/10.32614/CRAN.package.RNAseqQC>, Rpackage version
0.2.1, <https://CRAN.R-project.org/package=RNAseqQC>.

## Acknowledgements

This package was developed as part of an assessment for 2025 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `QCseq` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the [GitHub
issues](https://github.com/anjalisilva/TestingPackage/issues).
