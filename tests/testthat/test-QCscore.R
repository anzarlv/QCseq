library(QCseq)

test_that("QCscore runs correctly on example data", {
  data("sample_scRNAseq")

  qc_results <- QCscore(sample_scRNAseq)

  # Output checks
  testthat::expect_s3_class(qc_results, "data.frame")
  testthat::expect_true(all(c("Cell", "n_genes", "n_transcripts",
                              "percent_mit", "percent_ribo",
                              "qc_score_rank") %in%
                              colnames(qc_results)))

  # QC score range should be between 0 and 1
  testthat::expect_true(all(qc_results$qc_score_rank >= 0, na.rm = TRUE))
  testthat::expect_true(all(qc_results$qc_score_rank <= 1, na.rm = TRUE))

  # QC score should correlate positively with n_genes (rough sanity check)
  cor_val <- cor(qc_results$n_genes, qc_results$qc_score_rank, use =
                   "complete.obs")
  testthat::expect_true(!is.na(cor_val))
  testthat::expect_true(cor_val > 0)
})

test_that("QCscore handles invalid input", {
  # Non-Seurat object
  fake_matrix <- matrix(rpois(25, 5), nrow = 5)
  testthat::expect_error(QCscore(fake_matrix))

  # Missing input
  testthat::expect_error(QCscore())
})

# [END]
