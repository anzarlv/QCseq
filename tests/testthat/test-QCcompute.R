library(QCseq)

test_that("QCcompute runs correctly on example data", {
  data("sample_scRNAseq")

  qc_results <- QCcompute(sample_scRNAseq)

  # Basic checks on output
  testthat::expect_s3_class(qc_results, "data.frame")
  testthat::expect_true(nrow(qc_results) > 0)
  testthat::expect_true(all(c("Cell", "nGenes", "nTranscipts",
                              "percent_mit", "percent_ribo") %in%
                              colnames(qc_results)))

  # Values should be numeric and non-negative
  testthat::expect_true(all(qc_results$nGenes >= 0))
  testthat::expect_true(all(qc_results$nTranscipts >= 0))
  testthat::expect_true(all(qc_results$percent_mit >= 0, na.rm = TRUE))
  testthat::expect_true(all(qc_results$percent_ribo >= 0, na.rm = TRUE))
})


test_that("QCcompute handles invalid input", {
  # Non-Seurat object should throw an error
  fake_matrix <- matrix(rpois(25, 5), nrow = 5)
  testthat::expect_error(QCcompute(fake_matrix))

  # Missing input should throw an error
  testthat::expect_error(QCcompute())
})

# [END]
