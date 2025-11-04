library(QCseq)

test_that("QCsummary runs correctly and returns a character summary", {
  data("sample_scRNAseq")

  qc_summary <- QCsummary(sample_scRNAseq)

  # Should return a single string
  testthat::expect_type(qc_summary, "character")
  testthat::expect_length(qc_summary, 1)

  # Should mention all key QC metrics
  expected_terms <- c("nGenes", "nTranscipts", "percent_mit",
                      "percent_ribo", "qc_score_rank")
  for (term in expected_terms) {
    testthat::expect_true(any(grepl(term, qc_summary)),
                          info = paste("Missing metric:", term))
  }

  # Should include percentage and quality label
  testthat::expect_true(grepl("% of cells", qc_summary))
  testthat::expect_true(grepl("quality profile", qc_summary))
})

# ChatGPT (OpenAI, 2025) helped with this test below
test_that("QCsummary reflects changes in threshold", {
  data("sample_scRNAseq")

  text_low_thresh  <- QCsummary(sample_scRNAseq, threshold = 0.2)
  text_high_thresh <- QCsummary(sample_scRNAseq, threshold = 0.8)

  # The reported percentage of low-quality cells should differ
  percent_pattern <- "(\\d+\\.?\\d*)% of cells"
  pct_low  <- as.numeric(sub(percent_pattern, "\\1", regmatches(text_low_thresh,
                                                                regexpr(percent_pattern, text_low_thresh))))
  pct_high <- as.numeric(sub(percent_pattern, "\\1", regmatches(text_high_thresh,
                                                                regexpr(percent_pattern, text_high_thresh))))

  testthat::expect_true(is.finite(pct_low))
  testthat::expect_true(is.finite(pct_high))
  testthat::expect_true(pct_high != pct_low)
})

test_that("QCsummary handles invalid input", {
  fake_matrix <- matrix(rpois(25, 5), nrow = 5)

  testthat::expect_error(QCsummary(fake_matrix))
  testthat::expect_error(QCsummary())
})

# [END]

