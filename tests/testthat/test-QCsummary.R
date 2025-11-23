library(QCseq)

test_that("QCsummary runs correctly and returns a character summary", {
  data("sample_scRNAseq")

  qc_summary <- QCsummary(sample_scRNAseq)

  # Should return a single string
  testthat::expect_type(qc_summary, "character")
  testthat::expect_length(qc_summary, 1)

  # Should mention all key QC metrics
  expected_terms <- c("n_genes", "n_transcripts", "percent_mit",
                      "percent_ribo", "qc_score_rank")
  for (term in expected_terms) {
    testthat::expect_true(any(grepl(term, qc_summary)),
                          info = paste("Missing metric:", term))
  }

  # Should include quality label
  testthat::expect_true(grepl("quality profile", qc_summary))
})

# ChatGPT (OpenAI, 2025) helped with this test below
test_that("QCsummary reflects changes in threshold", {
  data("sample_scRNAseq")

  # Use extreme thresholds so the quality label is very likely to change
  text_low_thresh  <- QCsummary(sample_scRNAseq,
                                lower_threshold = 0,
                                upper_threshold = 0)
  text_high_thresh <- QCsummary(sample_scRNAseq,
                                lower_threshold = 1,
                                upper_threshold = 1)

  # Extract quality labels from the sentence
  pattern <- "dataset has a ([a-z]+) quality profile"
  label_low  <- sub(pattern, "\\1", regmatches(text_low_thresh,
                                               regexpr(pattern, text_low_thresh)))
  label_high <- sub(pattern, "\\1", regmatches(text_high_thresh,
                                               regexpr(pattern, text_high_thresh)))

  testthat::expect_true(label_low %in% c("good", "moderate", "poor"))
  testthat::expect_true(label_high %in% c("good", "moderate", "poor"))
  testthat::expect_true(label_low != label_high)
})

test_that("QCsummary handles invalid input", {
  fake_matrix <- matrix(rpois(25, 5), nrow = 5)

  testthat::expect_error(QCsummary(fake_matrix))
  testthat::expect_error(QCsummary())
})

# [END]
