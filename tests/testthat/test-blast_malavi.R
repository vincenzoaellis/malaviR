has_decipher <- requireNamespace("DECIPHER", quietly = TRUE) &&
  utils::packageVersion("DECIPHER") >= "3.0.0" &&
  requireNamespace("Biostrings", quietly = TRUE)

test_that("blast_malavi errors helpfully when DECIPHER >= 3.0 is unavailable", {
  skip_if(has_decipher)
  expect_error(blast_malavi("ACGTACGTACGT"), "DECIPHER")
})

test_that("blast_malavi finds an exact MalAvi sequence as its top hit", {
  skip_if_not(has_decipher, "needs DECIPHER >= 3.0 and Biostrings")

  ## use a real MalAvi sequence (gaps removed) as the query -> self-match at 100%
  aln <- extract_alignment()
  query <- gsub("-", "", paste(as.character(aln[1, ]), collapse = ""))

  res <- blast_malavi(query, top_n = 3)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 3)
  expect_true(all(c("Lineage", "ProportionMatch", "PercentMatch", "AlignmentLength",
                    "Matches", "Mismatches", "Score", "QueryGapLength",
                    "ReferenceLineageLength", "ReferenceFullLength") %in% names(res)))
  expect_identical(res$Lineage[1], rownames(aln)[1])
  expect_equal(res$PercentMatch[1], 100)
  expect_true(res$ReferenceFullLength[1] > 0)
})

test_that("blast_malavi warns and returns NA row when there are no hits", {
  skip_if_not(has_decipher, "needs DECIPHER >= 3.0 and Biostrings")
  expect_warning(out <- blast_malavi(paste(rep("A", 80), collapse = "")),
                 "No hits found")
  expect_s3_class(out, "data.frame")
})

test_that("blast_malavi rejects invalid query characters", {
  skip_if_not(has_decipher, "needs DECIPHER >= 3.0 and Biostrings")
  expect_error(blast_malavi("not a sequence!"), "invalid DNA")
})
