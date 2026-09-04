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
                    "ReferenceGapLength", "ReferenceLineageLength",
                    "ReferenceFullLength") %in% names(res)))
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

test_that("top_n = 1 returns the lineage the query came from", {
  skip_if_not(has_decipher, "needs DECIPHER >= 3.0 and Biostrings")

  ## The regression this guards against. Until version 1.1.2 top_n was applied to
  ## the SearchIndex k-mer score before anything was aligned, and that score does
  ## not rank an exact match first: blast_malavi(<full SGS1>, top_n = 1) returned
  ## P_PADOM07 at 99.776% and never reported the 100% self-match. The old self-hit
  ## test passed only because it used top_n = 3 on the first sequence in the file.
  ## Several lineages, because the failure depended on which references the index
  ## happened to score above the true one.
  aln  <- extract_alignment()
  seqs <- toupper(vapply(seq_len(nrow(aln)),
                         function(i) paste(as.character(aln[i, ]), collapse = ""),
                         character(1)))
  ## complete, unambiguous sequences only: a partial one is identical to several
  ## other lineages over its shorter window, so "the" top hit is genuinely a tie
  complete <- which(!grepl("[^ACGT]", seqs) & nchar(seqs) == 479L)
  set.seed(42)
  for (i in sample(complete, 4)) {
    res <- blast_malavi(seqs[i], top_n = 1)
    expect_identical(res$Lineage[1], rownames(aln)[i])
    expect_equal(res$PercentMatch[1], 100)
  }
})

test_that("a query with an insertion still ranks its source lineage first, and the columns add up", {
  skip_if_not(has_decipher, "needs DECIPHER >= 3.0 and Biostrings")

  ## Two things at once. The insertion pushes the source lineage down the index
  ## ranking (it was omitted entirely at top_n = 2), and it opens gaps in the
  ## REFERENCE, whose lengths used to be dropped from the result -- so the columns
  ## did not reconcile and the user could not see where the extra alignment
  ## columns came from.
  aln  <- extract_alignment()
  seqs <- toupper(vapply(seq_len(nrow(aln)),
                         function(i) paste(as.character(aln[i, ]), collapse = ""),
                         character(1)))
  i    <- which(!grepl("[^ACGT]", seqs) & nchar(seqs) == 479L)[1]
  full <- seqs[i]
  query <- paste0(substr(full, 1, 240), "ACG", substr(full, 241, 479))

  res <- blast_malavi(query, top_n = 2)
  expect_identical(res$Lineage[1], rownames(aln)[i])
  expect_true(all(res$Matches + res$Mismatches + res$QueryGapLength +
                    res$ReferenceGapLength == res$AlignmentLength))
  expect_gt(res$ReferenceGapLength[1], 0)
})
