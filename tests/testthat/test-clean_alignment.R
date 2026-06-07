## a small alignment where seq1 and seq3 are identical (a repeated haplotype)
make_alignment <- function() {
  seq1 <- c("a", "t", "c", "g", "a", "t", "a", "c", "c", "g", "a")
  seq2 <- c("a", "c", "g", "g", "a", "g", "g", "t", "t", "g", "c")
  seq3 <- c("a", "t", "c", "g", "a", "t", "a", "c", "c", "g", "a")
  ape::as.DNAbin(rbind(seq1, seq2, seq3))
}

test_that("clean_alignment finds repeated haplotypes and de-duplicates", {
  res <- clean_alignment(make_alignment())
  expect_named(res, c("synonymies", "kept", "dropped", "alignment_clean"))
  expect_s3_class(res$alignment_clean, "DNAbin")

  ## seq1 and seq3 form one group; default keeps the alphabetically first (seq1)
  expect_equal(sort(res$synonymies$lineage), c("seq1", "seq3"))
  expect_identical(res$kept, "seq1")
  expect_identical(res$dropped, "seq3")
  expect_equal(nrow(res$alignment_clean), 2)
})

test_that("clean_alignment lets the user choose which lineage to keep", {
  res <- clean_alignment(make_alignment(), keep = "seq3")
  expect_identical(res$kept, "seq3")
  expect_identical(res$dropped, "seq1")
})

test_that("clean_alignment select = 'random' keeps a group member, reproducibly", {
  ## the only repeated group is {seq1, seq3}; a random pick must be one of them
  set.seed(42)
  res <- clean_alignment(make_alignment(), select = "random")
  expect_length(res$kept, 1)
  expect_true(res$kept %in% c("seq1", "seq3"))
  expect_setequal(c(res$kept, res$dropped), c("seq1", "seq3"))

  ## same seed -> same choice
  set.seed(42)
  res2 <- clean_alignment(make_alignment(), select = "random")
  expect_identical(res$kept, res2$kept)

  ## keep still overrides the random rule
  res3 <- clean_alignment(make_alignment(), select = "random", keep = "seq3")
  expect_identical(res3$kept, "seq3")
})

test_that("clean_alignment errors when nothing repeats or input is wrong", {
  res <- clean_alignment(make_alignment())
  expect_error(clean_alignment(res$alignment_clean), "no repeated haplotypes")
  expect_error(clean_alignment(matrix(1:4)), "DNAbin")
})

test_that("overlap method collapses a partial into a more complete sequence", {
  full    <- c("a", "t", "c", "g", "a", "t", "a", "c", "c", "g", "a")
  other   <- c("a", "c", "g", "g", "a", "g", "g", "t", "t", "g", "c")
  partial <- c("a", "t", "c", "g", "-", "-", "-", "c", "c", "g", "a")  # subset of full
  aln <- ape::as.DNAbin(rbind(full = full, other = other, partial = partial))

  ## strict: partial differs from full (gaps) -> no repeats
  expect_error(clean_alignment(aln, method = "strict"), "no repeated haplotypes")

  ## overlap: partial is contained in full -> grouped, full (more complete) kept
  res <- clean_alignment(aln, method = "overlap")
  expect_identical(res$kept, "full")
  expect_identical(res$dropped, "partial")
  expect_equal(nrow(res$alignment_clean), 2)
})
