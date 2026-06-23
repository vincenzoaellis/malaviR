test_that("pairwise deletion skips ambiguous positions in either sequence", {
  ref <- c(
    R1 = "ACGTACGTAC",   # identical to the query
    R2 = "ACGTACGTAG",   # one real difference, at position 10
    R3 = "ACGTACGTAN"    # identical to R1 except an N at position 10
  )
  q <- "ACGTACGTAC"
  res <- pairwise_deletion_distance(q, reference = ref)

  expect_named(res, c("lineage", "distance", "n_comparable"))
  ## ordered nearest first; among the two distance-0 hits the most complete
  ## (R1, 10 comparable) precedes the N-bearing R3 (9 comparable)
  expect_equal(res$lineage, c("R1", "R3", "R2"))
  expect_equal(res$distance, c(0, 0, 1))
  expect_equal(res$n_comparable, c(10, 9, 10))
})

test_that("an N in the query is also skipped, not counted as a difference", {
  ref <- c(R1 = "ACGTACGTAC")
  res <- pairwise_deletion_distance("ACGTACGTAN", reference = ref)
  expect_equal(res$distance, 0)
  expect_equal(res$n_comparable, 9)   # position 10 dropped (query N)
})

test_that("top_n limits the number of references returned, nearest first", {
  ref <- c(R1 = "ACGTACGTAC", R2 = "ACGTACGTAG", R3 = "ACGTACGTTT")
  res <- pairwise_deletion_distance("ACGTACGTAC", reference = ref, top_n = 1L)
  expect_equal(nrow(res), 1L)
  expect_equal(res$lineage, "R1")
  expect_equal(res$distance, 0)
})

test_that("a query not aligned to the reference length errors clearly", {
  ref <- c(R1 = "ACGTACGTAC")
  expect_error(pairwise_deletion_distance("ACGT", reference = ref),
               "4 bp but the reference alignment is 10 bp")
})

test_that("query must be a single sequence", {
  ref <- c(R1 = "ACGTACGTAC")
  expect_error(pairwise_deletion_distance(c("ACGTACGTAC", "ACGTACGTAG"),
                                          reference = ref),
               "single aligned sequence")
})

test_that("against the bundled MalAvi data a lineage matches itself at distance 0", {
  skip_on_cran()
  aln <- extract_alignment()
  q <- paste(as.character(aln[1, ]), collapse = "")
  res <- pairwise_deletion_distance(q)
  expect_equal(nrow(res), nrow(aln))          # all references by default
  expect_equal(res$distance[1], 0)            # the nearest is an exact match
  expect_true(res$n_comparable[1] > 0)
})
