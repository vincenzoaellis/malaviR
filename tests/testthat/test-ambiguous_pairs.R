test_that("ambiguous_pairs returns the documented structure", {
  ap <- ambiguous_pairs()
  expect_named(ap, c("summary", "by_genus", "pairs"))
  expect_s3_class(ap$summary, "data.frame")
  expect_s3_class(ap$pairs, "data.frame")
  expect_named(ap$pairs,
    c("lineage_a", "lineage_b", "genus_a", "genus_b", "same_genus",
      "n_comparable", "a_private", "b_private"))
  expect_equal(nrow(ap$summary), 1L)
  expect_equal(ap$summary$n_ambiguous_pairs, nrow(ap$pairs))
})

test_that("a synthetic alignment isolates the mutual-partial pair only", {
  ## A: fully determined.
  ## B: A with an N at the last position  -> contained in A (synonymy, excluded).
  ## C: determined at the last position but N at the first.
  ## E: a real one-base neighbor of A at the last position (distance 1, excluded).
  ## B and C are mutually partial and agree on their overlap -> the only pair.
  aln <- c(
    A = "ACGTACGTAC",
    B = "ACGTACGTAN",
    C = "NCGTACGTAC",
    E = "ACGTACGTAG"
  )
  ap <- ambiguous_pairs(alignment = aln)

  expect_equal(ap$summary$n_ambiguous_pairs, 1L)
  got <- sort(c(ap$pairs$lineage_a, ap$pairs$lineage_b))
  expect_equal(got, c("B", "C"))
  expect_equal(ap$pairs$n_comparable, 8L)   # positions 2:9 determined in both
  expect_equal(ap$pairs$a_private, 1L)
  expect_equal(ap$pairs$b_private, 1L)
})

test_that("containments are excluded (they are synonymy, not ambiguous pairs)", {
  ## D is strictly contained in A: identical wherever D is determined, and A has
  ## no ambiguous position for D to be privately determined at.
  aln <- c(A = "ACGTACGTAC", D = "ACGTACGTNN")
  ap <- ambiguous_pairs(alignment = aln)
  expect_equal(ap$summary$n_ambiguous_pairs, 0L)
})

test_that("a fully-determined alignment has no ambiguous pairs", {
  aln <- c(A = "ACGTACGTAC", B = "ACGTACGTAG", C = "TCGTACGTAC")
  ap <- ambiguous_pairs(alignment = aln)
  expect_equal(ap$summary$n_partial, 0L)
  expect_equal(ap$summary$n_ambiguous_pairs, 0L)
})

test_that("true one-base neighbors are not reported", {
  ## both partial, but they conflict at a position determined in both (pos 5):
  ## distance 1, so not an ambiguous pair despite the mutual ambiguity elsewhere.
  aln <- c(A = "ACGTACGTAN", B = "NCGTGCGTAC")
  ap <- ambiguous_pairs(alignment = aln)
  expect_equal(ap$summary$n_ambiguous_pairs, 0L)
})

test_that("non-overlapping mutually-partial sequences are excluded", {
  ## P determined only at pos1, Q only at pos10 -> they never overlap on a
  ## determined base (comparable = 0), so this is not a meaningful pair.
  aln <- c(P = "ANNNNNNNNN", Q = "NNNNNNNNNA")
  ap <- ambiguous_pairs(alignment = aln, min_comparable = 1L)
  expect_equal(ap$summary$n_ambiguous_pairs, 0L)
})

test_that("min_comparable filters out weakly-overlapping pairs", {
  ## P and Q agree on a single overlapping determined position (pos5) and are
  ## mutually partial elsewhere -> a pair at min_comparable = 1, gone at 2.
  aln <- c(P = "ANNNANNNNN", Q = "NNNNANNNNA")
  ap1 <- ambiguous_pairs(alignment = aln, min_comparable = 1L)
  ap2 <- ambiguous_pairs(alignment = aln, min_comparable = 2L)
  expect_equal(ap1$summary$n_ambiguous_pairs, 1L)
  expect_equal(ap1$pairs$n_comparable, 1L)
  expect_equal(ap2$summary$n_ambiguous_pairs, 0L)
})

test_that("known ambiguous pairs and excluded containments on the bundled MalAvi data", {
  skip_on_cran()
  ap <- ambiguous_pairs()

  ## a validated genuine ambiguous pair: agree across 477 positions, each with a
  ## single private determined base (see results/notes/malaviASV_plan.md)
  has_pair <- function(p, x, y) {
    any((grepl(x, p$lineage_a) & grepl(y, p$lineage_b)) |
        (grepl(y, p$lineage_a) & grepl(x, p$lineage_b)))
  }
  expect_true(has_pair(ap$pairs, "H_DENPEN02", "H_PASILI01"))

  ## the pairs the plan once cited as "ambiguous" are actually containments that
  ## synonymy_report collapses, so they must NOT appear here
  expect_false(has_pair(ap$pairs, "SEIAUR01", "CARCAR11"))
  expect_false(has_pair(ap$pairs, "SEIAUR01", "TABI08"))

  ## sanity: every reported pair really is distance 0 with both sides private
  expect_true(all(ap$pairs$a_private > 0))
  expect_true(all(ap$pairs$b_private > 0))
  expect_true(all(ap$pairs$n_comparable >= ap$summary$min_comparable))
})
