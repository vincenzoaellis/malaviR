## A tiny in-frame reference alignment (length 12 = 4 codons) so translation,
## site profiling, and nearest-lineage logic can be checked deterministically.
##   ref1: ATG TTT GGG CCC  -> M F G P
##   ref2: ATG TTC GGG CCC  -> M F G P (synonymous variant at position 6)
##   ref3: ATG TTT GGA CCC  -> position-9 variant
make_ref <- function() {
  ref1 <- strsplit("atgtttgggccc", "")[[1]]
  ref2 <- strsplit("atgttcgggccc", "")[[1]]
  ref3 <- strsplit("atgtttggaccc", "")[[1]]
  ape::as.DNAbin(rbind(ref1 = ref1, ref2 = ref2, ref3 = ref3))
}

test_that("build_malavi_site_profile summarizes sites correctly", {
  prof <- build_malavi_site_profile(make_ref())
  expect_equal(nrow(prof), 12)
  expect_equal(prof$codon_position, rep(1:3, 4))

  ## position 1 is invariant (all A); position 6 varies (T in ref1/ref3, C in ref2)
  expect_true(prof$invariant[1])
  expect_equal(prof$major_base[1], "A")
  expect_false(prof$invariant[6])
  expect_equal(sort(strsplit(prof$observed_alleles[6], "")[[1]]), c("C", "T"))
})

test_that("lineage_qc calls an exact match a known lineage", {
  qc <- lineage_qc("atgtttgggccc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_s3_class(qc, "malavi_lineage_qc")
  expect_equal(qc$call, "known_lineage")
  expect_equal(qc$summary$call, "known_lineage")
  expect_equal(qc$summary$nearest_distance, 0)
  expect_true("exact_match_to_known_lineage" %in% qc$flags)
})

test_that("lineage_qc reports distance to the nearest lineage", {
  ## one base off ref1 (position 12 C->A) -> nearest distance 1
  qc <- lineage_qc("atgtttgggcca", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_equal(qc$summary$nearest_distance, 1)
  expect_true("near_known_lineage" %in% qc$flags)
})

test_that("lineage_qc detects a stop codon under genetic code 4", {
  ## TAA at codon 2 is a stop under code 4
  qc <- lineage_qc("atgtaagggccc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_true("contains_stop_codon" %in% qc$flags)
  expect_equal(qc$summary$n_stop_codons, 1)
  expect_equal(qc$call, "invalid_or_strong_warning")

  ## TGA is tryptophan (not a stop) under code 4 -> no stop flag
  qc2 <- lineage_qc("atgtgagggccc", make_ref(),
                    expected_length = 12, chimera_check = FALSE)
  expect_false("contains_stop_codon" %in% qc2$flags)
})

test_that("lineage_qc returns an invalid_sequence result for wrong length", {
  qc <- lineage_qc("atgtttgggcc", make_ref(),    # 11 bp
                   expected_length = 12, chimera_check = FALSE)
  expect_equal(qc$call, "invalid_sequence")
  expect_equal(qc$score, 0)
  expect_true(any(grepl("wrong_length", qc$flags)))
})

test_that("lineage_qc honors a user-set rare_base_frequency threshold", {
  ## a base that is real but uncommon at its site should be counted as "rare"
  ## only once the rare_base_frequency cutoff is raised above its frequency
  q <- "atgttcgggccc"   # the position-6 minority base (C)
  n_low  <- lineage_qc(q, make_ref(), expected_length = 12,
                       rare_base_frequency = 0.001, chimera_check = FALSE)$counts[["n_rare_site_bases"]]
  n_high <- lineage_qc(q, make_ref(), expected_length = 12,
                       rare_base_frequency = 0.99, chimera_check = FALSE)$counts[["n_rare_site_bases"]]
  expect_gt(n_high, n_low)
})

test_that("lineage_qc keeps full category counts and details on request", {
  qc <- lineage_qc("atgtttgggcca", make_ref(), expected_length = 12,
                   chimera_check = FALSE)
  expect_true(all(c("n_invariant_site_changes", "n_rare_site_bases",
                    "n_transversions") %in% names(qc$counts)))
  expect_null(qc$translation)            # details off by default

  qc2 <- lineage_qc("atgtttgggcca", make_ref(), expected_length = 12,
                    chimera_check = FALSE, details = TRUE)
  expect_false(is.null(qc2$translation))
})

test_that("an exact match outranks an N-containing near-twin at distance 0", {
  ## Regression for the .qc_nearest tie-break bug: distance is computed with
  ## pairwise deletion (a reference N is skipped, not a mismatch), so a reference
  ## that is identical to the query except for one N ties a true exact match at
  ## distance 0. The N-twin is placed FIRST here so that the old arbitrary
  ## tie-break (alignment order) would wrongly report it; the fix breaks ties by
  ## most comparable positions, so the genuine full-length exact match wins.
  ## (Mirrors the real case: a perfect P_SEIAUR01 ASV was mislabeled P_CARCAR11,
  ## whose MalAvi reference carries a single N.)
  ntwin <- strsplit("atgtttgggccn", "")[[1]]   # identical to query but N at pos 12
  clean <- strsplit("atgtttgggccc", "")[[1]]   # the genuine exact match
  ref   <- ape::as.DNAbin(rbind(ntwin = ntwin, clean = clean))

  qc <- lineage_qc("atgtttgggccc", ref, expected_length = 12, chimera_check = FALSE)
  expect_equal(qc$summary$nearest_distance, 0)
  expect_true("exact_match_to_known_lineage" %in% qc$flags)
  expect_equal(qc$summary$nearest_lineage, "clean")   # not the N-twin
})

test_that("lineage_qc rejects unsupported genetic codes", {
  expect_error(lineage_qc("atgtttgggccc", make_ref(), genetic_code = 1),
               "genetic_code = 4")
})
