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

## thresholds tuned to the tiny reference (expected length 12)
small_thresholds <- function() {
  th <- default_lineage_qc_thresholds()
  th$expected_length <- 12
  th
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
                   thresholds = small_thresholds(), chimera_check = FALSE)
  expect_s3_class(qc, "malavi_lineage_qc")
  expect_equal(qc$summary$call, "known_lineage")
  expect_equal(qc$summary$nearest_distance, 0)
  expect_true("exact_match_to_known_lineage" %in% qc$flags)
})

test_that("lineage_qc reports distance to the nearest lineage", {
  ## one base off ref1 (position 12 C->A) -> nearest distance 1
  qc <- lineage_qc("atgtttgggcca", make_ref(),
                   thresholds = small_thresholds(), chimera_check = FALSE)
  expect_equal(qc$summary$nearest_distance, 1)
  expect_true("near_known_lineage" %in% qc$flags)
})

test_that("lineage_qc detects a stop codon under genetic code 4", {
  ## TAA at codon 2 is a stop under code 4
  qc <- lineage_qc("atgtaagggccc", make_ref(),
                   thresholds = small_thresholds(), chimera_check = FALSE)
  expect_true("contains_stop_codon" %in% qc$flags)
  expect_equal(qc$summary$n_stop_codons, 1)
  expect_equal(qc$summary$call, "invalid_or_strong_warning")

  ## TGA is tryptophan (not a stop) under code 4 -> no stop flag
  qc2 <- lineage_qc("atgtgagggccc", make_ref(),
                    thresholds = small_thresholds(), chimera_check = FALSE)
  expect_false("contains_stop_codon" %in% qc2$flags)
})

test_that("lineage_qc returns an invalid_sequence result for wrong length", {
  qc <- lineage_qc("atgtttgggcc", make_ref(),    # 11 bp
                   thresholds = small_thresholds(), chimera_check = FALSE)
  expect_equal(qc$call, "invalid_sequence")
  expect_equal(qc$overall_score, 0)
  expect_true(any(grepl("wrong_length", qc$flags)))
})

test_that("lineage_qc honors a user-set rare_base_frequency threshold", {
  ## a base that is real but uncommon at its site should be counted as "rare"
  ## only once the rare_base_frequency cutoff is raised above its frequency
  th_low  <- small_thresholds(); th_low$rare_base_frequency  <- 0.001
  th_high <- small_thresholds(); th_high$rare_base_frequency <- 0.99
  q <- "atgttcgggccc"   # the position-6 minority base (C)
  n_low  <- lineage_qc(q, make_ref(), thresholds = th_low,  chimera_check = FALSE)$summary$n_rare_site_bases
  n_high <- lineage_qc(q, make_ref(), thresholds = th_high, chimera_check = FALSE)$summary$n_rare_site_bases
  expect_gt(n_high, n_low)
})

test_that("lineage_qc rejects unsupported genetic codes", {
  expect_error(lineage_qc("atgtttgggccc", make_ref(), genetic_code = 1),
               "genetic_code = 4")
})
