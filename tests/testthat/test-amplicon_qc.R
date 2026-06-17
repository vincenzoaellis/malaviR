## Reuse the tiny in-frame reference from the lineage_qc tests.
make_ref_amp <- function() {
  ref1 <- strsplit("atgtttgggccc", "")[[1]]
  ref2 <- strsplit("atgttcgggccc", "")[[1]]
  ref3 <- strsplit("atgtttggaccc", "")[[1]]
  ape::as.DNAbin(rbind(ref1 = ref1, ref2 = ref2, ref3 = ref3))
}

small_lineage_thresholds <- function() {
  th <- default_lineage_qc_thresholds()
  th$expected_length <- 12
  th
}

test_that("amplicon_qc returns the expected structure and per-variant calls", {
  ## two abundant known lineages plus a rare one-base derivative of the first
  variants <- data.frame(
    sequence = c("atgtttgggccc", "atgttcgggccc", "atgtttgggcca"),
    count    = c(10000, 4000, 5),
    stringsAsFactors = FALSE
  )
  aqc <- amplicon_qc(variants, make_ref_amp(),
                     lineage_qc_thresholds = small_lineage_thresholds(),
                     chimera_check = FALSE)

  expect_s3_class(aqc, "malavi_amplicon_qc")
  expect_true(all(c("relative_frequency", "lineage_call", "amplicon_call",
                    "amplicon_flags", "nearest_malavi_lineage") %in% names(aqc)))
  expect_equal(nrow(aqc), 3)

  ## the rare one-off derivative should be flagged as a likely amplicon artifact
  rare <- aqc[aqc$count == 5, ]
  expect_true(grepl("oneoff_from_much_more_abundant_variant", rare$amplicon_flags))
  expect_equal(rare$amplicon_call, "possible_amplicon_artifact")

  ## the abundant known lineages should pass
  expect_true(all(aqc$amplicon_call[aqc$count > 1000] == "passes"))
})

test_that("amplicon_qc errors on missing columns", {
  expect_error(
    amplicon_qc(data.frame(seq = "atgtttgggccc", count = 1),
                make_ref_amp(), chimera_check = FALSE),
    "sequence_col"
  )
  expect_error(
    amplicon_qc(data.frame(sequence = "atgtttgggccc", n = 1),
                make_ref_amp(), chimera_check = FALSE),
    "count_col"
  )
})

test_that("amplicon_qc computes relative frequencies within samples", {
  variants <- data.frame(
    sequence = c("atgtttgggccc", "atgttcgggccc",
                 "atgtttgggccc", "atgtttggaccc"),
    count    = c(75, 25, 90, 10),
    sample   = c("A", "A", "B", "B"),
    stringsAsFactors = FALSE
  )
  aqc <- amplicon_qc(variants, make_ref_amp(), sample_col = "sample",
                     lineage_qc_thresholds = small_lineage_thresholds(),
                     chimera_check = FALSE)
  ## frequencies are within-sample: each sample's relative_frequency sums to 1
  by_sample <- tapply(aqc$relative_frequency, aqc$sample, sum)
  expect_equal(as.numeric(by_sample), c(1, 1))
})
