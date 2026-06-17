## Reuse the tiny in-frame reference from the lineage_qc tests.
make_ref_amp <- function() {
  ref1 <- strsplit("atgtttgggccc", "")[[1]]
  ref2 <- strsplit("atgttcgggccc", "")[[1]]
  ref3 <- strsplit("atgtttggaccc", "")[[1]]
  ape::as.DNAbin(rbind(ref1 = ref1, ref2 = ref2, ref3 = ref3))
}

test_that("amplicon_qc returns the expected structure and per-variant calls", {
  ## two abundant known lineages plus a rare one-base derivative of the first
  variants <- data.frame(
    sequence = c("atgtttgggccc", "atgttcgggccc", "atgtttgggcca"),
    count    = c(10000, 4000, 5),
    stringsAsFactors = FALSE
  )
  aqc <- amplicon_qc(variants, make_ref_amp(), chimera_check = FALSE)

  expect_s3_class(aqc, "malavi_amplicon_qc")
  expect_true(all(c("relative_frequency", "lineage_call", "amplicon_call",
                    "amplicon_flags", "nearest_malavi_lineage",
                    "n_nonsynonymous") %in% names(aqc)))
  expect_equal(nrow(aqc), 3)

  ## the rare one-off derivative should be flagged as a likely amplicon artifact
  rare <- aqc[aqc$count == 5, ]
  expect_true(grepl("oneoff_from_much_more_abundant_variant", rare$amplicon_flags))
  expect_equal(rare$amplicon_call, "possible_amplicon_artifact")

  ## the abundant known lineages should pass
  expect_true(all(aqc$amplicon_call[aqc$count > 1000] == "passes"))
})

test_that("amplicon_qc errors on a malformed two-column input", {
  ## fewer than two columns
  expect_error(
    amplicon_qc(data.frame(sequence = "atgtttgggccc"),
                make_ref_amp(), chimera_check = FALSE),
    "ASV"
  )
  ## columns in the wrong order/type (counts first, sequences second)
  expect_error(
    amplicon_qc(data.frame(count = 1, sequence = "atgtttgggccc"),
                make_ref_amp(), chimera_check = FALSE),
    "look wrong"
  )
})

test_that("amplicon_qc errors when a sequence is not the reference length", {
  variants <- data.frame(
    sequence = c("atgtttgggccc", "atgtttgggcc"),   # second sequence is 11 bp
    count    = c(100, 5),
    stringsAsFactors = FALSE
  )
  expect_error(
    amplicon_qc(variants, make_ref_amp(), chimera_check = FALSE),
    "reference alignment length"
  )
})

test_that("amplicon_qc computes relative frequencies within the pool", {
  variants <- data.frame(
    sequence = c("atgtttgggccc", "atgttcgggccc", "atgtttggaccc"),
    count    = c(75, 20, 5),
    stringsAsFactors = FALSE
  )
  aqc <- amplicon_qc(variants, make_ref_amp(), chimera_check = FALSE)
  ## relative frequencies are read fractions over the whole pool, summing to 1
  expect_equal(sum(aqc$relative_frequency), 1)
  expect_equal(nrow(aqc), 3)
})

test_that("amplicon_qc accepts a two-column frame by position, not by name", {
  ## column names are irrelevant; the first column is sequences, second counts
  variants <- data.frame(
    asv   = c("atgtttgggccc", "atgttcgggccc"),
    reads = c(100, 40),
    stringsAsFactors = FALSE
  )
  aqc <- amplicon_qc(variants, make_ref_amp(), chimera_check = FALSE)
  expect_equal(nrow(aqc), 2)
  expect_true(all(c("sequence", "count") %in% names(aqc)))
})
