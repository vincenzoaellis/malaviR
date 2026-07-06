## The metadata attribute (attr(x, "malavi_meta")) and the S3 print methods added
## for reproducibility must not disturb the existing return structures: names(),
## element access, and identical() of the data must all be unchanged.

test_that("extract_table stamps the resolved MalAvi version without changing the data", {
  skip_on_cran()
  tab <- extract_table("Grand Lineage Summary")
  meta <- attr(tab, "malavi_meta")
  expect_equal(meta$malavi_version, malavi_version())
  ## the descriptive vs snake_case names must still give identical objects,
  ## i.e. the stamp does not encode the (differing) table label
  expect_identical(extract_table("references"),
                   extract_table("Table of References"))
})

test_that("bundled-data outputs carry a malavi_version stamp", {
  skip_on_cran()
  expect_equal(attr(extract_alignment(), "malavi_meta")$malavi_version,
               malavi_version())
  expect_equal(attr(lineage_studies(), "malavi_meta")$malavi_version,
               malavi_version())
})

test_that("a user-supplied reference is stamped with NA version, not the bundled one", {
  ref <- c(R1 = "ACGTACGTAC")
  res <- pairwise_deletion_distance("ACGTACGTAC", reference = ref)
  meta <- attr(res, "malavi_meta")
  expect_true(is.na(meta$malavi_version))
  expect_equal(meta$method, "pairwise_deletion")
})

test_that("lineage_qc records genetic code, expected length, and version", {
  skip_on_cran()
  aln <- extract_alignment()
  seq <- paste(as.character(aln[1, ]), collapse = "")
  qc <- lineage_qc(seq)
  meta <- attr(qc, "malavi_meta")
  expect_equal(meta$genetic_code, 4)
  expect_equal(meta$expected_length, 479)
  expect_equal(meta$malavi_version, malavi_version())
  ## print mentions the heuristic caveat
  expect_output(print(qc), "heuristic")
})

test_that("bundled-data list outputs record the MalAvi version (not NA)", {
  skip_on_cran()
  ## regression: these functions reassign their alignment/species argument
  ## internally, so the bundled-vs-supplied flag must be captured up front
  expect_equal(attr(synonymy_report(method = "strict"), "malavi_meta")$malavi_version,
               malavi_version())
  expect_equal(attr(match_taxonomy(), "malavi_meta")$malavi_version,
               malavi_version())
})

test_that("list outputs gain an S3 class and a working print method", {
  skip_on_cran()
  sr <- synonymy_report(method = "strict")
  expect_s3_class(sr, "malavi_synonymy_report")
  expect_named(sr, c("summary", "by_genus", "synonymies"))   # names unchanged
  expect_output(print(sr), "synonymy report")

  ap <- ambiguous_pairs()
  expect_s3_class(ap, "malavi_ambiguous_pairs")
  expect_named(ap, c("summary", "by_genus", "pairs"))
  expect_output(print(ap), "ambiguous pairs")

  tm <- match_taxonomy(c("Turdus merula", "Anas sp."))
  expect_s3_class(tm, "malavi_taxonomy_match")
  expect_named(tm, c("key", "differences"))
  expect_output(print(tm), "match_type")
})

test_that("clean_alignment output is classed and prints without disturbing names", {
  aln <- extract_alignment()
  res <- clean_alignment(aln, method = "strict")
  expect_s3_class(res, "malavi_alignment_clean")
  expect_named(res, c("synonymies", "kept", "dropped", "alignment_clean"))
  expect_output(print(res), "cleaned alignment")
})
