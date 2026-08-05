test_that("malavi_issues prints a heading naming the release, and nothing else", {
  out <- capture.output(malavi_issues())
  expect_equal(out[1],
               paste0("Known issues in the current MalAvi data release; ",
                      "MalAvi version ", malavi_version()))
  expect_equal(out[2], "")
})

test_that("the four issues recorded for the bundled release are still found", {
  capture.output(res <- malavi_issues())
  expect_s3_class(res, "data.frame")
  expect_named(res, c("title", "text"))
  expect_equal(nrow(res), 4)
})

test_that("the text is derived from the release, not stored", {
  ## This is the property the whole design exists for: the sentence must be
  ## written from what the check found in the loaded release, so it cannot go
  ## stale the way hard-coded text would.
  capture.output(res <- malavi_issues())
  txt <- paste(res$text, collapse = " ")

  expect_match(txt, "L_MEAPI23 has the Leucocytozoon prefix \\(L_\\) in the alignment")
  expect_match(txt, "but in the Grand Lineage Summary table it is listed as Plasmodium")
  expect_match(txt, "7 lineages in the Grand Lineage Summary table do not have")
  expect_match(txt, "ACCNIS06, ACCNIS07, ACCNIS08, ACCNIS09, BUBT1, BUTLIV02, CIAE08")
  expect_match(txt, "Lineage TUPHI01 is listed as two morphospecies")
  expect_match(txt, "Haemoproteus minutus and H\\. asymmetricus")
  expect_match(txt, "alignment but not the Grand Lineage Summary Table: SETCAE02")
  expect_match(txt, "not the alignment: PLOCUC14, PLOCUC15, PLOCUC16, SETCAE01")
})

test_that("the counts in the prose match the data, not a remembered number", {
  ## the "7 lineages" in the text must be a real count of "N/A" genus rows
  gls      <- extract_table("Grand Lineage Summary")
  n_no_gen <- length(unique(gls$LINEAGE_NAME[!is.na(gls$GENUS_NAME) &
                                               gls$GENUS_NAME == "N/A"]))
  capture.output(res <- malavi_issues())
  expect_match(res$text[2], paste0("^", n_no_gen, " lineages"))
})

test_that("the report carries no registry bookkeeping", {
  txt <- paste(capture.output(malavi_issues()), collapse = "\n")

  ## a list to read, not a table to program against
  expect_false(grepl("PRESENT", txt))
  expect_false(grepl("not_checked", txt))
  expect_false(grepl("reported by", txt))
  expect_false(grepl("recorded in", txt))
  expect_false(grepl("n_affected", txt))
  expect_false(grepl("Suggested workaround", txt))
  expect_false(grepl("extract_table\\(", txt))
})

test_that("the retired issues stay retired", {
  ## Sequence-name shape and cross-genus ambiguous pairs were dropped on
  ## 2026-08-05: the first is a well-known convention rather than a fault, and
  ## the second reported weak evidence alongside strong.
  txt <- paste(capture.output(malavi_issues()), collapse = "\n")
  expect_false(grepl("two shapes", txt))
  expect_false(grepl("ambiguous", txt))
})

test_that("an issue the release no longer shows drops out of the report", {
  ## Simulate a fixed release by checking the registry contract directly: a
  ## check returning nothing must contribute no text.
  registry <- malaviR:::.malavi_issue_registry()
  expect_true(all(vapply(registry, function(i) is.function(i$check), logical(1))))
  expect_true(all(vapply(registry, function(i) is.function(i$describe), logical(1))))
  expect_true(all(vapply(registry, function(i) nzchar(i$title), logical(1))))
})
