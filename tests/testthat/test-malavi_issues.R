## The property the whole design exists for: every sentence malavi_issues() prints is
## written from what its check found in the loaded release, so an issue a release fixes
## drops out by itself. These tests therefore derive their expectations from the same
## release rather than remembering what the 2026-03-23 release showed -- the earlier
## version of this file did remember (SETCAE02, PLOCUC14-16, seven N/A-genus lineages)
## and every one of those went stale the day release 2026-09-15 fixed them.

test_that("malavi_issues prints a heading naming the release, and nothing else", {
  out <- capture.output(malavi_issues())
  expect_equal(out[1],
               paste0("Known issues in the current MalAvi data release; ",
                      "MalAvi version ", malavi_version()))
  expect_equal(out[2], "")
})

test_that("the result is one row per issue the release actually shows", {
  capture.output(res <- malavi_issues())
  expect_s3_class(res, "data.frame")
  expect_named(res, c("title", "text"))

  registry <- malaviR:::.malavi_issue_registry()
  ctx      <- malaviR:::.malavi_issue_context("latest")
  found    <- vapply(registry,
                     function(i) length(unique(as.character(i$check(ctx)))) > 0,
                     logical(1))
  titles   <- vapply(registry, function(i) i$title, character(1))

  ## exactly the issues whose check found something, in registry order
  expect_equal(res$title, titles[found])
  expect_equal(nrow(res), sum(found))
})

test_that("the text is derived from the release, not stored", {
  registry <- malaviR:::.malavi_issue_registry()
  ctx      <- malaviR:::.malavi_issue_context("latest")
  capture.output(res <- malavi_issues())

  for (issue in registry) {
    affected <- unique(as.character(issue$check(ctx)))
    if (length(affected) == 0) {
      ## an issue the release no longer shows contributes no text at all
      expect_false(issue$title %in% res$title)
      next
    }
    txt <- res$text[res$title == issue$title]
    expect_length(txt, 1)
    ## the sentence names what the check found; the first affected lineage is the
    ## cheapest thing to assert without re-implementing the describe() prose
    first <- clean_names(affected[1])
    expect_true(grepl(first, txt, fixed = TRUE),
                info = paste(issue$title, "should name", first))
  }
})

test_that("the morphospecies issue reads its species from the table", {
  ## TUPHI01 has been two morphospecies since long before this package; if a future
  ## release resolves it, this test skips rather than fails
  gls <- extract_table("Grand Lineage Summary")
  dup <- sort(unique(gls$LINEAGE_NAME[duplicated(gls$LINEAGE_NAME)]))
  skip_if(length(dup) == 0, "no lineage is listed as more than one morphospecies")
  capture.output(res <- malavi_issues())
  txt <- res$text[res$title == "One lineage listed as more than one morphospecies"]
  for (lin in dup) {
    sp <- unique(gls$SPECIES_NAME[gls$LINEAGE_NAME == lin])
    sp <- sp[!is.na(sp) & nzchar(sp)]
    expect_match(txt, paste0("Lineage ", lin, " is listed as "))
    for (s in sp) {
      ## the second and later species are abbreviated to the genus initial
      expect_true(grepl(s, txt, fixed = TRUE) ||
                    grepl(sub("^(\\w)\\w+ ", "\\1. ", s), txt, fixed = TRUE),
                  info = paste(lin, "should name", s))
    }
  }
})

test_that("the counts in the prose match the data, not a remembered number", {
  gls      <- extract_table("Grand Lineage Summary")
  n_no_gen <- length(unique(gls$LINEAGE_NAME[!is.na(gls$GENUS_NAME) &
                                               gls$GENUS_NAME == "N/A"]))
  capture.output(res <- malavi_issues())
  txt <- res$text[res$title == "No parasite genus listed"]
  if (n_no_gen == 0) {
    expect_length(txt, 0)
  } else {
    expect_match(txt, paste0("^", n_no_gen, " lineage"))
  }
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

test_that("every registry entry has a check, a describe and a title", {
  registry <- malaviR:::.malavi_issue_registry()
  expect_true(all(vapply(registry, function(i) is.function(i$check), logical(1))))
  expect_true(all(vapply(registry, function(i) is.function(i$describe), logical(1))))
  expect_true(all(vapply(registry, function(i) nzchar(i$title), logical(1))))
})
