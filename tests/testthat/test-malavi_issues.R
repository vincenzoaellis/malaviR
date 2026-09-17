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

## ---- the genus-versus-neighbours check ---------------------------------------

## a hand-made alignment: `base` and three one-substitution relatives form a
## cluster; `far` and one relative differ from the cluster at every fourth
## position (about 25 % divergence, the order of a genus difference)
.genus_fixture <- function(seed = 1) {
  set.seed(seed)
  base <- sample(c("A", "C", "G", "T"), 479, replace = TRUE)
  mutate <- function(x, at) { x[at] <- ifelse(x[at] == "A", "C", "A"); x }
  far <- base
  far[seq(1, 479, by = 4)] <- ifelse(far[seq(1, 479, by = 4)] == "G", "T", "G")
  list(base = base, far = far, mutate = mutate)
}

test_that("the genus check finds a lineage planted inside a cluster of another genus", {
  f <- .genus_fixture()
  charmat <- rbind(P1 = f$base, P2 = f$mutate(f$base, 10), P3 = f$mutate(f$base, 20),
                   WRONG = f$mutate(f$base, 30), H1 = f$far, H2 = f$mutate(f$far, 40))
  genus <- c("Plasmodium", "Plasmodium", "Plasmodium", "Haemoproteus",
             "Haemoproteus", "Haemoproteus")
  out <- malaviR:::.malavi_genus_outliers(charmat, genus)
  expect_equal(out$lineage, "WRONG")
  expect_equal(out$genus, "Haemoproteus")
  expect_equal(out$neighbor_genus, "Plasmodium")
  expect_equal(out$n_neighbors, 3L)
  expect_equal(out$n_neighbors_other, 3L)
  expect_equal(out$nearest_mismatches, 1L)
  expect_equal(out$nearest_comparable, 479L)
  ## the two real Haemoproteus have only one neighbour each: below the minimum
  expect_false(any(c("H1", "H2") %in% out$lineage))
})

test_that("the genus check needs enough neighbours, enough overlap and a known genus", {
  f <- .genus_fixture(2)
  ## two neighbours only: not reported
  charmat <- rbind(P1 = f$base, P2 = f$mutate(f$base, 10), WRONG = f$mutate(f$base, 30))
  out <- malaviR:::.malavi_genus_outliers(charmat, c("Plasmodium", "Plasmodium", "Haemoproteus"))
  expect_equal(nrow(out), 0L)
  ## three neighbours, but the lineage is mostly undetermined: not reported
  short <- f$mutate(f$base, 30); short[1:200] <- "N"
  charmat <- rbind(P1 = f$base, P2 = f$mutate(f$base, 10), P3 = f$mutate(f$base, 20), WRONG = short)
  out <- malaviR:::.malavi_genus_outliers(charmat, c(rep("Plasmodium", 3), "Haemoproteus"))
  expect_equal(nrow(out), 0L)
  ## an undetermined stretch that still leaves enough overlap is counted right
  short <- f$mutate(f$base, 30); short[1:100] <- "-"
  charmat <- rbind(P1 = f$base, P2 = f$mutate(f$base, 10), P3 = f$mutate(f$base, 20), WRONG = short)
  out <- malaviR:::.malavi_genus_outliers(charmat, c(rep("Plasmodium", 3), "Haemoproteus"))
  expect_equal(out$lineage, "WRONG")
  expect_equal(out$nearest_comparable, 379L)
  ## an unknown genus neither reports nor counts as a neighbour
  charmat <- rbind(P1 = f$base, P2 = f$mutate(f$base, 10), P3 = f$mutate(f$base, 20), WRONG = f$mutate(f$base, 30))
  out <- malaviR:::.malavi_genus_outliers(charmat, c(rep("Plasmodium", 3), NA))
  expect_equal(nrow(out), 0L)
  out <- malaviR:::.malavi_genus_outliers(charmat, c("Plasmodium", "Plasmodium", NA, "Haemoproteus"))
  expect_equal(nrow(out), 0L)
})

test_that("on the release, every reported genus outlier really sits among another genus", {
  ctx <- malaviR:::.malavi_issue_context("latest")
  out <- malaviR:::.malavi_genus_outliers_ctx(ctx)
  skip_if(nrow(out) == 0, "no lineage contradicts its nearest sequences in this release")
  expect_true(all(out$genus != out$neighbor_genus))
  expect_true(all(out$n_neighbors >= malaviR:::.MALAVI_GENUS_MIN_NEIGHBORS))
  expect_true(all(out$n_neighbors_other >= malaviR:::.MALAVI_GENUS_MIN_SHARE * out$n_neighbors))
  expect_true(all(out$nearest_mismatches <= malaviR:::.MALAVI_GENUS_MAX_MISMATCH))
  expect_true(all(out$nearest_comparable >= malaviR:::.MALAVI_GENUS_MIN_COMPARABLE))
  ## the sentence states the count, the rule once, and every lineage by name;
  ## the per-lineage detail stays in the table, not the sentence
  registry <- malaviR:::.malavi_issue_registry()
  issue <- registry[[which(vapply(registry, function(i) i$title, character(1)) ==
                             "Parasite genus contradicts the nearest sequences")]]
  txt <- issue$describe(issue$check(ctx), ctx)
  expect_true(startsWith(txt, paste0(nrow(out), " lineage")))
  expect_true(grepl("within five mismatches", txt, fixed = TRUE))
  expect_true(grepl(paste0("Lineages are: ", paste(sort(out$lineage), collapse = ", "), "."),
                    txt, fixed = TRUE))
  expect_false(grepl("the nearest is", txt, fixed = TRUE))
  ## the second call reads the cache rather than recomputing
  expect_identical(malaviR:::.malavi_genus_outliers_ctx(ctx), out)
})

test_that("the genus-outlier table stored in the bundle is what a fresh computation gives", {
  ## the release build stores the table so malavi_issues() stays quick; this
  ## recomputes it from the bundled alignment (about half a minute) so a stale
  ## table cannot survive a rebuild of the alignment without one
  skip_if(identical(Sys.getenv("MALAVI_SKIP_SLOW"), "true"), "slow test skipped")
  bundle <- malaviR:::.malavi_load("latest")
  skip_if(is.null(bundle$genus_outliers), "this bundle carries no stored genus-outlier table")
  fresh <- malaviR:::.malavi_genus_outliers_bundle(bundle)
  expect_equal(bundle$genus_outliers, fresh)
})
