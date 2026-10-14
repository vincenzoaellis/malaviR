test_that("malavi_version reports the bundled release", {
  v <- malavi_version("all")
  expect_type(v, "character")
  expect_gte(length(v), 1)
  expect_identical(malavi_version(), v[1])
  expect_match(v[1], "^[0-9]{4}-[0-9]{2}-[0-9]{2}$")
})

test_that("extract_table returns tables and errors on a bad name", {
  h <- extract_table("Hosts and Sites Table")
  expect_s3_class(h, "data.frame")
  expect_true("LINEAGE_NAME" %in% names(h))

  ## descriptive name and snake_case key are equivalent
  expect_identical(extract_table("references"), extract_table("Table of References"))

  ## "all" returns the seven tables (five, plus the Lineage History since 2026-09-25
  ## and the Genomic Resources since malaviR 1.3.4)
  expect_length(extract_table("all"), 7)

  ## bad name is an error, not a message (issue #1)
  expect_error(extract_table("not a real table"), "choose one")

  ## a non-scalar 'table' is rejected with a clear message
  expect_error(extract_table(c("references", "vector_data")), "single table name")
})

test_that("bundled tables carry no blank-cell export artifact", {
  ## A malformed spreadsheet export can make blank cells read back as the
  ## literal first-column name (e.g. "LINEAGE_NAME") across many columns. Guard
  ## against ever shipping that: no value should equal its own table's first
  ## column name. (See data-raw/process_release.R, which strips it.)
  tables <- extract_table("all")
  for (nm in names(tables)) {
    df <- tables[[nm]]
    col1 <- names(df)[1]
    bad <- sum(vapply(df, function(col) sum(col == col1, na.rm = TRUE), integer(1)))
    expect_identical(bad, 0L, info = paste("artifact cells in", nm))
  }
})

test_that("extract_alignment returns a DNAbin and subsets by genus", {
  a <- extract_alignment()
  expect_s3_class(a, "DNAbin")

  p <- extract_alignment(genus = "Plasmodium")
  expect_true(all(grepl("^P_", rownames(p))))
  expect_lt(nrow(p), nrow(a))

  ph <- extract_alignment(genus = c("Plasmodium", "Haemoproteus"))
  expect_true(all(grepl("^(P|H)_", rownames(ph))))
  expect_gt(nrow(ph), nrow(p))
})

test_that('extract_table("all") returns the same tables as the single-table calls', {
  ## Until 1.1.2 the "all" branch handed back the raw bundled tables while every
  ## single-table call tidied whitespace, so the same table differed in 3,143
  ## cells depending on how it was asked for -- and a host-name join built from
  ## the list silently missed the rows with a stray newline.
  all5 <- extract_table("all")
  for (nm in names(all5)) {
    expect_identical(all5[[nm]], extract_table(nm), info = nm)
  }
})

test_that("the whitespace tidy also removes non-breaking spaces", {
  ## [[:space:]] does not match U+00A0, so 22 GENBANK_ACC values in the
  ## 2026-03-23 release kept an invisible leading non-breaking space that
  ## trimws() could not take off and that breaks an accession join.
  nbsp <- " "
  g <- extract_table("Grand Lineage Summary")
  expect_false(any(grepl(nbsp, g$GENBANK_ACC, fixed = TRUE), na.rm = TRUE))

  ## and directly, so the test does not depend on this release still having them
  df <- data.frame(x = c(paste0(nbsp, "AB123"), "  spaced\tout\n "),
                   stringsAsFactors = FALSE)
  expect_equal(malaviR:::.clean_table_ws(df)$x, c("AB123", "spaced out"))
})


test_that("the Lineage History table is always there, with its ten columns", {
  ## Releases bundled before 2026-09-25 carry no such table; it is returned empty
  ## rather than failing, so code written against it runs on every release.
  h <- extract_table("Lineage History")
  expect_s3_class(h, "data.frame")
  expect_identical(names(h), c("LINEAGE_NAME", "EVENT", "OLD_NAME", "NEW_NAME",
                               "OLD_SEQUENCE", "NEW_SEQUENCE", "DATE", "REASON",
                               "SOURCE", "RELEASE"))
  expect_identical(extract_table("lineage_history"), h)
  expect_true(all(h$EVENT %in% c("sequence_updated", "renamed", "merged_into", "removed")))
  expect_true("lineage_history" %in% names(extract_table("all")))
})

test_that("the Genomic Resources table is always there, with its eleven columns", {
  ## Release format since 2026-09-25, first filled 2026-10-14; an older bundle
  ## returns it empty (attribute malavi_absent) rather than failing.
  g <- extract_table("Genomic Resources")
  expect_s3_class(g, "data.frame")
  expect_identical(names(g), c("LINEAGE_NAME", "RESOURCE_TYPE", "ACCESSION", "DATABASE",
                               "LENGTH_BP", "DESCRIPTION", "HOST_SPECIES",
                               "REFERENCE_NAME", "DOI", "SOURCE", "NOTES"))
  expect_identical(extract_table("genomic_resources"), g)
  expect_true(all(g$RESOURCE_TYPE %in% c("mitogenome", "mitogenome_partial", "nuclear_genome",
                                         "apicoplast_genome", "transcriptome", "other")))
  expect_true("genomic_resources" %in% names(extract_table("all")))
})
