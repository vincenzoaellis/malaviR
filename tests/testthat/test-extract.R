test_that("malavi_versions and malavi_version report the bundled release", {
  v <- malavi_versions()
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

  ## "all" returns the five tables
  expect_length(extract_table("all"), 5)

  ## bad name is an error, not a message (issue #1)
  expect_error(extract_table("not a real table"), "choose one")
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
