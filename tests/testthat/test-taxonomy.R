test_that("match_taxonomy matches exact names and flags unmatchable ones", {
  res <- match_taxonomy(c("Turdus merula", "Cyanistes caeruleus",
                          "Anas sp.", "Luscinia luscinia x megarhynchos"))
  expect_named(res, c("key", "differences"))
  expect_s3_class(res$key, "data.frame")

  mt <- res$key$match_type[match(c("Turdus merula", "Anas sp.",
                                   "Luscinia luscinia x megarhynchos"),
                                 res$key$malavi_species)]
  expect_equal(mt, c("exact", "generic", "generic"))

  ## an exact match carries a phylogeny tip label and family
  tm <- res$key[res$key$malavi_species == "Turdus merula", ]
  expect_equal(tm$ebird_species, "Turdus merula")
  expect_false(is.na(tm$ott_name))

  ## differences excludes the exact matches
  expect_false(any(res$differences$match_type == "exact"))
})

test_that("match_taxonomy recovers a genus reassignment within the host family", {
  res <- match_taxonomy("Anas clypeata", family = "Anatidae",
                        order = "Anseriformes")
  row <- res$key[1, ]
  expect_equal(row$ebird_species, "Spatula clypeata")
  expect_equal(row$match_type, "reassigned:family")
})

test_that("match_taxonomy leaves an epithet ambiguous within family unmatched", {
  ## 'americana' is shared by Mareca americana and Aythya americana in Anatidae;
  ## Querquedula is an old duck genus not in the override/legacy keys
  res <- match_taxonomy("Querquedula americana", family = "Anatidae",
                        order = "Anseriformes")
  expect_equal(res$key$match_type[1], "none")
})

test_that("match_taxonomy applies a maintainer manual override", {
  res <- match_taxonomy("Anas americana")
  row <- res$key[1, ]
  expect_equal(row$ebird_species, "Mareca americana")
  expect_equal(row$match_type, "manual")
})

test_that("match_taxonomy flags 'spp' as generic", {
  res <- match_taxonomy("Somateria spp")
  expect_equal(res$key$match_type[1], "generic")
})

test_that("match_taxonomy bridges a name via the legacy hand-curated key", {
  expect_message(res <- match_taxonomy("Buarremon basilicus"), "legacy")
  row <- res$key[res$key$malavi_species == "Buarremon basilicus", ]
  expect_equal(row$ebird_species, "Arremon torquatus")
  expect_equal(row$match_type, "legacy")
})

test_that("clootl_taxonomy_version returns a year", {
  yr <- clootl_taxonomy_version()
  expect_true(is.numeric(yr))
  expect_gte(yr, 2021)
})

test_that("bundled taxonomy dataset has the expected shape", {
  taxonomy <- NULL
  data("taxonomy", package = "malaviR", envir = environment())
  expect_s3_class(taxonomy, "data.frame")
  expect_true(all(c("malavi_species", "ebird_species", "ott_name",
                    "order", "family", "match_type") %in% names(taxonomy)))
})
