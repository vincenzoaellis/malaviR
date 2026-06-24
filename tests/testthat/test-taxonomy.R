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
  expect_message(res <- match_taxonomy("Icterus chrysocephalus"), "legacy")
  row <- res$key[res$key$malavi_species == "Icterus chrysocephalus", ]
  expect_equal(row$ebird_species, "Icterus cayanensis")
  expect_equal(row$match_type, "legacy")
})

test_that("match_taxonomy recovers a same-genus gender shift before the family step", {
  ## Saxicola maura -> Saxicola maurus is a pure gender change within one genus.
  ## The same-genus step must win over the family/order pool (MalAvi files this
  ## host under Turdidae, which previously pulled it to Turdus torquatus).
  res <- match_taxonomy("Saxicola torquata")
  row <- res$key[1, ]
  expect_equal(row$ebird_species, "Saxicola torquatus")
  expect_equal(row$match_type, "reassigned:genus")
})

test_that("match_taxonomy disambiguates a synonym shared by two eBird species", {
  ## Howard & Moore "Trochalopteron cachinnans" is carried by both Montecincla
  ## cachinnans and M. jerdoni; the epithet must decide, not row order.
  res <- match_taxonomy("Trochalopteron cachinnans")
  row <- res$key[1, ]
  expect_equal(row$ebird_species, "Montecincla cachinnans")
  expect_equal(row$match_type, "synonym:HowardMoore")
})

test_that("match_taxonomy does not make a cross-family epithet match on a mislabeled family", {
  ## MalAvi files Oriolus brachyrhynchus under Corvidae, which once forced a
  ## false stem match to Corvus brachyrhynchos. A manual override now fixes it,
  ## and in no case should an oriole resolve to a crow.
  res <- match_taxonomy("Oriolus brachyrhynchus")
  row <- res$key[1, ]
  expect_equal(row$ebird_species, "Oriolus brachyrynchus")
  expect_false(identical(row$ebird_species, "Corvus brachyrhynchos"))
})

test_that("every reassigned row preserves the host's epithet (gender-relaxed)", {
  ## reassigned:* matches are defined by epithet agreement, so the bundled
  ## crosswalk should never carry a reassigned row whose epithet stem differs --
  ## a guard against a future resolver change reintroducing epithet collisions.
  taxonomy <- NULL
  data("taxonomy", package = "malaviR", envir = environment())
  re <- taxonomy[grepl("^reassigned", taxonomy$match_type), ]
  mal <- malaviR:::.epithet_stem(malaviR:::.epithet(re$malavi_species))
  eb  <- malaviR:::.epithet_stem(malaviR:::.epithet(re$ebird_species))
  expect_equal(mal, eb)
})

test_that(".audit_taxonomy surfaces only weak reassignments and legacy rows", {
  taxonomy <- NULL
  data("taxonomy", package = "malaviR", envir = environment())
  aud <- malaviR:::.audit_taxonomy(taxonomy)
  expect_true(all(c("malavi_species", "ebird_species", "match_type", "reason")
                  %in% names(aud)))
  expect_true(all(aud$reason %in% c("legacy", "weak_reassignment")))
  ## all legacy rows are reported
  expect_equal(sum(aud$reason == "legacy"), sum(taxonomy$match_type == "legacy"))
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
