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
  ## Icterus chrysocephalus used to be the example here. Since 1.1.2 splits
  ## clootl's semicolon-joined synonym cells it resolves through the synonym step
  ## instead (BirdLife lists it under Icterus cayanensis), which is the better
  ## answer -- so a name the legacy key alone still reaches is used instead.
  expect_message(res <- match_taxonomy("Luscinia cyanura"), "legacy")
  row <- res$key[res$key$malavi_species == "Luscinia cyanura", ]
  expect_equal(row$ebird_species, "Tarsiger cyanurus")
  expect_equal(row$match_type, "legacy")
})

test_that("a name inside a joined synonym cell no longer needs the legacy bridge", {
  ## The other half of the change above: Icterus chrysocephalus sits inside
  ## clootl's BirdLife cell "Icterus chrysocephalus;Icterus cayanensis", which the
  ## old exact-string compare could not see. It reached the same species only by
  ## falling through to the legacy key; now the synonym step resolves it directly
  ## and says that eBird lumps it.
  res <- match_taxonomy("Icterus chrysocephalus")
  row <- res$key[res$key$malavi_species == "Icterus chrysocephalus", ]
  expect_equal(row$ebird_species, "Icterus cayanensis")
  expect_equal(row$match_type, "synonym:BirdLife-lump")
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

test_that(".audit_taxonomy surfaces only low-evidence rows", {
  taxonomy <- NULL
  data("taxonomy", package = "malaviR", envir = environment())
  aud <- malaviR:::.audit_taxonomy(taxonomy)
  expect_true(all(c("malavi_species", "ebird_species", "match_type", "reason")
                  %in% names(aud)))
  expect_true(all(aud$reason %in% c("legacy", "retired_genus", "weak_reassignment")))
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

test_that("match_taxonomy warns when it is asked to resolve names without family/order", {
  ## The family/order-constrained epithet match cannot run without them and used to skip
  ## in silence, so 173 of MalAvi's own host names looked unmatched in a hand-rolled call
  ## while the same function resolves every one of them when called with no arguments.
  expect_warning(res <- match_taxonomy("Grus leucogeranus"), "family/order")
  expect_equal(res$key$match_type, "none")
})

test_that("supplying family and order resolves the reassignment, and does not warn", {
  expect_silent(res <- match_taxonomy("Grus leucogeranus",
                                      family = "Gruidae", order = "Gruiformes"))
  expect_equal(res$key$match_type, "reassigned:family")
  expect_equal(res$key$ebird_species, "Leucogeranus leucogeranus")
})

test_that("a name resolved by the legacy bridge does not trigger the warning", {
  ## REGRESSION: the warning first fired at step 5, before the legacy bridge at step 6 had
  ## had its turn, so it complained about names that went on to resolve perfectly well.
  expect_message(res <- match_taxonomy("Luscinia cyanura"), "legacy")
  expect_equal(res$key$match_type, "legacy")
})

test_that("the false species matches corrected in 1.1.2 stay corrected", {
  ## Eight host names in the shipped crosswalk resolved to the wrong bird, all by
  ## the same route: .epithet_reassign() builds its candidate pool from MalAvi's
  ## FAMILY_NAME, and when MalAvi files a genus under a family clootl now fills
  ## with different birds, a lone same-epithet member of that pool wins. Tiaris
  ## obscura, a Peruvian grassquit, came back as Akialoa obscura -- an extinct
  ## Hawaiian honeycreeper -- so three infection records landed on a tip no
  ## grassquit ever occupied. Two further rows were concept errors from the
  ## synonym and legacy steps. Locked here because these are data, not logic:
  ## a rebuild that loses a manual override would silently reintroduce them.
  expected <- c(
    "Alethe poliocephala"    = "Chamaetylas poliocephala",
    "Loxigilla violacea"     = "Melopyrrha violacea",
    "Thraupis cyanocephala"  = "Sporathraupis cyanocephala",
    "Tiaris obscura"         = "Asemospiza obscura",
    "Turdoides gularis"      = "Argya gularis",
    "Hemispingus frontalis"  = "Sphenopsis frontalis",
    "Phaethornis baroni"     = "Phaethornis longirostris",
    "Vermivora pinus"        = "Vermivora cyanoptera",
    "Kittacincla malabarica" = "Copsychus malabaricus",
    "Zosterops capensis"     = "Zosterops virens",
    ## correct cross-family transfers that the genus guard declines, so they are
    ## carried by manual rows and would otherwise regress to "none"
    "Pitohui ferrugineus"    = "Pseudorectes ferrugineus",
    "Pitohui incertus"       = "Pseudorectes incertus")

  got <- taxonomy$ebird_species[match(names(expected), taxonomy$malavi_species)]
  expect_equal(got, unname(expected))
})

test_that("the genus guard declines an epithet match landing outside the genus's family", {
  ## The mechanism itself, on the real clootl snapshot. clootl files Tiaris in
  ## Thraupidae; Akialoa obscura is Fringillidae, so the family-pool epithet match
  ## must be refused even though MalAvi labels Tiaris as Fringillidae.
  ref <- malaviR:::clootl_ref
  ref$latin_family <- sub(" .*$", "", ref$FAMILY)

  refused <- malaviR:::.epithet_reassign("Tiaris obscura", "Fringillidae",
                                         "Passeriformes", ref)
  expect_true(is.na(refused$ebird))

  ## and it does not get in the way of an ordinary genus transfer within the
  ## family clootl already files the genus in
  ok <- malaviR:::.epithet_reassign("Dendroica virens", "Parulidae",
                                    "Passeriformes", ref)
  expect_equal(ok$ebird, "Setophaga virens")
})

test_that("synonym lookup sees names inside semicolon-joined clootl cells", {
  ## clootl records every name of an authority that splits an eBird species in one
  ## cell, joined by ";". A plain == could never match inside those, so 153 MalAvi
  ## host names were invisible to the synonym step; Phaethornis baroni fell through
  ## to the family pool and came back as Metallura baroni.
  ref <- malaviR:::clootl_ref
  ref$latin_family <- sub(" .*$", "", ref$FAMILY)
  syn <- malaviR:::.syn_table(ref)

  ## the exploded table splits joined cells and keeps the unjoined ones intact
  expect_false(any(grepl(";", syn$synonym, fixed = TRUE)))
  expect_true(any(syn$lump))

  res <- malaviR:::.syn_resolve("Phaethornis baroni", ref, syn)
  expect_equal(res$ebird, "Phaethornis longirostris")
  ## flagged as a lump: BirdLife treats baroni as a species, eBird does not
  expect_match(res$type, "-lump$")

  ## a name in an unjoined cell keeps the plain label
  plain <- syn[!syn$lump, ][1, ]
  expect_false(grepl("-lump", malaviR:::.syn_resolve(plain$synonym, ref, syn)$type))
})

test_that("the audit flags genus transfers clootl cannot check", {
  ## The genus guard is blind whenever clootl has retired the MalAvi genus: there
  ## is no home family to compare the candidate against. Hemispingus frontalis ->
  ## Crithagra frontalis hid in exactly that gap, and was never audited because
  ## the old criterion exempted exact-epithet matches.
  audit <- malaviR:::.audit_taxonomy(taxonomy)
  expect_true("retired_genus" %in% audit$reason)
  ## Hemispingus is now carried by a manual row, so it should NOT be in the audit
  expect_false("Hemispingus frontalis" %in% audit$malavi_species)
  ## but the standard retired-genus transfers are, so the gap stays visible
  expect_true("Dendroica virens" %in% audit$malavi_species)
})
