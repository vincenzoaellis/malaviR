test_that("lineage_plausibility returns one row per detection with the expected columns", {
  out <- lineage_plausibility("SGS1", host = "Parus major", country = "Sweden")
  expect_s3_class(out, "malavi_plausibility")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1)
  expect_named(out, c("lineage", "call", "flags", "n_studies", "n_host_records",
                      "n_countries", "host", "host_family", "host_order",
                      "host_recorded", "host_genus_recorded",
                      "host_family_recorded", "host_order_recorded",
                      "n_records_host", "country", "country_recorded",
                      "n_records_country", "region", "region_recorded"))
})

test_that("a well-known lineage in a well-known host and country is not flagged", {
  ## SGS1 in Parus major in Sweden is about as well recorded as MalAvi gets
  out <- lineage_plausibility("SGS1", host = "Parus major", country = "Sweden")
  expect_equal(out$call, "previously_recorded")
  expect_equal(out$flags, "")
  expect_true(out$host_recorded)
  expect_true(out$country_recorded)
  expect_true(out$n_records_host > 0)
  expect_true(out$n_studies > 1)
})

test_that("a host from an unrecorded family is called new_host_family", {
  ## TUPHI01 is a Haemoproteus of thrushes; Parus major is a tit
  out <- lineage_plausibility("TUPHI01", host = "Parus major")
  expect_equal(out$call, "new_host_family")
  expect_false(out$host_recorded)
  expect_false(out$host_family_recorded)
  expect_match(out$flags, "new_host_family")
  ## the host itself is well known in MalAvi, just not for this lineage
  expect_equal(out$host_family, "Paridae")
  expect_false(grepl("host_not_in_malavi", out$flags))
})

test_that("an unrecorded region is called new_location", {
  out <- lineage_plausibility("SGS1", host = "Parus major", region = "Antarctica")
  expect_equal(out$call, "new_location")
  expect_false(out$region_recorded)
  expect_match(out$flags, "new_region")

  ## and a recorded one is not
  eur <- lineage_plausibility("SGS1", host = "Parus major", region = "Europe")
  expect_true(eur$region_recorded)
  expect_equal(eur$call, "previously_recorded")
})

test_that("region names are normalized and unknown ones are rejected", {
  variants <- lineage_plausibility(rep("SGS1", 3),
                                   region = c("Europe", "europe", "EUROPE"))
  expect_true(all(variants$region == "EUROPE"))
  expect_true(all(variants$region_recorded))

  ## MalAvi spells it HAWAI; users will write Hawaii
  expect_equal(lineage_plausibility("SGS1", region = "Hawaii")$region, "HAWAI")

  expect_error(lineage_plausibility("SGS1", region = "Narnia"), "Unknown region")
  expect_error(lineage_plausibility("SGS1", region = "Narnia"), "EUROPE")
})

test_that("a lineage with no host records is reported as unknown, without noise", {
  out <- lineage_plausibility("NOTALINEAGE1", host = "Parus major",
                              country = "Sweden")
  expect_equal(out$call, "lineage_not_in_malavi")
  expect_equal(out$flags, "lineage_not_in_malavi")
  ## the novelty flags would be trivially true and are deliberately suppressed
  expect_false(grepl("new_host", out$flags))
  expect_true(is.na(out$n_studies))
})

test_that("a host that has never been screened is flagged as such", {
  out <- lineage_plausibility("SGS1", host = "Notabird species")
  expect_match(out$flags, "host_not_in_malavi")
  expect_true(is.na(out$host_family))
  expect_false(out$host_recorded)
})

test_that("arguments recycle to length 1 or n, and refuse anything else", {
  out <- lineage_plausibility(c("SGS1", "GRW04"), host = "Parus major")
  expect_equal(nrow(out), 2)
  expect_equal(out$host, rep("Parus major", 2))

  out2 <- lineage_plausibility(c("SGS1", "GRW04"),
                               host = c("Parus major", "Acrocephalus arundinaceus"))
  expect_equal(nrow(out2), 2)
  expect_equal(out2$host, c("Parus major", "Acrocephalus arundinaceus"))

  expect_error(
    lineage_plausibility(c("SGS1", "GRW04", "COLL2"),
                         host = c("Parus major", "Acrocephalus arundinaceus")),
    "length 1 or the same length"
  )
})

test_that("optional arguments left out give NA columns, not false novelty", {
  out <- lineage_plausibility("SGS1")
  expect_true(is.na(out$host))
  expect_true(is.na(out$host_recorded))
  expect_true(is.na(out$country_recorded))
  expect_true(is.na(out$region_recorded))
  ## nothing was asked about, so nothing can be new
  expect_equal(out$call, "previously_recorded")
  expect_equal(out$flags, "")
})

test_that("full alignment names are accepted as lineage names", {
  bare <- lineage_plausibility("SGS1", host = "Parus major")
  full <- lineage_plausibility("P_SGS1", host = "Parus major")
  expect_equal(full$lineage, "SGS1")
  expect_equal(full$call, bare$call)

  ## including the morphospecies-suffixed form
  suffixed <- lineage_plausibility("H_COLL2_Haemoproteus_pallidus")
  expect_equal(suffixed$lineage, "COLL2")
})

test_that("lineage_plausibility stamps the version and prints", {
  out <- lineage_plausibility("SGS1", host = "Parus major")
  expect_equal(attr(out, "malavi_meta")$malavi_version, malavi_version())
  expect_output(print(out), "host and biogeographic plausibility")
  expect_output(print(out), "SGS1")
  expect_error(lineage_plausibility(character(0)), "empty")
})
