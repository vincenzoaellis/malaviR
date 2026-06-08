test_that("clean_names strips the genus prefix and any morpho-species suffix", {
  labs <- c("H_COLL2_Haemoproteus_pallidus", "P_GRW04_Plasmodium_relictum",
            "L_CIAE02")
  expect_equal(clean_names(labs), c("COLL2", "GRW04", "CIAE02"))
})

test_that("clean_names keep.genus returns expanded genus and cleaned lineage", {
  out <- clean_names(c("H_COLL2_Haemoproteus_pallidus", "L_CIAE02"),
                     keep.genus = TRUE)
  expect_s3_class(out, "data.frame")
  expect_equal(out$parasiteGenus, c("Haemoproteus", "Leucocytozoon"))
  expect_equal(out$Lineage_Name, c("COLL2", "CIAE02"))
})

test_that("clean_names handles mixed-length labels without recycling errors", {
  ## a 2-part label mixed with 4-part labels used to rely on rbind recycling
  labs <- c("L_CIAE02", "H_COLL2_Haemoproteus_pallidus", "P_SGS1")
  expect_equal(clean_names(labs), c("CIAE02", "COLL2", "SGS1"))
})
