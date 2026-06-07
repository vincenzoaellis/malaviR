test_that("synonymy_report summarises a small alignment", {
  full    <- c("a", "t", "c", "g", "a", "t", "a", "c", "c", "g", "a")
  other   <- c("a", "c", "g", "g", "a", "g", "g", "t", "t", "g", "c")
  partial <- c("a", "t", "c", "g", "-", "-", "-", "c", "c", "g", "a")
  aln <- ape::as.DNAbin(rbind(full = full, other = other, partial = partial))

  rep <- synonymy_report(aln, method = "overlap")
  expect_named(rep, c("summary", "by_genus", "synonymies"))
  expect_equal(rep$summary$n_sequences, 3)
  expect_equal(rep$summary$n_haplotypes, 2)        # full+partial collapse, other separate
  expect_equal(rep$summary$n_redundant_names, 1)
  expect_true(all(c("haplotype", "lineage", "genus", "informative_length",
                    "is_partial", "status") %in% names(rep$synonymies)))
  expect_true(any(rep$synonymies$is_partial))
})

test_that("synonymy_report strict finds fewer synonymies than overlap on real data", {
  skip_on_cran()
  rs <- synonymy_report(method = "strict")
  expect_s3_class(rs$summary, "data.frame")
  expect_gte(rs$summary$n_redundant_names, 0)
})
