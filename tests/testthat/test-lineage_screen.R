## Synthetic 12 bp, frame-1 alignment with controlled singleton substitutions.
## Consensus is ATG TTT GGG CCC (M F G P). One sequence carries a singleton
## non-synonymous change (pos 4 t->a: codon 2 TTT->ATT, F->I) and another a
## singleton synonymous change (pos 6 t->c: codon 2 TTT->TTC, still F). Names use
## the MalAvi genus-prefix scheme so lineage/genus parsing works.
make_ref_screen <- function() {
  cons <- "atgtttgggccc"
  non  <- "atgattgggccc"   # singleton nonsynonymous (pos 4)
  syn  <- "atgttcgggccc"   # singleton synonymous (pos 6)
  c(P_REF01 = cons, P_REF02 = cons, P_REF03 = cons, P_REF04 = cons,
    P_REF05 = cons, P_NON01 = non, P_SYN01 = syn)
}

test_that("lineage_screen counts and classifies singleton substitutions", {
  res <- lineage_screen(reference = make_ref_screen(), studies = FALSE)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 7)
  expect_true(all(res$parasite_genus == "Plasmodium"))

  non <- res[res$lineage == "NON01", ]
  expect_equal(non$n_singleton_substitutions, 1)
  expect_equal(non$n_singleton_nonsynonymous, 1)
  expect_equal(non$n_singleton_synonymous, 0)

  syn <- res[res$lineage == "SYN01", ]
  expect_equal(syn$n_singleton_substitutions, 1)
  expect_equal(syn$n_singleton_nonsynonymous, 0)
  expect_equal(syn$n_singleton_synonymous, 1)

  ## the consensus sequences carry no singleton substitutions
  cons <- res[grepl("^REF", res$lineage), ]
  expect_true(all(cons$n_singleton_substitutions == 0))
})

test_that("lineage_studies counts distinct studies per lineage", {
  studies <- lineage_studies()
  expect_s3_class(studies, "data.frame")
  expect_true(all(c("lineage", "parasite_genus", "n_studies",
                    "n_host_records", "n_countries") %in% names(studies)))
  expect_false("references" %in% names(studies))

  ## SGS1 (P. relictum) is one of the most widely reported lineages
  sgs1 <- studies[studies$lineage == "SGS1", ]
  expect_equal(nrow(sgs1), 1)
  expect_gt(sgs1$n_studies, 50)
  expect_lte(sgs1$n_studies, sgs1$n_host_records)
  expect_gt(sgs1$n_countries, 1)              # found across many countries
  expect_lte(sgs1$n_countries, sgs1$n_host_records)
  expect_equal(sgs1$parasite_genus, "Plasmodium")

  ## the references list is opt-in
  with_refs <- lineage_studies(references = TRUE)
  expect_true("references" %in% names(with_refs))
})

test_that("lineage_screen can restrict the screen to one genus", {
  full <- lineage_screen(studies = FALSE)
  plas <- lineage_screen(genus = "Plasmodium", studies = FALSE)

  ## only Plasmodium lineages remain, and there are fewer than the full pool
  expect_true(all(plas$parasite_genus == "Plasmodium"))
  expect_lt(nrow(plas), nrow(full))
  expect_equal(nrow(plas), sum(full$parasite_genus == "Plasmodium", na.rm = TRUE))

  ## an unknown genus value is rejected by match.arg
  expect_error(lineage_screen(genus = "Babesia", studies = FALSE))
})

test_that("lineage_screen joins study counts to the bundled alignment", {
  res <- lineage_screen()
  expect_true(all(c("n_studies", "in_hosts_table", "n_singleton_nonsynonymous")
                  %in% names(res)))
  ## SGS1 is present, has a sequence, and is reported by many studies
  sgs1 <- res[res$lineage == "SGS1", ]
  expect_equal(nrow(sgs1), 1)
  expect_true(sgs1$in_hosts_table)
  expect_gt(sgs1$n_studies, 50)
})
