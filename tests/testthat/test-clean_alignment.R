context("clean alignments")


## make a test alignment with one repeated haplotype
library(ape)
seq1 <- c("a", "t", "c", "g", "a", "t", "a", "c", "c", "g", "a")
seq2 <- c("a", "c", "g", "g", "a", "g", "g", "t", "t", "g", "c")
seq3 <- c("a", "t", "c", "g", "a", "t", "a", "c", "c", "g", "a")
align <- as.DNAbin(rbind(seq1, seq2, seq3))

test_that("clean_alignment returns a DNAbin alignment", {
  out <- clean_alignment(alignment = align)
  expect_is(out$alignment_clean, "DNAbin")
})

test_that("clean_alignment throws an error if alignment has no repeated haplotypes", {
  out <- clean_alignment(alignment = align)
  expect_error(clean_alignment(out$alignment_clean), "The alignment has no repeated haplotypes")
})

test_that("clean_alignment evaluates haplotype format correctly", {
  out.wide <- clean_alignment(alignment = align, haplotype_format_wide = TRUE)
  out.long <- clean_alignment(alignment = align, haplotype_format_wide = FALSE)

  # convert wide to long format using the base function reshape to check that the
  # tidyverse way in clean_alignment is performing correctly
  manual.long <- reshape(out.wide$repeated_haplotypes, varying = c("V1", "V2"),
                         v.names = "Lineage_Name", timevar = "lin_number",
                         times = c("V1", "V2"),
                         direction = "long", idvar = "haplotype")
  manual.long.n <- data.frame(manual.long[complete.cases(manual.long), ])
  row.names(manual.long.n) <- 1:dim(manual.long.n)[1]
  manual.long.n$Lineage_Name <- as.character(manual.long.n$Lineage_Name)
  expect_equal(out.long$repeated_haplotypes, manual.long.n)
})

