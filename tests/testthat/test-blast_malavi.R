context("BLAST MalAvi")

## define a sequence to blast to MalAvi
ACAGR1 <- "GCAACTGGTGCTTCATTTGTATTTATT
TTAACTTATTTACATATTTTAAGAGGATTAAATTATTC
ATATTCATATTTACCTTTATCATGGATATCTGGATTAA
TAATATTTTTAATATCTATAGTAACAGCTTTTATGGGT
TACGTATTACCTTGGGGTCAAATGAGTTTCTGGGGTGC
TACCGTAATAACTAATTTATTATATTTTATACCTGGAC
TAGTTTCATGGATATGTGGTGGATATCTTGTAAGTGAC
CCAACCTTAAAAAGATTCTTTGTACTACATTTTACATT
TCCTTTTATAGCTTTATGTATTGTATTTATACATATAT
TCTTTCTACATTTACAAGGTAGCACAAATCCTTTAGGG
TATGATACAGCTTTAAAAATACCCTTCTATCCAAATCT
TTTAAGTCTTGATATTAAAGGATTTAATAATGTATTAG
TATTATTTTTAGCACAAAGTTTATTTGGAATACT"

test_that("blast_malavi returns a data frame when given real sequence", {

  skip_on_cran() # bc downloading data can take a while

  df <- blast_malavi(ACAGR1)
  expect_is(df, "data.frame")
})

test_that("blast_malavi returns the correct number of hits", {

  skip_on_cran() # bc downloading data can take a while

  df <- blast_malavi(ACAGR1, hits = 2)
  expect_that(dim(df)[1], equals(2))
})

test_that("blast_malavi produces warning and data frame when a match is not found", {

  skip_on_cran() # bc downloading data can take a while

  not_real <- "AAAAAAAAAAAAATTTTTTTTTTTTTGGGGGGGGGGGGGGCCCCCCCCCCCC"
  expect_warning(df <- blast_malavi(not_real), "No hits found: check your input sequence")
  expect_is(df, "data.frame")
})

