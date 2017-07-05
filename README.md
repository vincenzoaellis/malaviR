
<!-- README.md is generated from README.Rmd. Please edit that file -->
malaviR
=======

The goal of this package is to be an R interface to the global avian haemosporidian database MalAvi (<http://mbio-serv2.mbioekol.lu.se/Malavi/>). The package includes functions for downloading data from the MalAvi website directly into your R environment. You can also use malaviR to BLAST your haemosporidian DNA sequences against the MalAvi database programmatically, which should facilitate comparisons with MalAvi.

The package also includes a key for linking the host taxonomic classifications in MalAvi with the avian taxonomic classifications found on <http://birdtree.org/>.

Installation
------------

You can install malaviR from github with:

``` r
# install.packages("devtools")
devtools::install_github("vincenzoaellis/malaviR")
```

Download tables from MalAvi
---------------------------

There are nine tables that can be downloaded from the MalAvi website summarizing the dataset. Here's how you would download the table of references from the MalAvi database:

``` r
library(malaviR)
refs <- extract_table("Table of References")
head(refs) # check it out
#>   #no             Reference_Name year
#> 1   1           Alley et al 2010 2010
#> 2   2         Argilla et al 2013 2013
#> 3   3           Aysul et al 2013 2013
#> 4   4     Baillie & Brunton 2011 2011
#> 5   5 Balasubramaniam et al 2013 2013
#> 6   6           Banda et al 2012 2012
#>                                                                                                                                            title
#> 1          Concurrent avian malaria and avipox virus infection in translocated South Island Saddlebacks (Philesturnus carunculatus carunculatus)
#> 2 High prevalence of Leucocytozoon spp. in the endangered yellow-eyed penguin (Megadyptes antipodes) in the sub-Antarctic regions of New Zealand
#> 3                                      Detection and molecular characterization of a Haemoproteus lineage in a tawny owl (Strix Aluco) in Turkey
#> 4                     Diversity, distribution and biogeographical origins of Plasmodium parasites from New Zealand bellbird (Anthornis melanura)
#> 5                                                           Prevalence and diversity of avian haematozoa in three speciesof Australian passerine
#> 6                                                                                A cluster of avian malaria cases in a kiwi management programme
#>                                           journal               volume
#> 1                  New Zealand Veterinary Journal           58:218-223
#> 2                                    Parasitology          140:672-682
#> 3 Ankara Universitesi Veteriner Fakultesi Dergisi           60:179-183
#> 4                                    Parasitology        138:1843-1851
#> 5                                             Emu          113:353-358
#> 6                  New Zealand Veterinary Journal DOI:10.1080/00480169
#>       study_type
#> 1 Single species
#> 2 Single species
#> 3 Single species
#> 4 Single species
#> 5    Two Species
#> 6     Veterinary
```

The `extract_table()` help file lists all nine tables you can download directly into `R`. Or you can specify `"all"` and get them all as a list.

BLAST a sequence to MalAvi
--------------------------

You can BLAST your own sequences against the MalAvi database using the `blast_malavi()` function. Your input sequence just needs to be specified as a character string.

``` r
## define a sequence. This is the Plasmodium parasite ACAGR1
ACAGR1 <- "GCAACTGGTGCTTCATTTGTATTTATTTTAACTTATTTACATATTTTAAGAGGATTAAATTATTCATATTCATATTTACCTTTATCATGGATATCTGGATTAATAATATTTTTAATATCTATAGTAACAGCTTTTATGGGTTACGTATTACCTTGGGGTCAAATGAGTTTCTGGGGTGCTACCGTAATAACTAATTTATTATATTTTATACCTGGACTAGTTTCATGGATATGTGGTGGATATCTTGTAAGTGACCCAACCTTAAAAAGATTCTTTGTACTACATTTTACATTTCCTTTTATAGCTTTATGTATTGTATTTATACATATATTCTTTCTACATTTACAAGGTAGCACAAATCCTTTAGGGTATGATACAGCTTTAAAAATACCCTTCTATCCAAATCTTTTAAGTCTTGATATTAAAGGATTTAATAATGTATTAGTATTATTTTTAGCACAAAGTTTATTTGGAATACT"

## BLAST it against the MalAvi database and save the top five hits to a data frame
hits <- blast_malavi(ACAGR1)
#> Submitting with 'sequence'
hits # check it out
#>    Lineage Score Identities  Gaps    Strand Coverage Perfect.Match
#> 1   ACAGR1   885    479/479 0/479 Plus/Plus      479           Yes
#> 2     SGS1   880    478/479 0/479 Plus/Plus      479            No
#> 3  CXPIP20   876    476/477 0/477 Plus/Plus      477            No
#> 4     YWT4   874    477/479 0/479 Plus/Plus      479            No
#> 5 SERCAN01   874    477/479 0/479 Plus/Plus      479            No
```
