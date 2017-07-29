
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

There are nine tables that can be downloaded from the MalAvi website summarizing the dataset. Here's how you would download the table linking morphological species to genetic lineages from the MalAvi database:

``` r
library(malaviR)
morph <- extract_table("Morpho Species Summary")
head(morph) # check it out
#>   #no Lineage_Name         genus                   species
#> 1   1     ACCFRA01 Leucocytozoon       Leucocytozoon toddi
#> 2   2      ACCOP01 Leucocytozoon     Leucocytozoon mathisi
#> 3   3       ACNI04 Leucocytozoon     Leucocytozoon mathisi
#> 4   4      ALARV01  Haemoproteus Haemoproteus tartakovskyi
#> 5   5      ALARV02  Haemoproteus Haemoproteus tartakovskyi
#> 6   6      ALARV03  Haemoproteus Haemoproteus tartakovskyi
#>           Reference_Name comment
#> 1 Barraclough et al 2008        
#> 2  Valkiunas et al 2010a        
#> 3  Valkiunas et al 2010a        
#> 4 Zehtindjiev et al 2012        
#> 5 Zehtindjiev et al 2012        
#> 6 Zehtindjiev et al 2012
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

Download the MalAvi sequence alignments
---------------------------------------

Using the `extract_alignment()` function, you can download all of the sequences on MalAvi, or the complete or nearly complete sequences only, or the sequences associated with morphological species only. These alignments will appear as objects of the class `DNAbin`, defined in the `ape` package.

``` r
## download all sequences
all.seqs <- extract_alignment("all seqs")
all.seqs #check it out
#> 2771 DNA sequences in binary format stored in a matrix.
#> 
#> All sequences of same length: 479 
#> 
#> Labels:
#> H_ACAED01
#> P_ACAGR1
#> H_ACAGR2
#> H_ACATEN01
#> L_ACCBRE01
#> L_ACCBRE02
#> ...
#> 
#> Base composition:
#>     a     c     g     t 
#> 0.294 0.133 0.135 0.438
```

The `extract_alignment()` help file lists the names of the sequence alignments that the function understands.

Taxonomic key for host species
------------------------------

The package includes a taxonomic key linking (nearly all) host species names from MalAvi with the species names used in the phylogenetic analysis found on <http://birdtree.org/>. The latter species names are in a column labeled `Jetz.species` and have an underscore between the genus and species as they would appear if you were to download trees from the website. The key is stored in a data object called `taxonomy` which can be called directly once the package is loaded.

``` r
## check out taxonomic host key
head(taxonomy)
#>              species       Jetz.species match
#> 1 Nectarinia sperata Nectarinia_sperata   yes
#> 2    Prinia inornata    Prinia_inornata   yes
#> 3 Acrocephalus aedon Acrocephalus_aedon   yes
#> 4   Emberiza elegans   Emberiza_elegans   yes
#> 5   Emberiza pusilla   Emberiza_pusilla   yes
#> 6        Parus major        Parus_major   yes
```
