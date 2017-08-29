
<!-- README.md is generated from README.Rmd. Please edit that file -->
malaviR
=======

This package is an `R` interface to the global avian haemosporidian database MalAvi (<http://mbio-serv2.mbioekol.lu.se/Malavi/>). The package includes functions for downloading data from the MalAvi website directly into your `R` environment and some basic utility functions for manipulating those data. Furthermore, you can use `malaviR` to BLAST your haemosporidian DNA sequences against the MalAvi database programmatically, which should facilitate comparisons with MalAvi.

The package also includes a key for linking the host taxonomic classifications in MalAvi with the avian taxonomic classifications found on <http://birdtree.org/>.

Installation
------------

You can install `malaviR` from github with:

``` r
# install.packages("devtools")
devtools::install_github("vincenzoaellis/malaviR")
```

Then you can load it in your `R` session with:

``` r
library(malaviR)
```

Download tables from MalAvi
---------------------------

There are nine tables that can be downloaded from the MalAvi website that summarize the database. Here's how you would download the table linking morphological species to genetic lineages:

``` r
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

You can BLAST your own sequences against the MalAvi database using the `blast_malavi()` function, which leverages the MalAvi website's already existing BLAST capabilities. Your input sequence just needs to be specified as a character string.

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

The `taxonomy` help file notes when the comparison was made (I will update it occasionally) and lists the multiple MalAvi host species that correspond to single host species in the phylogenetic analysis (there are several).

Identifying sister lineages from a given node
---------------------------------------------

Phylogenetic analysis of the MalAvi lineages is difficult at the moment due to the limited sequence data available for most lineages. Often, the best one can do is to study sister lineages in the phylogeny that are linked by well supported nodes. The function `sister_taxa()` gives a list of lineages on either side of a specific node in a phylogeny.

``` r
## simulate a phylogenetic tree with 10 taxa using the rtree() function in the ape package
library(ape)
tree <- rtree(n=10)

## the node labels of the tree can then be examined (not run)
plot(tree)
nodelabels()
```

![](README-example%205-1.png)

``` r

## the root of the tree has node label "11" and we can extract sister lineages from the root;
## the sister taxa are grouped into two clades with arbitrary labels of "1" and "2"
sis.tax.df <- sister_taxa(tree, 11)
sis.tax.df # check it out
#>    ancestral.node sister.clade taxa
#> 1              11            1   t1
#> 2              11            1   t5
#> 3              11            1   t8
#> 4              11            1   t4
#> 5              11            1   t3
#> 6              11            1  t10
#> 7              11            1   t7
#> 8              11            1   t2
#> 9              11            2   t6
#> 10             11            2   t9
```

This function can be used to identify lineages for further analysis or for visualization purposes.
