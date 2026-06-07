# malaviR

An R interface to [MalAvi](http://mbio-serv2.mbioekol.lu.se/Malavi/), the public
database of avian haemosporidian (malaria and related) parasite mtDNA cytochrome
*b* lineages.

> **Note on this release.** MalAvi is no longer hosted as a queryable website, so
> the functions that used to download data from the MalAvi server no longer work.
> This version of `malaviR` instead **bundles a dated snapshot of the MalAvi
> database** (tables and the sequence alignment) inside the package, and provides
> a local, BLAST-like search to replace the old online BLAST. I used Claude Code
> (Opus 4.8) to rewrite the package. I believe everything is working as expected,
> but please open an issue if you find mistakes.

## Installation

```r
# install.packages("remotes")
remotes::install_github("vincenzoaellis/malaviR")
```

The core data functions need only CRAN packages. Two optional features need extra
packages:

```r
# for blast_malavi() — local BLAST-like search (needs R >= 4.4)
# install.packages("BiocManager")
BiocManager::install(c("DECIPHER", "Biostrings"))

# clootl is only needed to *rebuild* the bundled taxonomy, not to use it
install.packages("clootl")
```

## What's in the package

The bundled database is identified by its release date (e.g. `2026-03-23`).

```r
library(malaviR)

malavi_version()        # release date of the bundled MalAvi snapshot
malavi_version("all")   # all bundled releases
```

### Data tables and the alignment

```r
hosts <- extract_table("Hosts and Sites Table")   # one of five MalAvi tables
aln   <- extract_alignment()                       # the cyt-b alignment (DNAbin)
plas  <- extract_alignment(genus = "Plasmodium")   # subset by parasite genus
```

### Local BLAST-like search

`blast_malavi()` searches a query sequence against the bundled MalAvi alignment
using a pre-built [DECIPHER](https://decipher.codes/) index, returning the most
similar lineages. (Requires DECIPHER >= 3.0 and Biostrings.)

```r
query <- gsub("-", "", paste(as.character(aln[1, ]), collapse = ""))
blast_malavi(query, top_n = 5)
```

### Repeated haplotypes ("synonymies")

Different MalAvi names are sometimes assigned to the same sequence, which inflates
estimates of parasite diversity — see
[Tamayo-Quintero et al. 2025](https://doi.org/10.1371/journal.ppat.1012911).
`synonymy_report()` quantifies the problem and points to the lineages to check;
`clean_alignment()` produces a de-duplicated alignment, letting you choose which
name to keep. By default it keeps the most complete sequence in each group
(deterministic); use `select = "random"` for a quick random pick, or `keep =` to
override specific groups (see `?clean_alignment`).

```r
synonymy_report()$summary             # how many names share a haplotype
res <- clean_alignment(aln, method = "overlap")   # keeps most complete per group
head(res$synonymies)

set.seed(1)
res_rand <- clean_alignment(aln, select = "random")   # quick random pick
```

### Host taxonomy

`match_taxonomy()` aligns MalAvi host species names to the modern
[clootl](https://github.com/McTavishLab/clootl) (eBird) avian taxonomy and flags
names that don't match (synonyms, hybrids, `sp.`, etc.). The pre-built crosswalk
for all MalAvi hosts is bundled as the `taxonomy` dataset.

```r
match_taxonomy(c("Turdus merula", "Anas sp."))$key
data(taxonomy)
```

### Phylogenies

`sister_taxa()` returns the tips descending from each child clade at a node, and
`clean_names()` strips the genus prefix from alignment tip labels.

## Updating the bundled database

When a new MalAvi release arrives (a `MalAvi_<date>.zip`), maintainers drop it in
`data-raw/` and run `data-raw/process_release.R` to regenerate the bundled data
and BLAST index, then commit and push.

## Citation

If you use `malaviR`, please cite the package and:

* Ellis VA, Bensch S (2018). Host specificity of avian haemosporidian parasites
  is unrelated among sister lineages but shows phylogenetic signal across larger
  clades. *International Journal for Parasitology* **48**: 897–902.
  <https://doi.org/10.1016/j.ijpara.2018.05.005>
* Bensch S, Hellgren O, Pérez-Tris J (2009). MalAvi: a public database of malaria
  parasites and related haemosporidians in avian hosts based on mitochondrial
  cytochrome *b* lineages. *Molecular Ecology Resources* **9**: 1353–1358.

```r
citation("malaviR")
```

## Issues

Please report problems at
<https://github.com/vincenzoaellis/malaviR/issues>.
