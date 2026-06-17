# malaviR

<!-- badges: start -->
[![r-universe status](https://vincenzoaellis.r-universe.dev/badges/malaviR)](https://vincenzoaellis.r-universe.dev/malaviR)
[![r-universe downloads](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fvincenzoaellis.r-universe.dev%2Fapi%2Fpackages%2FmalaviR&query=%24._downloads.count&label=downloads&color=blue)](https://vincenzoaellis.r-universe.dev/malaviR)
<!-- badges: end -->

An R interface to [MalAvi](https://wimanet-science.github.io/web/malavi/), the public
database of avian haemosporidian (malaria and related) parasite mtDNA cytochrome
*b* lineages.

> **Note on this release.** MalAvi is no longer hosted at a permanent web address (although it will be soon). The
> functions in the first version of `malaviR`  downloaded data from the MalAvi server / web address. Since that is no longer an option,
> I wanted to update `malaviR`, so that it would be useful again. Staffan Bensch (the creator and maintainer of MalAvi) has been emailing out the latest versions of the database to users, so what I've done is to bundle those database files with `malaviR`. For BLAST functions,  I'm using the code from a [Shiny app I developed](https://wimanet-science.github.io/web/malavi/blast/), which uses `DECIPHER` to create a BLAST-like functionality. There are also a few new helper functions here. This is something that I've wanted to do for a while, but have not had the time. So I decided to experiment with Claude Code
> (Opus 4.8) as a helper to more quickly rewrite the `malaviR` package. It was incredibly fast, and I think all the functions behave properly. But I'm still testing them and working on the description files. If you spot anything that should be changed or fixed, please do let me know either by opening an issue here or email me directly (vaellis@udel.edu).

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
```

## What's in the package

The bundled database is identified by its release date (e.g. the most recent one as of the update of this package is `2026-03-23`).

```r
library(malaviR)

malavi_version()        # release date of the bundled MalAvi database that the package is currently using
malavi_version("all")   # all bundled releases...I will keep some older versions of the database here and this is how you can see them
```

### Data tables and the alignment

```r
hosts_dat <- extract_table("Hosts and Sites Table")   # one of the five MalAvi tables you can access
aln   <- extract_alignment()                       # the cyt b alignment (it's stored as a DNAbin object)
plas  <- extract_alignment(genus = "Plasmodium")   # you can filter the alignment by parasite genus
```

### BLAST-like search

`blast_malavi()` searches a query sequence that the user provides against the MalAvi alignment
using a pre-built [DECIPHER](https://decipher.codes/) index (Requires DECIPHER >= v3.0 and Biostrings.)

```r
query <- gsub("-", "", paste(as.character(aln[1, ]), collapse = "")) # here we just select the first sequence in the MalAvi alignment to BLAST, but this is designed for you putting in your own sequences as a character string.
blast_malavi(query, top_n = 5) # here's your BLAST-like output
```

### Lineage and amplicon quality control checks

These are new functions that I'm still trying to get right. Please treat them as **experimental**.

`lineage_qc()` is a check of whether a MalAvi cyt b sequence (like one you get out of a Sanger sequence or elsewhere) looks plausible or not based on the larger database. It works by flagging strange or surprising features about a sequence including length (too short), gaps/ambiguities, stop codons (translated in frame under the protozoan mitochondrial
genetic code...code 4), distance to the nearest lineage in the MalAvi alignment, mutations at invariant or rarely varying
sites (it checks across the full MalAvi alignment), nonsynonymous/second-position/transversion changes,
and a sliding-window chimera checker (basically checking if part of the sequence matches one MalAvi lineage and another part matches a different sequence...could be indicative of a sequence chimera). Then it computes a rather arbitrary
`score` from 0 (suspicious) to 1 (expected based on the MalAvi alignment). This is creating flags or warnings for you to investigate, not telling you your sequence is necessarily wrong. `amplicon_qc()` is meant for working with cleaned amplicon sequence variants (ASVs) that you would get from sequencing the MalAvi region with short-read deep sequencing (2 x 300bp). After you go through a normal pipeline like dada2 or vsearch, my experience is that you will still have many ASVs and it's hard to know what's real. This will flag relatively rare ASVs and any that are genetically close (e.g., 1 bp different) from a very common ASV in the pool.

```r
seq <- paste(as.character(aln[1, ]), collapse = "")   # your own sequence here (should be aligned to MalAvi already)...this just grabs one of the existing MalAvi lineages
lineage_qc(seq)                                       # see the report and investigate any flags

variants <- data.frame(sequence = c(seq_a, seq_b), count = c(10000, 5)) # the sequences need to be real and 479bp length aligned to MalAvi...what you'd get out of your amplicon seq project
amplicon_qc(variants)                                 # check out the report...see if there are any suspicious ASVs in there that you will consider not analyzing further (there will be)
```

### Screening the whole database (studies vs. non-synonymous mutations)

Staffan Bensch (the MalAvi curator) pointed out to me that lineages reported by only a single
study may be more likely to carry non-synonymous changes in *cytb* than lineages found by
multiple studies. That pattern would be consistent with some single-study lineages being
sequencing errors. Two functions help you investigate this.

`lineage_studies()` counts how many distinct studies report each lineage (from the references
in the Hosts and Sites table — basically just a simple helper). `lineage_screen()` counts each
lineage's **singleton** substitutions: bases that the lineage *alone* carries at a well-covered
site, classified as synonymous, non-synonymous, or stop-codon-creating. A lineage reported by a
single study that *also* carries singleton non-synonymous changes is the kind of thing I'd be
wary of.

```r
library(dplyr)

lineage_studies() %>%                                  # one row per lineage
  select(lineage, n_studies, n_host_records, n_countries)

# do single-study lineages carry more singleton non-synonymous changes?
lineage_screen() %>%
  filter(in_hosts_table) %>%
  group_by(single_study = n_studies == 1) %>%
  summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#> # A tibble: 2 x 3
#>   single_study     n mean_nonsyn
#>   <lgl>        <int>       <dbl>
#> 1 FALSE         1179      0.0102
#> 2 TRUE          3571      0.0443   # ~4x higher among single-study lineages
```

You can sharpen the screen by restricting it to one parasite genus. The three haemosporidian
genera are deeply divergent, so judging a singleton against only its own genus (rather than the
pooled database) is more meaningful — and the pattern gets stronger:

```r
lineage_screen(genus = "Plasmodium") %>%
  filter(in_hosts_table) %>%
  group_by(single_study = n_studies == 1) %>%
  summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#> 1 FALSE          360      0.0139
#> 2 TRUE          1045      0.175    # single-study Plasmodium lineages stand out even more
```

You can also focus on a phylogenetic group of your choosing. Here we take SGS1 (a widespread
*P. relictum* lineage) and every lineage within 3 bp of it, then run the same comparison inside
that little clade:

```r
aln  <- extract_alignment()
m    <- toupper(as.character(aln))                       # one row per lineage
name <- sub("^[A-Za-z]_([^_]+).*$", "\\1", rownames(m))  # bare lineage names
sgs1 <- m[match("SGS1", name), ]

is_base  <- function(x) x %in% c("A", "C", "G", "T")     # number of base differences from SGS1
n_diff   <- apply(m, 1, function(s) sum(is_base(s) & is_base(sgs1) & s != sgs1))
sgs1_grp <- name[n_diff <= 3]                            # SGS1 + its close relatives (62 lineages)

lineage_screen() %>%
  filter(in_hosts_table, lineage %in% sgs1_grp) %>%
  group_by(single_study = n_studies == 1) %>%
  summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#> 1 FALSE           10      0
#> 2 TRUE            39      0.0769   # within the SGS1 group it's the single-study lineages again
```

### Repeated haplotypes ("synonymies")

Some incomplete MalAvi sequences (i.e., < 479 bp) match longer sequences but retain different
lineage names. This can inflate estimates of parasite lineage diversity as pointed out recently
([Tamayo-Quintero et al. 2025](https://doi.org/10.1371/journal.ppat.1012911)).
`synonymy_report()` quantifies the problem and identifies the overlapping lineages ("synonymies");
`clean_alignment()` produces a de-duplicated alignment, letting you choose which
name to keep (this function was present already in the old `malaviR`, but it has been rewritten). By default `clean_alignment()` keeps the most complete sequence (i.e., the longest ignoring Ns) in each group (`method = "overlap"`), but you can also choose which of the sequences in the overlap groups you want to keep (`keep = `). In the old version, you could randomly select from the overlapping haplotypes. That's probably not very useful, but in case you are nostalgic about it, I've included a random selector in the new version (`select = "random"`). See `?clean_alignment` for more details.

```r
synonymy_report()$summary             # how many names share a haplotype (i.e., shorter sequences that match longer sequences completely)
res <- clean_alignment(aln, method = "overlap")   # keeps most complete (i.e., longest) lineage per haplotype group (synonymy group)
```

### Host taxonomy

`match_taxonomy()` aligns MalAvi host species names to the taxonomy provided with the
[clootl](https://github.com/eliotmiller/clootl) package (i.e., eBird taxonomy) and it flags
names that don't match (synonyms, hybrids, `sp.`, etc.). The pre-built taxonomic key linking MalAvi to the `clootl` names is kept in the package as the `taxonomy` dataset. So you won't need to use `match_taxonomy` necessarily, but I keep it here because it should be useful updating the taxonomic key with each new update of MalAvi.

```r
match_taxonomy(c("Turdus merula", "Anas sp."))$key
data(taxonomy)
```

### Phylogenies

`sister_taxa()` returns the sister tips descending from nodes in a phylogeny. It was in the old version of `malaviR` and I've kept it here unchanged even though it more or less overlaps with functions in other phylogenetics R packges.

`clean_names()` strips the genus prefix from alignment tip labels. This was also in the older version of `malaviR` and can be useful for linking the alignment to the tables (alignment uses the genus prefix, tables do not).

## Updating the bundled database

When a new MalAvi release arrives (a `MalAvi_<date>.zip`), I will put it in
`data-raw/` and run `data-raw/process_release.R` to regenerate the bundled data
and BLAST index and then I'll push it to github.

## Citation

If you use `malaviR`, please cite the package and:

* Ellis VA, Bensch S (2018). Host specificity of avian haemosporidian parasites
  is unrelated among sister lineages but shows phylogenetic signal across larger
  clades. *International Journal for Parasitology* **48**: 897–902.
  <https://doi.org/10.1016/j.ijpara.2018.05.005>
* Bensch S, Hellgren O, Pérez-Tris J (2009). MalAvi: a public database of malaria
  parasites and related haemosporidians in avian hosts based on mitochondrial
  cytochrome *b* lineages. *Molecular Ecology Resources* **9**: 1353–1358.

The first citation (Ellis and Bensch 2018) was what I originally wrote `malaviR` for. The second (Bensch et al. 2009) is, of course, the citation for the MalAvi database.

```r
citation("malaviR")
```

## Issues

Please report problems at
<https://github.com/vincenzoaellis/malaviR/issues>.
