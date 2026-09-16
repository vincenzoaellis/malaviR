# malaviR

<!-- badges: start -->
[![r-universe status](https://vincenzoaellis.r-universe.dev/badges/malaviR)](https://vincenzoaellis.r-universe.dev/malaviR)
[![GitHub release](https://img.shields.io/github/v/tag/vincenzoaellis/malaviR?label=release&sort=semver&color=blue)](https://github.com/vincenzoaellis/malaviR/releases)
<!-- badges: end -->

An R interface to [MalAvi](https://malavi-db.github.io/), the public
database of avian haemosporidian (malaria and related) parasite mtDNA cytochrome
*b* lineages. This is a rebuilt version of malaviR made with Claude Code. If you spot
anything that should be changed or fixed please let me know by opening an issue here
or emailing me directly at vaellis@udel.edu.

## Installation

```r
# install.packages("remotes")
remotes::install_github("vincenzoaellis/malaviR")
```

The core functions need only CRAN packages. Two optional features need extra
packages:

```r
# for blast_malavi() — local BLAST-like search (needs R >= 4.4)
# install.packages("BiocManager")
BiocManager::install(c("DECIPHER", "Biostrings"))
```

## What's in the package

The MalAvi database is identified by its release date.

```r
library(malaviR)

malavi_version()        # release date of the MalAvi database that the package is currently using
malavi_version("all")   # every release bundled with this version of the package
```

### Data tables and the alignment

```r
hosts_dat <- extract_table("Hosts and Sites Table")   # one of the five MalAvi tables you can access
aln   <- extract_alignment()                       # the cyt b alignment (it's stored as a DNAbin object...it's just one alignment, no more long vs. all seqs. But also see the sections on synonymies and ambiguous pairs below for some new controls for filtering the alignment)
plas  <- extract_alignment(genus = "Plasmodium")   # you can filter the alignment by parasite genus
```

### BLAST-like search

`blast_malavi()` searches a query sequence that the user provides against the MalAvi alignment
using a pre-built [DECIPHER](https://decipher.codes/) index (Requires DECIPHER >= v3.0 and Biostrings.)

```r
query <- gsub("-", "", paste(as.character(aln[1, ]), collapse = "")) # here we just select the first sequence in the MalAvi alignment, but this is designed thinking about unique sequences as character strings.
blast_malavi(query, top_n = 5) # here's your BLAST-like output
```

### Lineage quality control check

This is a new function that I'm still trying to get right. Please treat it as **experimental**.

`lineage_qc()` is a check of whether a MalAvi cyt b sequence (like one you get out of a Sanger sequence or elsewhere) looks plausible or not based on the larger database. It works by flagging strange or surprising features about a sequence including length (we expect 479bp), gaps/ambiguities, stop codons (translated in frame under the protozoan mitochondrial
genetic code...code 4), distance to the nearest lineage in the MalAvi alignment, mutations at invariant or rarely varying
sites, nonsynonymous/second-position/transversion changes,
and a sliding-window chimera checker (basically checking if part of the sequence matches one MalAvi lineage and another part matches a different sequence).

It reports all of that as **counts**, one row per sequence in `$summary`, so you can `rbind` a set of them into a table and sort or model it however suits your question. It does not combine them into a single number: a sequence can be unusual because it is wrong or because it is new, and no summary of these counts separates those two. Deciding is your job, and where you draw the line depends on what you are screening.

(Working with denoised amplicon sequence variants [ASVs] from short-read deep sequencing of the MalAvi region — quantifying lineages per sample, reconciling mixed infections, and flagging the rare 1-bp error variants of an abundant ASV — is handled by the companion **malaviASV** package [in development], which builds on `lineage_qc()`.)

```r
seq <- paste(as.character(aln[1, ]), collapse = "")   # your own sequence here (should be aligned to MalAvi already)
lineage_qc(seq)                                       # see the report and investigate any flags
```

One thing worth knowing: if you hand `lineage_qc()` a sequence that is *already* in MalAvi, it will find an exact match to itself and always come back `known_lineage`, no matter what the sequence looks like. To QC a lineage that's already in the database, hold it out of the reference (`lineage_qc(seq, reference = aln[-i, ])`), or use `lineage_screen()`, which is built to do exactly that for the whole database at once.

### Screening the whole database (studies vs. non-synonymous mutations)

Staffan Bensch pointed out to me that lineages reported by only a single
study may be more likely to carry non-synonymous changes in cytb than lineages found by
multiple studies. That pattern would be consistent with some single-study lineages being
sequencing errors. Two functions help you investigate this.

`lineage_studies()` counts how many distinct studies report each lineage (from the references
in the Hosts and Sites table; I know that's easy to do yourself, but it can be nice to have the helper, also it reminds you to do it). `lineage_screen()` counts each lineage's **singleton** substitutions: bases that the lineage *alone* carries, classified as synonymous, non-synonymous, or stop-codon-creating. A lineage reported by a single study that *also* carries singleton non-synonymous changes is the pattern Staffan was referring to.

```r
library(dplyr)

# do single-study lineages carry more singleton non-synonymous changes?
lineage_screen() %>%
  filter(in_hosts_table) %>%
  mutate(single_study = if_else(n_studies == 1, TRUE, FALSE)) %>%
  group_by(single_study) %>%
  summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#> # A tibble: 2 x 3
#>   single_study     n mean_nonsyn
#>   <lgl>        <int>       <dbl>
#> 1 FALSE         1179      0.0102
#> 2 TRUE          3571      0.0443   # ~4x higher among single-study lineages
```

You can also restrict to one parasite genus:

```r
lineage_screen(genus = "Plasmodium") %>%
  filter(in_hosts_table) %>%
  mutate(single_study = if_else(n_studies == 1, TRUE, FALSE)) %>%
  group_by(single_study) %>%
  summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#> 1 FALSE          360      0.0139
#> 2 TRUE          1045      0.175    
```

You can also focus on a phylogenetic group. For example, here we take SGS1 (a
*P. relictum* lineage) and every lineage within 3 bp of it, then run the same comparison inside
that clade. We measure genetic distances within *Plasmodium* only because it helps with speed.
`clean_names()` turns the alignment's genus-prefixed tip labels
into the regular lineage names, and `ape::dist.dna()` is used for calculating the genetic distances:

```r
library(ape)

aln <- extract_alignment(genus = "Plasmodium")
rownames(aln) <- clean_names(rownames(aln))                 

d <- dist.dna(aln, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
near_sgs1 <- names(which(d["SGS1", ] <= 3 / 479))  # SGS1 + lineages within 3 bp of it

lineage_screen() %>%
  filter(in_hosts_table, lineage %in% near_sgs1) %>%
  mutate(single_study = if_else(n_studies == 1, TRUE, FALSE)) %>%
  group_by(single_study) %>%
  summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#> 1 FALSE           10      0
#> 2 TRUE            35      0.0857   # similar pattern
```

### Repeated haplotypes ("synonymies")

Some incomplete MalAvi sequences (i.e., < 479 bp) match longer sequences but retain different
lineage names. This can inflate estimates of parasite lineage diversity as pointed out recently
([Tamayo-Quintero et al. 2025](https://doi.org/10.1371/journal.ppat.1012911)).
`synonymy_report()` quantifies the problem and identifies the overlapping lineages ("synonymies");
`clean_alignment()` produces a de-duplicated alignment, letting you choose which
name to keep (this function was present already in the old `malaviR`, but it has been rewritten). By default `clean_alignment()` keeps the most complete sequence (i.e., the longest ignoring Ns) in each group (`method = "overlap"`), but you can also choose which of the sequences in the overlap groups you want to keep (`keep = `). In the old version, you could randomly select from the overlapping haplotypes. That's probably not very useful, but in case you want it, I've included a random selector in the new version (`select = "random"`). See `?clean_alignment` for more details.

```r
synonymy_report()$summary # how many names share a haplotype (i.e., shorter sequences that match longer sequences completely)
res <- clean_alignment(aln, method = "overlap") # keeps most complete (i.e., longest) lineage per haplotype group (synonymy group)
```

Collapsing assumes that the shorter sequence matches the longer one at the positions at which its undefined (N or gap). That could be true, but of course we don't know for sure. Good to keep that in mind.

### Ambiguous pairs (a different problem from synonymies)

A **synonymy** (above) is one lineage *contained in* another: wherever the shorter one has a base, the longer one agrees, and the longer one fills in the rest. Keeping the most complete sequence means you are not throwing out any observed nucleotides. So for building an alignment, or for assigning a lineage name to a sequence, collapsing the synonymies and choosing the longer one seems reasonable.

An **ambiguous pair** also agrees wherever both lineages are determined, but here *each* lineage has a base where the *other* has an N or a gap. Neither is contained in the other, so there is no most-complete sequence to keep, and dropping either would discard a nucleotide from the overall dataset. An example should help:

```
SYNONYMY (containment) — synonymy_report() flags it, consider collapsing the lineages by selecting one (probably the longer one)
                                          alignment position 336
  P_SEIAUR01… (479/479 nucleotides known)         T        <- determined at pos 336
  P_CARCAR11  (478/479 nucleotides known)         N        <- undetermined at same position
  
  So, CARCAR11 is contained in SEIAUR01; it carries no base SEIAUR01 lacks.

AMBIGUOUS PAIR (mutually partial) — ambiguous_pairs() flags it
                                   pos 1     pos 404
  H_DENPEN02                         G         N      <- determined at 1, undetermined at 404
  H_PASILI01                         -         G      <- undetermined at 1, determined at 404
  
  So, each is determined where the other is not and neither contains the other.
```

```r
ambiguous_pairs()$summary      
head(ambiguous_pairs()$pairs)   # the pairs
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

`sister_taxa()` returns the sister tips descending from nodes in a phylogeny. It was in the old version of `malaviR` and I've kept it here unchanged even though it more or less overlaps with functions in other phylogenetics R packages.

`clean_names()` strips the genus prefix from alignment tip labels. This was also in the older version of `malaviR` and can be useful for linking the alignment to the tables (alignment uses the genus prefix, tables do not).

### Known problems in the MalAvi data

`malaviR` ships the MalAvi tables **verbatim** — I don't quietly correct them, because then the package and the database would disagree and you'd have no way to tell. But problems do turn up in the data, usually when someone is building something downstream and hits one. `malavi_issues()` is where the package remembers them.

```r
malavi_issues()                          # the known issues in the bundled release
malavi_issues(version = malavi_version())  # or a particular bundled release
```

It prints a short, plain list: a heading naming the release, then a title and one sentence per issue. The kinds of things it catches: a lineage whose alignment genus prefix disagrees with its `GENUS_NAME`, a `GENUS_NAME` that is the string `"N/A"` (which turns into a fake `N_` genus prefix), a lineage name that appears twice in the Grand Lineage Summary because it has been reported under two morphospecies, and lineages that are in the alignment but not the table or the other way round.

The part that matters is that each issue is **re-derived from the release you have loaded** rather than stored as text. The counts and lineage names you read are the ones found in your release, so the wording can't drift out of step with the data — and an issue that a future MalAvi release fixes stops being found and simply drops off the list, without me editing anything.

If you find something wrong in the data, please open an issue — I'd rather have it recorded here than have everyone rediscover it. Actual corrections go upstream, into a future MalAvi release.

## Important assumptions

A few important points:

- Most sequence functions use the MalAvi database reported by
  `malavi_version()`. To make that clear, that version (and the method used)
  is included as an attribute (`attr(x, "malavi_meta")`), so a saved result records how it was made.
- Pairwise-deletion distances (`pairwise_deletion_distance()`, `lineage_qc()`,
  `ambiguous_pairs()`) compare only positions where **both** sequences have a
  defined base (A/C/G/T). Ambiguous bases, `N`s, and gaps are ignored, so a
  `distance` should always be considered in the context of `n_comparable` bases.
- `N`s added by `frame_to_malavi()` mark positions that were **not observed**
  because they fall outside the sequence window; they are not biological
  ambiguity calls.
- `synonymy_report(method = "overlap")` and `clean_alignment(method = "overlap")`
  group sequences that are **compatible over their observed positions**; missing
  positions are unknown. Don't consider collapsed sequences as necessarily identical (we just don't know).
- `ambiguous_pairs()` returns mutually partial pairs (neither contains the other...see above).
  These should be **reviewed**, not necessarily collapsed.
- `lineage_qc()` returns counts, not a verdict. A sequence can be unusual because it is
  wrong or because it is new, and nothing in the output separates those two — deciding
  is your job, and where you draw the line depends on what you are screening. Also note
  that a lineage already in MalAvi will always match itself, so screening one
  is only informative if you hold it out of the reference (see above).
- The MalAvi tables are shipped **verbatim**. Known problems in them are reported by
  `malavi_issues()`, not silently patched.
- In `match_taxonomy()`, only `match_type == "exact"` is a perfect match.
  You should review the non-exact rows (`$differences`). I've done this, but extra sets of eyes are always helpful. Also note that `generic`/`none` are unresolved.

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
