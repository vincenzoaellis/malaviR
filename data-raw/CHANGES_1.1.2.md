# What changed in malaviR 1.1.2, and why

A review record for Vincenzo. Every change here came from the independent code and
biology review of 1.1.1 (commit `e597ac9`) and is one of the eight confirmed bugs in
section A of that review. Nothing in section B (interpretation) or section C
(housekeeping) is done yet.

The per-row taxonomy diff is also in `taxonomy_changes_1.1.2.csv` next to this file,
with the `ott_name` and `family` columns as well, if you want to work with it in R.

Reproduce the diff at any time by pulling the old crosswalk out of git:

```sh
git show e597ac9:data/taxonomy.rda > /tmp/taxonomy_1.1.1.rda
```
```r
old <- new.env(); load("/tmp/taxonomy_1.1.1.rda", envir = old)   # 1.1.1
new <- new.env(); load("data/taxonomy.rda",       envir = new)   # 1.1.2
i <- which(old$taxonomy$ebird_species != new$taxonomy$ebird_species)
```

Commits: `2526152`, `51fcab0`, `02af0b8`, `07e910b`, `da91459`.
`R CMD check` on a clean tarball: 0 errors, 0 warnings, 0 notes.

---

## 1. Host taxonomy crosswalk (`data/taxonomy.rda`)

17 of the 2,339 rows changed, 60 host records in total. Ten are real species
corrections; seven keep the same species and only gain a better-sourced label.
No row became unmatched: `none` is still 0.

### Species corrections (10 names, 39 host records)

| MalAvi host name | was | now | old match_type | new match_type | records |
|---|---|---|---|---|---|
| *Alethe poliocephala* | *Turdus poliocephalus* | *Chamaetylas poliocephala* | `reassigned:family` | `manual` | 5 | 
| *Hemispingus frontalis* | *Crithagra frontalis* | *Sphenopsis frontalis* | `reassigned:family` | `manual` | 10 | 
| *Kittacincla malabarica* | *Copsychus albiventris* | *Copsychus malabaricus* | `synonym:HowardMoore` | `manual` | 3 | 
| *Loxigilla violacea* | *Euphonia violacea* | *Melopyrrha violacea* | `reassigned:family` | `manual` | 1 | 
| *Phaethornis baroni* | *Metallura baroni* | *Phaethornis longirostris* | `reassigned:family` | `synonym:BirdLife-lump` | 1 | 
| *Thraupis cyanocephala* | *Chlorophonia cyanocephala* | *Sporathraupis cyanocephala* | `reassigned:family` | `manual` | 13 | 
| *Tiaris obscura* | *Akialoa obscura* | *Asemospiza obscura* | `reassigned:family` | `manual` | 3 | 
| *Turdoides gularis* | *Mixornis gularis* | *Argya gularis* | `reassigned:family` | `manual` | 1 | 
| *Vermivora pinus* | *Setophaga pinus* | *Vermivora cyanoptera* | `reassigned:family` | `manual` | 1 | 
| *Zosterops capensis* | *Zosterops pallidus* | *Zosterops virens* | `legacy` | `manual` | 1 | 

### Same species, label changed (7 names, 21 host records)

| MalAvi host name | eBird species (unchanged) | old match_type | new match_type | records |
|---|---|---|---|---|
| *Accipiter tachiro* | *Aerospiza tachiro* | `synonym:HowardMoore` | `synonym:IOC-lump` | 2 |
| *Camaroptera brevicaudata* | *Camaroptera brachyura* | `legacy` | `synonym:IOC-lump` | 2 |
| *Icterus chrysocephalus* | *Icterus cayanensis* | `legacy` | `synonym:BirdLife-lump` | 1 |
| *Pitohui ferrugineus* | *Pseudorectes ferrugineus* | `reassigned:family` | `manual` | 1 |
| *Pitohui incertus* | *Pseudorectes incertus* | `reassigned:family` | `manual` | 2 |
| *Tangara cayana* | *Stilpnia cayana* | `synonym:HowardMoore` | `synonym:BirdLife-lump` | 10 |
| *Tangara ruficervix* | *Chalcothraupis ruficervix* | `synonym:HowardMoore` | `synonym:BirdLife-lump` | 3 |

### Why the ten were wrong

Eight of them came in by one route. `.epithet_reassign()` recovers a genus transfer by
matching the specific epithet inside the host's **MalAvi family**, and `FAMILY_NAME` is
the least maintained field in the release. Where MalAvi still files a genus under an old
family, the candidate pool becomes the *current* clootl membership of that family name,
the true species is not in it because clootl has moved it elsewhere, and any lone bird in
the pool with the same epithet wins. *Tiaris obscura*, a Peruvian grassquit, was matched
to *Akialoa obscura*, an extinct Hawaiian honeycreeper, which put three infection records
on a tip no grassquit ever occupied. This is the same mechanism the June 2026 Copilot
review found for *Oriolus*; only that one row was patched at the time.

Seven of the eight were exact-epithet hits, and `.audit_taxonomy()` treated an exact
epithet as strong support and exempted it — so they never reached
`data-raw/taxonomy_audit.csv` for review. An exact epithet is no evidence at all when the
pool is the wrong family.

The two remaining corrections are concept errors rather than resolver failures:

- ***Kittacincla malabarica***. clootl carries the Howard & Moore name only on *Copsychus
  albiventris*, an Andaman endemic, because H&M treats *albiventris* as a subspecies. The
  synonym step found one unambiguous candidate and took it; the ambiguity guard could not
  help because the ambiguity is not in the data. The three MalAvi records are from China,
  which is White-rumped Shama *C. malabaricus*.
- ***Zosterops capensis***. Under Clements, *capensis* is a subspecies of Cape White-eye
  *Z. virens*; *Z. pallidus* is Orange River White-eye. The 2012 Jetz key lumped them
  under *pallidus* and the legacy bridge carried that forward. Confirmed by Vincenzo
  2026-09-04: *Z. virens* is the current name.

*Phaethornis baroni* is the one that fixed itself. It is a Clements subspecies of
Long-billed Hermit, and clootl lists it in the BirdLife cell
`"Phaethornis longirostris;Phaethornis baroni"` — invisible to a whole-cell string
compare, so the synonym step failed and the family pool matched *Metallura baroni*
instead. Splitting those cells resolves it correctly, and the `-lump` label says the
MalAvi host concept is narrower than the eBird species. Confirmed by Vincenzo 2026-09-04.

### What now prevents a recurrence

1. **Genus guard** in `.epithet_reassign()`. Before accepting a match that changes the
   genus, look up which families clootl files the *MalAvi* genus in, and decline the
   match if the candidate is outside them. Deliberately keyed on the genus's home in
   clootl, not on MalAvi's family label, because the label is the untrustworthy part.
   This refuses four of the eight on its own.

   It also refuses *Pitohui ferrugineus* and *Pitohui incertus* to *Pseudorectes*, which
   are **correct** — clootl files *Pitohui* in Oriolidae and *Pseudorectes* in
   Pachycephalidae, a genuine cross-family transfer the guard cannot tell from a false
   one. Both are carried by manual rows so they survive. Confirmed by Vincenzo 2026-09-04.

2. **Split synonym cells.** clootl's `IOC_name`, `Birdlife_name` and `H_M_name` hold
   semicolon-joined lists wherever another authority splits an eBird species (115, 302
   and 86 such cells in the bundled snapshot). 153 MalAvi host names sit inside one.
   Every name is now searched, and a name found only inside a joined cell gets a `-lump`
   label. That is what relabels the seven rows in the second table above, and it is
   better information than they had: `synonym:BirdLife-lump` tells you eBird lumps the
   concept, where `legacy` and `synonym:HowardMoore` did not.

3. **A `retired_genus` audit reason.** The genus guard is blind when clootl no longer uses
   the MalAvi genus, which is where *Hemispingus frontalis* hid. `.audit_taxonomy()` now
   lists those 33 genus-changing rows whether or not the epithet matched exactly. The
   reviewer's suggestion — drop the exact-epithet exemption outright — would have put 151
   rows in the audit, which is too many to eyeball; retired-genus is the sharp version of
   the same idea. Most of the 33 are the standard transfers (*Dendroica* to *Setophaga*,
   *Megalaima* to *Psilopogon*).

   One blind spot remains and is not automatable: *Vermivora pinus* to *Setophaga pinus*
   is a false match **inside the right family**, so no structural rule catches it. It is
   fixed by a manual row.

### Manual overrides

`data-raw/manual_taxonomy.csv` went from 83 to 94 rows. Each new row carries a `note`
saying which bird it is and what the resolver did instead. The nine older manual rows that
existed only because of the semicolon bug were **kept**, not removed: they are
hand-verified, and `manual` is a more stable label than one that depends on clootl not
changing.

The clootl snapshot itself did not change (11,167 rows, taxonomy year 2025).

---

## 2. `lineage_qc()` no longer calls a thin match exact

Distances use pairwise deletion, so a reference sharing no determined position with the
query also scores 0. That 0 was read as an exact match. Two confirmed cases on real data:

- a 180 bp partial query with three real substitutions came back `known_lineage`, nearest
  `L_EOPER01` at distance 0 — over **6** shared positions. Its true relative, at distance
  3 over 180, was not reported at all.
- a query of nothing but `N`s came back `known_lineage`, score 0.82, nearest `H_ABSUP01`,
  which is simply row 1 of the alignment.

With 3,338 of the 5,365 bundled sequences partial, and partial queries screened since
1.1.0, this was ordinary use rather than a corner case.

A distance-0 agreement now has to rest on at least 60% of the query's determined
positions. The floor is **relative to the query**: the old `.QC_MIN_COMPARABLE <- 300` was
a fixed count, which a 180 bp query can never reach against any reference, so it stopped
discriminating on exactly the queries that need it. An absolute floor was considered and
rejected — the least-covered bundled reference has 133 determined positions, so only a
relative rule fixes the short-query case without breaking a full-length query against a
partial reference.

New: flags `matches_known_lineage_over_short_overlap` and
`no_comparable_reference_overlap`, and `n_comparable` in `summary` so the overlap a
distance rests on is visible. An all-N query now reports `NA` for both nearest lineage and
distance.

Checked for regression: all 150 sampled bundled sequences still return `known_lineage`
when screened against the full reference.

## 3. `blast_malavi()` ranks by the alignment

`top_n` was applied to the `SearchIndex` k-mer score before anything was aligned. That
score is a filter, not a ranking, and it does not put an exact match first: for a full
479 bp SGS1 query the index ranks SGS1 **second**, so `blast_malavi(sgs1, top_n = 1)`
returned *P_PADOM07* at 99.776% and never reported the 100% self-match. For a 200 bp
partial the source lineage ranked 49th.

Every index hit is now aligned before the cut. Cost is about 0.7 s against roughly 4 s for
the search itself. A `ReferenceGapLength` column was added so that
`Matches + Mismatches + QueryGapLength + ReferenceGapLength` equals `AlignmentLength` for
a query carrying an insertion.

Verified under R 4.5.3 / DECIPHER 3.6, because system R 4.3 skips these tests: 20 of 20
random bundled lineages self-hit at `top_n = 1`, and the insertion query ranks its source
first.

**One limit remains and is documented, not fixed.** A short query is often identical to
many lineages over its own length — a 200 bp fragment of SGS1 matches about 40 lineages at
100% — and `top_n` then shows an arbitrary few of them. Ranking cannot separate hits the
alignment cannot separate.

## 4. Four smaller faults

- **`extract_table("all")`** returned the raw bundled tables while every single-table call
  whitespace-tidied, so the same table differed in **3,143 cells** depending on how it was
  asked for. The tidy also used `[[:space:]]`, which does not match a non-breaking space:
  22 `GENBANK_ACC` values kept an invisible leading U+00A0 that `trimws()` cannot remove
  and that breaks an accession join. Both are fixed; `"all"` is now identical to the
  single-table call for all five tables.
- **`sister_taxa()`** read the two children of a node out of the edge matrix by position,
  in four near-identical blocks, so the third and later children of a polytomy were
  dropped silently. It is now one loop over the children: a polytomy returns every
  descendant as its own sister group, and a node with one child says so instead of failing
  with "incorrect number of dimensions". Bifurcating behavior is unchanged.
- **The invariant-site penalty could never fire.** `.qc_score_site()` required the query
  base to have been *observed* before calling a change an invariant-site change, but at an
  invariant site only one base has ever been seen, so the test was always FALSE. The
  weight-4 penalty and the `N_changes_at_invariant_sites` flag were unreachable from 1.1.0
  to 1.1.1, while the mutations table reported the same base as an invariant-site change.
  Seven sites are invariant in the pooled profile today; more in a genus-restricted one.
- **`frame_to_malavi()`** aborted the whole call on a single `NA` element, and only
  trimmed the ends although the documentation said whitespace is stripped, so a 478 bp
  sequence with one interior space was reported off-length and discarded. `pad_char` is
  now checked to be a single character.

## 5. Housekeeping

Version 1.1.1 to **1.1.2** with a NEWS entry. `inst/CITATION` said 1.0.0 and now says
1.1.2. British spellings converted (behaviour, neighbour, summarises, labelled) and the
README "packges" typo fixed.

---

## Still open

- **Section B of the review** — eight biology and interpretation items. The headline is
  B1: leave-one-out on 60 curated bundled lineages puts 26% of them in
  `strong_warning`/`possible_error`, and every lineage more than 20 bp from its nearest
  neighbor is `possible_error`, so the `lineage_qc()` composite score is effectively a
  divergence measure rather than a plausibility measure.
- **Section C** — documentation and NEWS housekeeping, including the 1.1.0 NEWS entry's
  omissions and some stale README numbers.
- **Downstream**: the malavi_rebuild site carries this crosswalk. See RUNBOOK §6.
