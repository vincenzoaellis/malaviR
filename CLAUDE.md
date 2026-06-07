# malaviR — modernization project

Working notes for Claude Code and for Vincenzo. This file records **what we are doing
and why**, the **decisions** we have made, the **plan**, and a **living status log**.
It is auto-loaded by Claude Code each session. It is excluded from the built R package
(see `.Rbuildignore`).

Repo: https://github.com/vincenzoaellis/malaviR  (default branch: `master`)
Local working copy: `/mnt/biostore-all/Vellis/malaviR` (== `/mnt/ellisbiostore/malaviR`,
and reachable from `~/ellisbiostore/malaviR`).

---

## 1. Purpose / context

`malaviR` is an R interface to **MalAvi**, the public database of avian haemosporidian
(malaria & relatives) mtDNA cytochrome *b* lineages. The package was written ~2017–2019
and last touched Oct 2019. It is **broken today** because every web-facing function
scraped/downloaded from the MalAvi web server at IP `130.235.244.92`, which **no longer
exists** as a live, queryable site (open GitHub issues #6 and #7 are both caused by this).

We are modernizing the package so it works offline from **bundled database snapshots**
that Vincenzo receives directly from the MalAvi maintainers, plus rebuilding the BLAST
feature locally and refreshing the host taxonomy crosswalk against a modern avian
taxonomy (clootl / eBird).

## 2. Key decisions (confirmed with Vincenzo)

- **Ship the latest DB release only** inside the installed package (~0.5 MB of tables +
  alignment; +~3 MB for the prebuilt BLAST index). Readers are written to auto-discover
  any number of bundled releases, so more can be added later without code changes.
- **Archive each original release `.zip` in `data-raw/`** (git-tracked, NOT shipped to
  installs) for provenance.
- **Database "version" = the date in the file names** (e.g. `2026-03-23`). Simple.
- **Recurring update workflow:** drop a new `MalAvi_<date>.zip` into `data-raw/`, run
  `data-raw/process_release.R`, commit, push.
- **New version number: `1.0.0`** (from 0.2.0) — first stable, working, redesigned release.
- **Citation:** package should direct users to cite **Ellis VA, Bensch S. 2018. Host
  specificity of avian haemosporidian parasites is unrelated among sister lineages but
  shows phylogenetic signal across larger clades. Int. J. Parasitol. 48: 897–902.
  doi:10.1016/j.ijpara.2018.05.005** (first literature use of the package), alongside the
  original MalAvi paper (Bensch, Hellgren & Pérez-Tris 2009, Mol. Ecol. Resour. 9:
  1353–1358) and the package itself. Build an `inst/CITATION`; fix README's stale 0.1.0.
- **Coding style:** keep Vincenzo's existing style — **use `%>%` (magrittr), not `|>`**;
  base-R + tidyverse mix as before; roxygen2 docs; one dplyr verb per line after a pipe.
  Emphasis on **clarity and simplicity**.
- **New dependencies:** DECIPHER + Biostrings (Bioconductor, for offline BLAST) go in
  **Suggests** with a friendly runtime check so the core data functions still install
  without Bioconductor. clootl (CRAN, for taxonomy) — Suggests (used in a data-raw build
  step) unless a runtime function needs it.
- Package installs/builds/tests run **on a compute node** (currently srun job on
  biomix20), never the login node.

## 3. New data format (release = 6 files, dated)

A MalAvi release zip (`MalAvi_<date>.zip`) contains one folder with:

| File | Becomes |
|------|---------|
| `Hosts_and_Sites_<date>.xlsx`        | table `hosts_and_sites` (20 cols, ~18.5k rows) |
| `GrandLineageSummary_<date>.xlsx`    | table `grand_lineage_summary` (24 cols, ~5.4k rows) |
| `MorphoSpecies_<date>.xlsx`          | table `morpho_species` (5 cols) |
| `References_<date>.xlsx`             | table `references` (6 cols) |
| `VectorData_<date>.xlsx`             | table `vector_data` (6 cols) |
| `MalAvi_<date>.fas`                  | `alignment` (DNAbin, 5365 seqs × 479 bp) |

Note: the offline release has **5 tables** (vs. 9 scraped from the old website). The old
"Parasite Summary Per Host", "Table of Lineage Names", "Other Genes", and "Database
Summary Report" tables are not in the distributed release and are dropped.

### Storage layout in the package
- `inst/extdata/malavi_db_<date>.rds` — list(grand_lineage_summary, hosts_and_sites,
  morpho_species, references, vector_data, alignment). xz-compressed (~0.5 MB).
- `inst/extdata/malavi_blast_<date>.rds` — list(db = ungapped DNAStringSet,
  index = DECIPHER InvertedIndex). ~3 MB. Used only by `blast_malavi()`.
- `data-raw/MalAvi_<date>.zip` — archived original (provenance).
- `data-raw/process_release.R` — converts a zip → the two `.rds` files, prunes old ones.

## 4. Function plan (R/)

| Function | Status | Plan |
|----------|--------|------|
| `extract_table()` | rewrite | Offline reader: `extract_table(table, version = "latest")`. Reads bundled `malavi_db_*.rds`. Keeps the bad-name `stop()` (issue #1). Drops to the 5 available tables. `"all"` still supported. |
| `extract_alignment()` | rewrite | Offline reader: `extract_alignment(version = "latest", genus = c("all","Plasmodium","Haemoproteus","Leucocytozoon"), include_unknown=)`. Returns bundled DNAbin, optionally subset by parasite genus (folds in issue #3). |
| `malavi_version()` | rewrite | No more scraping a dead JS file. Returns the date(s) of the bundled release(s). |
| `malavi_versions()` | NEW | Lists all DB versions bundled in the install. |
| `blast_malavi()` | rewrite | Offline BLAST-like search using DECIPHER (ported from the Shiny app in `MalAviBLAST/app.R`): `SearchIndex()` against the prebuilt index, then `AlignPairs()`; returns a tidy hit table (Lineage, ProportionMatch, PercentMatch, AlignmentLength, Matches, Mismatches, Score, QueryGapLength, ReferenceLineageLength). Uses bundled `malavi_blast_*.rds`. Fixes issues #7. Requires DECIPHER/Biostrings (Suggests + runtime check). |
| `clean_alignment()` | redesign | Currently picks ONE lineage per identical-haplotype group **at random** (`sample_n(1)`). Redesign to (a) return a clear **synonymy table** of identical-sequence lineage groups, (b) give the user **control** over which lineage to keep (e.g. keep-by-rule or user-supplied choice) rather than random. Acknowledge Tamayo-Quintero et al. 2025 (the paper documenting this exact problem). |
| `clean_names()` | keep | Minor doc/style polish only. |
| `sister_taxa()` | keep | **Verified correct** on a known tree (two-clade, clade+tip, tip+clade, two-tip, and multi-node cases). Multi-node support already satisfies issue #2. Doc/style polish only. |

## 5. Taxonomy crosswalk (replaces 2017 `taxonomy` data object)

Old `data/taxonomy.rda`: 1457×3 (`species`, `Jetz.species`, `match`) mapped to the 2012
Jetz/birdtree.org phylogeny — **out of date**.

New plan: build a crosswalk from MalAvi host species (from `hosts_and_sites`
GENUS_NAME + SPECIES_NAME) to **clootl** species (modern, eBird-aligned avian taxonomy /
phylogeny). Provide:
- an updated bundled `taxonomy` data object (MalAvi host name → clootl/eBird name + match
  status), and
- a function that aligns MalAvi host names to clootl species and **reports the key plus
  any differences/mismatches** for the user to review.

## 6. Docs

- **Remove** the out-of-date vignette `vignettes/Getting_Started_with_malaviR.Rmd` (it
  downloads from the dead server during build → causes install failure, issue #6).
- **Refresh `README`** (landing page) for the new offline workflow.
- **Write ONE new basic vignette** covering normal use of the functions — must be fully
  **network-free** so the package builds/installs cleanly.

## 7. GitHub issue mapping

| Issue | State | Addressed by |
|-------|-------|--------------|
| #7 blast_malavi connection timeout | open | Offline DECIPHER BLAST (§4) |
| #6 install fails from GitHub (vignette downloads from dead server) | open | Remove vignette + offline data + network-free new vignette (§6) |
| #5 malavi_version regex fix | closed | `malavi_version` replaced by offline reporter |
| #4 add malavi_version | closed | preserved (offline form) |
| #3 subset extract_alignment by genus | closed | folded into new `extract_alignment` (§4) |
| #2 multiple nodes in sister_taxa | closed | preserved & verified working |
| #1 extract_table error on bad name | closed | preserved in rewrite |

## 8. Build / test / push workflow

- Bump `DESCRIPTION` Version; update `RoxygenNote` (currently ancient 6.1.1 → 7.x);
  add/adjust Imports & Suggests.
- `devtools::document()` to regenerate `NAMESPACE` + `man/`.
- `devtools::check()` on the compute node; rewrite tests in `tests/testthat/`.
- Commit in reviewable chunks with clear messages; push to `master`.
  (`gh` CLI not installed — use plain `git` over HTTPS; needs a PAT.)

## 9. Environment notes

- On compute node biomix20 (srun job), 2 CPU / 16 G. R 4.3.3.
- Biostrings: installed. DECIPHER: NOT yet. clootl: NOT yet. → install on the node.
- readxl + openxlsx: installed (used by process_release.R).

---

## STATUS LOG (newest first)

- **2026-06-07 (cont.)** — Implementation well underway:
  - Data layer DONE: `inst/extdata/malavi_db_2026-03-23.rds` + `malavi_blast_2026-03-23.rds`;
    `data-raw/process_release.R` validated end-to-end under R 4.5 (rebuilds index fresh).
  - Installed deps: system R 4.3.3 userlib has clootl + DECIPHER 2.30; created conda env
    `~/conda_envs/malaviR` (R 4.5.3 + DECIPHER 3.6) for index builds / blast testing,
    because `IndexSeqs`/`SearchIndex`/`AlignPairs` need DECIPHER >= 3.0 (R >= 4.4).
  - Functions rewritten/added & tested: extract_table, extract_alignment (+genus, issue #3),
    malavi_version, malavi_versions, blast_malavi (offline DECIPHER; verified self-hit 100%
    under R 4.5), clean_names (+docs), sister_taxa (+docs), clean_alignment (synonymy table,
    user-controlled keep, `method = "strict"|"overlap"`), match_taxonomy + clootl_taxonomy_version
    (clootl 2025 snapshot in R/sysdata.rda; data/taxonomy.rda rebuilt = 2339 hosts), and
    synonymy_report (quantifies synonymies paper-style, flags partials, by-genus).
  - DESCRIPTION -> 1.0.0; Imports trimmed to ape/dplyr/magrittr/utils; Suggests add
    DECIPHER/Biostrings/clootl/readxl. RoxygenNote 7.3.1. NAMESPACE/man regenerated.
  - Tests: 35 pass / 3 skip (blast skips without DECIPHER>=3) under system R 4.3.
  - TODO: README (incl. AI-acknowledgment sentence + Ellis&Bensch 2018 citation), one
    network-free vignette, inst/CITATION, remove old vignette, devtools::check, then
    PAUSE for review before pushing to GitHub.

- **2026-06-07** — Cloned repo, examined all source/data/tests/issues. Inspected the
  new release zip (`MalAvi_2026-03-23.zip`) and the Shiny BLAST app
  (`MalAviBLAST/app.R` + prebuilt `malavi_data.rds`, valid: db 5365 seqs + DECIPHER
  index). Verified `sister_taxa` correct. Confirmed all decisions above. Wrote this
  CLAUDE.md. **Next:** get plan approval, then install DECIPHER + clootl on the node and
  begin implementation starting with the data layer (`data-raw/process_release.R`).
