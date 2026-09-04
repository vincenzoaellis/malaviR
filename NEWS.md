# malaviR 1.2.0

**`lineage_qc()` no longer returns a plausibility score.** The `score` element, the
`exp(-penalty/10)` mapping and its weights, and the four calls that were score cutoffs
(`plausible_new_lineage`, `review`, `strong_warning`, `possible_error`) are all removed.

The reason is that the score did not measure what it claimed to. Holding out 60 curated
lineages already in MalAvi and re-screening them put 26% of them in `strong_warning` or
`possible_error`, and every lineage more than 20 bp from its nearest neighbor came back
`possible_error` -- including L_PIRIE01, a described *Leucocytozoon* species, at a score of
0. No penalty term was normalized by distance, so the composite was in effect a divergence
measure, and divergence is the one property a genuinely new lineage has. It could not
separate an artifact from a discovery while its call names implied it could.

Everything it was built from is kept, because each piece is a checkable fact. `summary` is
now one row carrying all of them -- `n_mutations` and its breakdown
(`n_nonsynonymous`, `n_second_position_changes`, `n_transversions`), how unusual the
query's bases are for their sites (`n_invariant_site_changes`, `n_bases_never_observed`,
`n_rare_site_bases`), `n_stop_codons`, `n_comparable` and `chimera_delta` -- so `rbind`
across a set of sequences gives a table to sort, filter and model. Judging them is left to
the user, because a within-sample ASV set, a batch of new deposits, and a re-check of the
database itself do not want the same thresholds. `counts` is unchanged for code written
against the old shape.

`call` keeps only the outcomes that rest on a single fact: `known_lineage` (an exact match
over enough of the query), `contains_stop_codon` (renamed from
`invalid_or_strong_warning`), `possible_chimera`, `invalid_sequence`, and
`no_exact_match` as the residual -- which is a statement about the reference, not a
verdict on the query.

# malaviR 1.1.2

Fixes from an independent code and biology review of 1.1.1. Two of them change
answers users act on, so results from 1.1.1 and 1.1.2 are not comparable.

**A match is only "exact" over enough shared sequence.** Distances are computed with
pairwise deletion, so a reference sharing no determined position with the query also
scores 0. That 0 was read as an exact match, and `lineage_qc()` reported the query as a
known lineage on the strength of it: a 180 bp partial query with three real substitutions
came back `known_lineage` against a reference sharing 6 positions with it, and a query of
nothing but `N`s matched the first row of the alignment. With 3,338 of the 5,365 bundled
sequences partial, and partial queries screened since 1.1.0, this was ordinary use.
A distance-0 agreement now needs to rest on at least 60% of the query's determined
positions, and references below that floor rank last whatever their distance. The floor is
relative to the query rather than the fixed 300 positions used before, which a 180 bp query
could never reach against any reference. Two new flags say what used to be called an exact
match — `matches_known_lineage_over_short_overlap` and `no_comparable_reference_overlap` —
and `summary` gains `n_comparable`, the overlap the distance was measured over.

**`blast_malavi()` ranks by the alignment, not by the index score.** `top_n` was applied to
the `SearchIndex` k-mer score before anything was aligned, and that score does not put an
exact match first: `blast_malavi(<full SGS1>, top_n = 1)` returned `P_PADOM07` at 99.776%
and never reported the 100% self-match. Every hit is now aligned before the cut. A new
`ReferenceGapLength` column makes the counts reconcile for a query carrying an insertion.

**Ten host species in `taxonomy` were matched to the wrong bird.** The family/order epithet
step builds its pool from MalAvi's `FAMILY_NAME`, so where MalAvi files a genus under an old
family a lone same-epithet member of the current family can win — *Tiaris obscura*, a
Peruvian grassquit, had been matched to *Akialoa obscura*, an extinct Hawaiian honeycreeper.
A genus-changing match is now checked against the families clootl files the MalAvi genus in.
Separately, clootl's synonym columns hold semicolon-joined lists wherever another authority
splits an eBird species, and whole-cell comparison could never see inside them; 153 MalAvi
host names sit in one. Those are now searched, and a name found only inside a joined cell
gets a `-lump` label, because the MalAvi host concept is then narrower than the eBird species
it maps to.

Smaller fixes: `extract_table("all")` skipped the whitespace tidy every single-table call
applies (3,143 cells), and the tidy left non-breaking spaces in place; `sister_taxa()`
dropped the third and later children of a polytomy; the invariant-site penalty in
`lineage_qc()` could never fire; `frame_to_malavi()` aborted on an `NA` element and did not
strip interior whitespace.

# malaviR 1.1.1

`match_taxonomy()` now warns when it is given species names without `family` and `order`.

The family/order-constrained epithet match — the step that recovers genus reassignments
such as *Grus leucogeranus* to *Leucogeranus leucogeranus* — constrains its candidate pool
by the host's family, falling back to its order. Without them it cannot run, and it
previously skipped in silence. Called with a bare vector of names, 173 of MalAvi's own
2,339 host binomials came back `none` that the same function resolves when called with no
arguments, and nothing explained the difference.

No matching behavior changed. Either pass `family =` and `order =` alongside `species`,
or call `match_taxonomy()` with no arguments to use MalAvi's host list, which carries them.
For MalAvi host names the resolved crosswalk is already shipped as `malaviR::taxonomy`.

# malaviR 1.1.0

Changes to `lineage_qc()`. All four were found on 2026-08-20 by one real community
submission whose sequences were slightly awkward rather than clean 479 bp barcodes.
Results from 1.0.0 and 1.1.0 are not comparable, which is why the version moved.

## A partial barcode is now screened instead of rejected

Any length other than 479 bp previously returned `invalid_sequence` with every metric
`NA` — and the reading-frame diagnosis written to explain that situation sat below the
length gate and was unreachable. Since 3,340 of MalAvi's 5,368 lineages cover only part
of the barcode window, and a primer-trimmed amplicon (478 or 476 bp) or a one-primer
read is what submitters commonly send, the function could not screen the majority of
real queries. A partial read of a lineage MalAvi *already held* came back invalid rather
than as an exact match.

A shorter query is now placed into the 479 bp frame and screened in full. Padded
positions score as unknown, so nothing is invented; the metrics simply cover fewer
positions. A query that cannot be placed convincingly is still `invalid_sequence`, and a
query *longer* than the frame is still refused, because placing it would discard real
bases. New flag: `placed_in_malavi_frame`.

## The nearest lineage is ranked by rate of mismatch, not by count

Ranking on the raw count let a reference overlapping the query in few positions win by
having less to disagree over. A real candidate was reported nearest to a reference with
23 mismatches over only 133 comparable positions (82.7% identity, and the least-covered
sequence in the bundled alignment), ahead of its true relative at 38 mismatches over 477
(92.0%) — naming a *Plasmodium* as the closest relative of a *Haemoproteus*.

An exact match still wins whatever it covers: never reporting a known lineage as new
outranks a tidier neighbor list.

## The chimera screen now tests for chimeras

`possible_chimera_or_mixed_template_pattern` fired on **23.2%** of ordinary lineages
(leave-one-out over the release). `parent_switches >= 2` was met by 98.4% of them, and
the approximate two-parent distance sat near zero, so the call collapsed to "the nearest
single lineage is 3 bp or more away" — a divergence measure, which is the one property a
genuinely new lineage has.

Replaced with a two-parent breakpoint scan: for every breakpoint, the best prefix parent
plus the best suffix parent, compared against the best single parent. Calibrated on 200
real lineages (leave-one-out) as negatives and 200 synthetic spliced chimeras as
positives — see `data-raw/chimera_v2_eval.R`. At the shipped threshold: **91.5% of true
chimeras detected at 2.0% false positives.**

New guards: the two parents must differ from each other (`chimera_min_parent_distance`,
default 5 bp), each segment must be at least `chimera_min_segment` (default 60 bp), and a
reference must cover most of the query. `chimera_delta_threshold` moves from 3 to 8.
`chimera_min_parent_switches` is still reported but no longer takes part in the call.

## The frame-shift message names the indel

A shifted reading frame has two causes and the test cannot separate them: a short
amplicon padded on the wrong end (a handling artifact — the sequence is fine), or an
indel (a sequencing error — the sequence needs re-reading). The message named only the
first, and confidently. It now gives both, and offers the one thing that discriminates:
stop codons surviving correct placement point to an indel.
