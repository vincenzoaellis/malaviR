# malaviR 1.1.1

`match_taxonomy()` now warns when it is given species names without `family` and `order`.

The family/order-constrained epithet match — the step that recovers genus reassignments
such as *Grus leucogeranus* to *Leucogeranus leucogeranus* — constrains its candidate pool
by the host's family, falling back to its order. Without them it cannot run, and it
previously skipped in silence. Called with a bare vector of names, 173 of MalAvi's own
2,339 host binomials came back `none` that the same function resolves when called with no
arguments, and nothing explained the difference.

No matching behaviour changed. Either pass `family =` and `order =` alongside `species`,
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
outranks a tidier neighbour list.

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
