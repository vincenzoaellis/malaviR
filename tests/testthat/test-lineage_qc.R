## A tiny in-frame reference alignment (length 12 = 4 codons) so translation,
## site profiling, and nearest-lineage logic can be checked deterministically.
##   ref1: ATG TTT GGG CCC  -> M F G P
##   ref2: ATG TTC GGG CCC  -> M F G P (synonymous variant at position 6)
##   ref3: ATG TTT GGA CCC  -> position-9 variant
make_ref <- function() {
  ref1 <- strsplit("atgtttgggccc", "")[[1]]
  ref2 <- strsplit("atgttcgggccc", "")[[1]]
  ref3 <- strsplit("atgtttggaccc", "")[[1]]
  ape::as.DNAbin(rbind(ref1 = ref1, ref2 = ref2, ref3 = ref3))
}

test_that("build_malavi_site_profile summarizes sites correctly", {
  prof <- build_malavi_site_profile(make_ref())
  expect_equal(nrow(prof), 12)
  expect_equal(prof$codon_position, rep(1:3, 4))

  ## position 1 is invariant (all A); position 6 varies (T in ref1/ref3, C in ref2)
  expect_true(prof$invariant[1])
  expect_equal(prof$major_base[1], "A")
  expect_false(prof$invariant[6])
  expect_equal(sort(strsplit(prof$observed_alleles[6], "")[[1]]), c("C", "T"))
})

test_that("lineage_qc calls an exact match a known lineage", {
  qc <- lineage_qc("atgtttgggccc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_s3_class(qc, "malavi_lineage_qc")
  expect_equal(qc$call, "known_lineage")
  expect_equal(qc$summary$call, "known_lineage")
  expect_equal(qc$summary$nearest_distance, 0)
  expect_true("exact_match_to_known_lineage" %in% qc$flags)
})

test_that("lineage_qc reports distance to the nearest lineage", {
  ## one base off ref1 (position 12 C->A) -> nearest distance 1
  qc <- lineage_qc("atgtttgggcca", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_equal(qc$summary$nearest_distance, 1)
  expect_true("near_known_lineage" %in% qc$flags)
})

test_that("lineage_qc detects a stop codon under genetic code 4", {
  ## TAA at codon 2 is a stop under code 4
  qc <- lineage_qc("atgtaagggccc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_true("contains_stop_codon" %in% qc$flags)
  expect_equal(qc$summary$n_stop_codons, 1)
  expect_equal(qc$call, "invalid_or_strong_warning")

  ## TGA is tryptophan (not a stop) under code 4 -> no stop flag
  qc2 <- lineage_qc("atgtgagggccc", make_ref(),
                    expected_length = 12, chimera_check = FALSE)
  expect_false("contains_stop_codon" %in% qc2$flags)
})

test_that("lineage_qc returns an invalid_sequence result for wrong length", {
  qc <- lineage_qc("atgtttgggcc", make_ref(),    # 11 bp
                   expected_length = 12, chimera_check = FALSE)
  expect_equal(qc$call, "invalid_sequence")
  expect_equal(qc$score, 0)
  expect_true(any(grepl("wrong_length", qc$flags)))
})

test_that("lineage_qc honors a user-set rare_base_frequency threshold", {
  ## a base that is real but uncommon at its site should be counted as "rare"
  ## only once the rare_base_frequency cutoff is raised above its frequency
  q <- "atgttcgggccc"   # the position-6 minority base (C)
  n_low  <- lineage_qc(q, make_ref(), expected_length = 12,
                       rare_base_frequency = 0.001, chimera_check = FALSE)$counts[["n_rare_site_bases"]]
  n_high <- lineage_qc(q, make_ref(), expected_length = 12,
                       rare_base_frequency = 0.99, chimera_check = FALSE)$counts[["n_rare_site_bases"]]
  expect_gt(n_high, n_low)
})

test_that("lineage_qc keeps full category counts and details on request", {
  qc <- lineage_qc("atgtttgggcca", make_ref(), expected_length = 12,
                   chimera_check = FALSE)
  expect_true(all(c("n_invariant_site_changes", "n_rare_site_bases",
                    "n_transversions") %in% names(qc$counts)))
  expect_null(qc$translation)            # details off by default

  qc2 <- lineage_qc("atgtttgggcca", make_ref(), expected_length = 12,
                    chimera_check = FALSE, details = TRUE)
  expect_false(is.null(qc2$translation))
})

test_that("an exact match outranks an N-containing near-twin at distance 0", {
  ## Regression for the .qc_nearest tie-break bug: distance is computed with
  ## pairwise deletion (a reference N is skipped, not a mismatch), so a reference
  ## that is identical to the query except for one N ties a true exact match at
  ## distance 0. The N-twin is placed FIRST here so that the old arbitrary
  ## tie-break (alignment order) would wrongly report it; the fix breaks ties by
  ## most comparable positions, so the genuine full-length exact match wins.
  ## (Mirrors the real case: a perfect P_SEIAUR01 ASV was mislabeled P_CARCAR11,
  ## whose MalAvi reference carries a single N.)
  ntwin <- strsplit("atgtttgggccn", "")[[1]]   # identical to query but N at pos 12
  clean <- strsplit("atgtttgggccc", "")[[1]]   # the genuine exact match
  ref   <- ape::as.DNAbin(rbind(ntwin = ntwin, clean = clean))

  qc <- lineage_qc("atgtttgggccc", ref, expected_length = 12, chimera_check = FALSE)
  expect_equal(qc$summary$nearest_distance, 0)
  expect_true("exact_match_to_known_lineage" %in% qc$flags)
  expect_equal(qc$summary$nearest_lineage, "clean")   # not the N-twin
})

test_that("the hardcoded genetic code is NCBI table 4, not table 5", {
  ## Regression: ATA/AGA/AGG were previously given their invertebrate
  ## mitochondrial (table 5) values M/S/S. Table 4 is the standard code with the
  ## single change TGA = W, so those three keep their standard meanings. Getting
  ## them wrong does not affect stop counts (none of the three is ever a stop),
  ## but it does corrupt reported translations and the synonymous /
  ## nonsynonymous classification of substitutions among them.
  code <- .qc_genetic_code_4()
  expect_equal(length(code), 64L)
  expect_equal(unname(code[["TGA"]]), "W")   # the one table-4 difference
  expect_equal(unname(code[["ATA"]]), "I")   # table 5 would say M
  expect_equal(unname(code[["AGA"]]), "R")   # table 5 would say S
  expect_equal(unname(code[["AGG"]]), "R")   # table 5 would say S
  ## TAA and TAG are the only stops under table 4
  expect_equal(sort(names(code)[code == "*"]), c("TAA", "TAG"))
})

test_that("the hardcoded genetic code matches Biostrings table 4 codon for codon", {
  ## The authoritative check, run whenever Biostrings is available: compare all
  ## 64 codons against NCBI table 4 rather than re-typing the table by hand.
  skip_if_not_installed("Biostrings")
  reference_code <- Biostrings::getGeneticCode("4")
  ours <- .qc_genetic_code_4()
  expect_equal(ours[sort(names(ours))],
               reference_code[sort(names(reference_code))],
               ignore_attr = TRUE)
})

test_that("an AGA/AGG substitution is scored synonymous, not nonsynonymous", {
  ## Consequence of the table-5 bug: AGA and AGG are both arginine under table 4,
  ## so a change between them is synonymous. Under the old (table 5) values they
  ## were both serine, which happened to agree here -- but AGA vs CGA (both R
  ## under table 4, R vs R under table 5 too) is not a discriminating case, so we
  ## use AGA -> CGA, which the old table scored S -> R, i.e. nonsynonymous.
  ref <- ape::as.DNAbin(rbind(r1 = strsplit("atgagaggg", "")[[1]]))
  qc  <- lineage_qc("atgcgaggg", ref, expected_length = 9, chimera_check = FALSE)
  expect_equal(nrow(qc$mutations), 1L)
  expect_equal(qc$mutations$nearest_aa, "R")
  expect_equal(qc$mutations$query_aa, "R")
  expect_true(qc$mutations$synonymous)
  expect_equal(qc$summary$n_nonsynonymous, 0L)
})

test_that("lineage_qc flags a likely frame shift when another frame is stop-free", {
  ## A right-length query that is stop-ridden in frame 1 but stop-free in frame 2
  ## or 3 has a shifted reading frame. Two causes do that -- a short amplicon
  ## padded on the wrong end (handling), or an indel (sequencing error) -- and
  ## this test does not distinguish them because the diagnostic cannot either;
  ## it only checks that the message names both. Here frames 1 and 2 both carry
  ## a stop and frame 3 does not, so the diagnosis should name frame 3.
  ##   frame 1: TAA CTA AGG GCC -> *  L R A   (stop)
  ##   frame 2: AAC TAA GGG CC  ->  N *  G    (stop)
  ##   frame 3: ACT AAG GGC C   ->  T K  G    (clean)
  qc <- lineage_qc("taactaagggcc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_true("contains_stop_codon" %in% qc$flags)
  expect_true("possible_frame_shift_check_padding" %in% qc$flags)
  expect_match(qc$message, "frame 3")
  expect_match(qc$message, "padded on the wrong end")
  ## An indel is one of the commonest faults sent back for resequencing, so the
  ## message must say the word rather than steering the reader toward the benign
  ## explanation.
  expect_match(qc$message, "indel")
})

test_that("lineage_qc does not claim a frame shift when every frame has a stop", {
  ## Stops in all three frames means the sequence is genuinely bad, not shifted,
  ## so the frame-shift flag and its message must stay absent.
  ##   frame 1: TTA GTT AGT TAG,  frame 2: TAG TTA GTT AG,  frame 3: AGT TAG TTA G
  qc <- lineage_qc("ttagttagttag", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_true("contains_stop_codon" %in% qc$flags)
  expect_false("possible_frame_shift_check_padding" %in% qc$flags)
  expect_null(qc$message)
})

test_that("a correctly framed query never gets the frame-shift flag", {
  ## The check only ever runs on a query that already has a frame-1 stop, so a
  ## clean sequence cannot pick the flag up.
  qc <- lineage_qc("atgtttgggccc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_false("possible_frame_shift_check_padding" %in% qc$flags)
  expect_null(qc$message)
})

test_that("a real MalAvi lineage padded on the wrong end is diagnosed as shifted", {
  ## The realistic case from the field: a primer-trimmed haemosporidian ASV is
  ## 478 bp covering frame positions 2-479, so padding it with N at the 3' end
  ## instead of the 5' end yields a 479 bp query that is out of frame. It should
  ## be rejected AND told why, rather than being reported as merely divergent.
  aln <- extract_alignment()
  ## use the first complete, unambiguous lineage: many MalAvi entries are partial
  ## (gaps/Ns), and those would confound the frame reading
  seqs <- toupper(apply(as.character(aln), 1, paste, collapse = ""))
  complete <- seqs[nchar(seqs) == 479 & !grepl("[^ACGT]", seqs)]
  skip_if_not(length(complete) > 0)
  full <- unname(complete[1])

  shifted <- paste0(substr(full, 2, 479), "N")   # dropped base 1, padded at 3'
  qc <- lineage_qc(shifted, allow_ambiguity = TRUE, chimera_check = FALSE)
  expect_true("possible_frame_shift_check_padding" %in% qc$flags)
  expect_match(qc$message, "frame 3")

  ## padding the same 478 bp on the correct (5') end puts it back in frame
  correct <- paste0("N", substr(full, 2, 479))
  qc2 <- lineage_qc(correct, allow_ambiguity = TRUE, chimera_check = FALSE)
  expect_false("possible_frame_shift_check_padding" %in% qc2$flags)
})

test_that("lineage_qc rejects unsupported genetic codes", {
  expect_error(lineage_qc("atgtttgggccc", make_ref(), genetic_code = 1),
               "genetic_code = 4")
})

test_that("a lineage already in the reference always matches itself (documented caveat)", {
  ## ?lineage_qc warns that screening a lineage that is in MalAvi is silently
  ## uninformative, because the lineage is in its own reference. Lock that
  ## behavior down so the documentation cannot drift from it.
  aln <- extract_alignment()
  ## SGS1's alignment name carries a morphospecies suffix, so match on the
  ## cleaned lineage rather than assuming the two-token form
  i    <- which(clean_names(rownames(aln)) == "SGS1")
  expect_length(i, 1)
  name  <- rownames(aln)[i]
  query <- paste(as.character(aln[i, ]), collapse = "")

  qc <- lineage_qc(query)
  expect_equal(qc$call, "known_lineage")
  expect_equal(qc$summary$nearest_lineage, name)
  expect_equal(qc$summary$nearest_distance, 0)

  ## holding it out is what makes the screen informative: the nearest lineage
  ## is now some other sequence
  loo <- lineage_qc(query, reference = aln[-i, ])
  expect_false(loo$summary$nearest_lineage == name)
})

test_that("the nearest lineage is ranked by rate of mismatch, not by count", {
  ## Measured on a real submission, 2026-08-20: a candidate was reported nearest to a
  ## reference with 23 mismatches over only 133 comparable positions (82.7% identity, and
  ## the least-covered sequence in the alignment), ahead of its true relative at 38
  ## mismatches over 477 (92.0%). That named a Plasmodium as the closest relative of a
  ## Haemoproteus. Ranking on the raw count is what did it.
  refcode <- rbind(
    ## a thinly covered reference: only 150 positions known, 20 of them disagreeing
    c(rep(1L, 130), rep(2L, 20), rep(0L, 329)),
    ## a fully covered reference: all 479 known, 40 disagreeing
    c(rep(1L, 439), rep(2L, 40))
  )
  qcode <- rep(1L, 479)
  got <- malaviR:::.qc_nearest(qcode, refcode, c("thin", "full"), top_n = 2L)

  expect_equal(got$lineage[1], "full")
  expect_gt(got$distance[1], got$distance[2])          # more mismatches, still nearer
  expect_gt(got$n_comparable[1], got$n_comparable[2])
})

test_that("an exact match outranks a better-covered near match when it covers the query", {
  ## Load-bearing: exact_match_to_known_lineage is decided from distance[1], and never
  ## reporting a known lineage as new outranks a tidier neighbor list. That protection
  ## is kept -- but only for an exact match that actually covers the query. Here the
  ## exact reference is determined at 400 of the query's 479 positions, well above the
  ## 60% floor, so it still beats a fully-overlapping reference that differs by one base.
  refcode <- rbind(
    c(rep(1L, 400), rep(0L, 79)),                      # exact over 400 of 479 positions
    c(rep(1L, 478), 2L)                                # 1 mismatch over all 479
  )
  qcode <- rep(1L, 479)
  got <- malaviR:::.qc_nearest(qcode, refcode, c("exact_covering", "near_full"), top_n = 2L)

  expect_equal(got$lineage[1], "exact_covering")
  expect_equal(got$distance[1], 0)
  expect_equal(got$n_comparable[1], 400)
})

test_that("an exact match over too little of the query does NOT outrank a real neighbor", {
  ## The other half of the rule above, and the bug it fixes. Distance uses pairwise
  ## deletion, so a reference determined at only 60 of the query's 479 positions has
  ## distance 0 there and, before version 1.1.2, sorted ahead of every reference that
  ## disagrees anywhere -- which lineage_qc() then reported as `known_lineage`. 60
  ## positions is 12.5% of the query, far below the 60% floor, so the genuine neighbor
  ## must win. The thin reference is still reported, just not first.
  refcode <- rbind(
    c(rep(1L, 60), rep(0L, 419)),                      # exact over only 60 positions
    c(rep(1L, 478), 2L)                                # 1 mismatch over all 479
  )
  qcode <- rep(1L, 479)
  got <- malaviR:::.qc_nearest(qcode, refcode, c("exact_thin", "near_full"), top_n = 2L)

  expect_equal(got$lineage[1], "near_full")
  expect_equal(got$distance[1], 1)
  expect_equal(got$lineage[2], "exact_thin")           # reported, not dropped
  expect_equal(got$n_comparable[2], 60)                # and its weak overlap is visible
})

test_that("a reference sharing no determined position with the query ranks last", {
  ## The extreme case of the same bug: pairwise deletion gives a reference that overlaps
  ## the query nowhere a distance of 0, because there is nothing to disagree about. It
  ## used to sort first. It must now sort last, behind a reference that genuinely differs.
  refcode <- rbind(
    c(rep(0L, 240), rep(1L, 239)),                     # determined only where the query is not
    c(rep(1L, 235), 2L, rep(0L, 243))                  # overlaps the query, 1 mismatch
  )
  qcode <- c(rep(1L, 240), rep(0L, 239))               # determined in positions 1-240 only
  got <- malaviR:::.qc_nearest(qcode, refcode, c("no_overlap", "overlapping"), top_n = 2L)

  expect_equal(got$lineage[1], "overlapping")
  expect_equal(got$lineage[2], "no_overlap")
  expect_equal(got$n_comparable[2], 0)
})

test_that("the overlap floor scales with the query, so it still bites on a partial query", {
  ## Before 1.1.2 the floor was a fixed 300 comparable positions, which a 180 bp partial
  ## query can never reach against any reference -- so the guard stopped discriminating on
  ## exactly the queries that need it. Expressed as a fraction of the query, it works at
  ## any length: for this 180 bp query the floor is 108 positions.
  qcode <- c(rep(0L, 299), rep(1L, 180))               # determined at positions 300-479
  refcode <- rbind(
    c(rep(0L, 473), rep(1L, 6)),                       # agrees, but over only 6 positions
    c(rep(1L, 479))                                    # full length, 0 mismatches... none here
  )
  refcode[2, c(310, 350, 400)] <- 2L                   # ...make it 3 real mismatches
  got <- malaviR:::.qc_nearest(qcode, refcode, c("agrees_over_6", "true_relative"), top_n = 2L)

  expect_equal(got$lineage[1], "true_relative")
  expect_equal(got$distance[1], 3)
  expect_equal(got$n_comparable[1], 180)
  expect_equal(got$lineage[2], "agrees_over_6")
})

test_that("a partial barcode is placed into the frame and screened, not rejected", {
  ## Until 2026-08-20 any length but 479 returned invalid_sequence with every metric NA.
  ## 3,340 of MalAvi's 5,368 lineages cover only part of the window, and a primer-trimmed
  ## amplicon or a one-primer read is what submitters actually send -- so the function
  ## could not screen the majority of real submissions. Worse, a partial read of a lineage
  ## MalAvi ALREADY HOLDS came back invalid rather than as an exact match, which is the
  ## one outcome lineage_qc exists to prevent.
  aln  <- extract_alignment()
  seqs <- toupper(vapply(seq_len(nrow(aln)),
                         function(i) paste(as.character(aln[i, ]), collapse = ""),
                         character(1)))
  full <- seqs[!grepl("[^ACGT]", seqs) & nchar(seqs) == 479L][1]
  name <- rownames(aln)[which(seqs == full)[1]]

  for (part in list(substr(full, 2, 479),    # 478 bp, the haem shape
                    substr(full, 2, 477),    # 476 bp, the leuc shape
                    substr(full, 1, 300))) { # a forward-primer-only read
    qc <- lineage_qc(part)
    expect_true("placed_in_malavi_frame" %in% qc$flags)
    expect_equal(qc$call, "known_lineage")
    expect_equal(qc$summary$nearest_distance, 0)
    expect_equal(qc$summary$nearest_lineage, name)
    expect_match(qc$message, "placed at frame position")
  }
})

test_that("a sequence that cannot be placed is still invalid_sequence", {
  ## Registration refuses rather than putting a query somewhere arbitrary.
  qc <- lineage_qc(paste(rep("ACGT", 30), collapse = ""))
  expect_equal(qc$call, "invalid_sequence")
  expect_true(is.na(qc$summary$nearest_distance))
})

test_that("a query longer than the frame is not placed", {
  ## Placing it would mean discarding real bases, which is a curator's decision.
  aln  <- extract_alignment()
  full <- toupper(paste(as.character(aln[1, ]), collapse = ""))
  qc <- lineage_qc(paste0(full, "ACGTACGTAC"))
  expect_equal(qc$call, "invalid_sequence")
  expect_false("placed_in_malavi_frame" %in% qc$flags)
})

test_that("a partial query with real substitutions is not called a known lineage", {
  ## The bug this guards against, end to end on the bundled alignment. Take the
  ## second half of a complete lineage, mutate three bases, hold the source lineage
  ## out of the reference, and screen it. Before version 1.1.2 the answer was
  ## `known_lineage` at distance 0 -- against a reference sharing 6 positions with the
  ## query, because pairwise deletion scores agreement over 6 positions the same as
  ## agreement over 478. 180 of the bundled sequences carry no base at all in the
  ## first 150 positions, so this is ordinary data, not a constructed reference.
  skip_on_cran()
  aln  <- extract_alignment()
  chars <- toupper(as.character(aln))
  determined <- rowSums(matrix(chars %in% c("A", "C", "G", "T"), nrow = nrow(aln)))
  i <- which(determined == 479)[1]

  ## positions 300-479 of a complete lineage, with three bases changed
  q <- rep("N", 479)
  q[300:479] <- chars[i, 300:479]
  q[c(310, 350, 400)] <- ifelse(q[c(310, 350, 400)] == "A", "G", "A")

  qc <- lineage_qc(paste(q, collapse = ""), reference = aln[-i, ],
                   chimera_check = FALSE)

  expect_false("exact_match_to_known_lineage" %in% qc$flags)
  expect_false(identical(qc$call, "known_lineage"))
  ## the reported neighbor is now one that actually shares the query's region
  expect_gte(qc$summary$n_comparable, 108)   # the 60% floor for a 180 bp query
  expect_gt(qc$summary$nearest_distance, 0)
})

test_that("an all-N query is not called a known lineage", {
  ## The extreme case: nothing can be compared, so every reference is at distance 0.
  ## Before 1.1.2 this returned call `known_lineage`, score 0.82, nearest H_ABSUP01 --
  ## alignment row 1, which simply happened to sort first among 5,365 ties.
  skip_on_cran()
  qc <- lineage_qc(paste(rep("N", 479), collapse = ""),
                   allow_ambiguity = TRUE, chimera_check = FALSE)

  expect_true("no_comparable_reference_overlap" %in% qc$flags)
  expect_false("exact_match_to_known_lineage" %in% qc$flags)
  expect_false(identical(qc$call, "known_lineage"))
  expect_true(is.na(qc$summary$nearest_distance))
  ## and no lineage name either: every reference tied at 0 comparable positions,
  ## so the "winner" was only the first row of the alignment
  expect_true(is.na(qc$summary$nearest_lineage))
  expect_equal(qc$summary$n_comparable, 0)
})

test_that("summary reports the overlap the nearest distance was measured over", {
  skip_on_cran()
  aln <- extract_alignment()
  seq <- paste(as.character(aln[1, ]), collapse = "")
  qc  <- lineage_qc(seq, chimera_check = FALSE)

  expect_true("n_comparable" %in% names(qc$summary))
  expect_equal(qc$summary$n_comparable, qc$nearest$n_comparable[1])
  expect_gt(qc$summary$n_comparable, 0)
})

test_that("a change at an invariant site is scored as one, not as an unobserved base", {
  ## .qc_score_site() required the query base to have been observed before it
  ## could call a change an invariant-site change. At an invariant site only one
  ## base has ever been seen, so that test was always FALSE: from 1.1.0 to 1.1.1
  ## the weight-4 penalty and the N_changes_at_invariant_sites flag could not fire
  ## at all, while the mutations table reported the same base as an invariant-site
  ## change. The two disagreed about the same position.
  aln <- ape::as.DNAbin(rbind(
    r1 = strsplit("aaacccgggttt", "")[[1]],
    r2 = strsplit("aaacccgggttt", "")[[1]],
    r3 = strsplit("aaacccgggtta", "")[[1]]))
  profile <- build_malavi_site_profile(aln)
  expect_true(profile$invariant[1])              # position 1 is invariant "A"

  query <- "gaacccgggttt"                        # G at the invariant first position
  ## the internal scorer works on the upper-cased query lineage_qc() hands it
  scored <- malaviR:::.qc_score_site(strsplit(toupper(query), "")[[1]], profile)
  expect_equal(scored$site_flags[1], "invariant_site_change")

  qc <- lineage_qc(query, reference = aln, expected_length = 12,
                   chimera_check = FALSE)
  expect_equal(unname(qc$counts["n_invariant_site_changes"]), 1L)
  expect_equal(unname(qc$counts["n_bases_never_observed"]), 0L)
  expect_true(any(grepl("changes_at_invariant_sites", qc$flags)))
})
