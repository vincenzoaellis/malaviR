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
  ## or 3 is almost always a short amplicon padded on the wrong end, not an
  ## aberrant sequence. Here frames 1 and 2 both carry a stop and frame 3 does
  ## not, so the diagnosis should name frame 3.
  ##   frame 1: TAA CTA AGG GCC -> *  L R A   (stop)
  ##   frame 2: AAC TAA GGG CC  ->  N *  G    (stop)
  ##   frame 3: ACT AAG GGC C   ->  T K  G    (clean)
  qc <- lineage_qc("taactaagggcc", make_ref(),
                   expected_length = 12, chimera_check = FALSE)
  expect_true("contains_stop_codon" %in% qc$flags)
  expect_true("possible_frame_shift_check_padding" %in% qc$flags)
  expect_match(qc$message, "frame 3")
  expect_match(qc$message, "padded on the wrong end")
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
