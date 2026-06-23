## Internal helpers for locating and loading the bundled MalAvi releases.
## Release data ships in inst/extdata as:
##   malavi_db_<date>.rds     (tables + alignment)
##   malavi_blast_<date>.rds  (DECIPHER db + index)
## Versions are the <date> tags in the file names. These helpers are not exported.

## List the version tags available for a given file prefix, newest first.
.malavi_versions <- function(prefix = "malavi_db_") {
  dir <- system.file("extdata", package = "malaviR")
  files <- list.files(dir, pattern = paste0("^", prefix, ".*\\.rds$"))
  if (length(files) == 0) return(character(0))
  vers <- sub(paste0("^", prefix), "", sub("\\.rds$", "", files))
  sort(vers, decreasing = TRUE)
}

## Resolve a requested version ("latest" or an explicit date tag) to a file path.
.malavi_file <- function(version = "latest", prefix = "malavi_db_") {
  vers <- .malavi_versions(prefix)
  if (length(vers) == 0) {
    stop("No bundled MalAvi data found in the package. Reinstall malaviR.", call. = FALSE)
  }
  if (identical(version, "latest")) {
    version <- vers[1]
  } else if (!version %in% vers) {
    stop("Version '", version, "' is not bundled. Available: ",
         paste(vers, collapse = ", "), call. = FALSE)
  }
  system.file("extdata", paste0(prefix, version, ".rds"), package = "malaviR")
}

## Load a bundled release object (a named list).
.malavi_load <- function(version = "latest", prefix = "malavi_db_") {
  readRDS(.malavi_file(version, prefix))
}

## ---------------------------------------------------------------------------
## Haplotype grouping, shared by clean_alignment() and synonymy_report().
##
## Each sequence is integer-coded (a/c/g/t -> 1:4, gaps/Ns/ambiguities -> 0,
## i.e. "missing"). Sequences are grouped into haplotypes by one of two methods:
##   "strict"  : byte-identical across the whole alignment (gaps included).
##   "overlap" : also collapse a partial sequence into any strictly more complete
##               sequence that contains it over the partial's informative
##               (non-gap/non-N) positions -- the synonymies highlighted by
##               Tamayo-Quintero et al. (2025).
## Returns a list: lineages, group (integer id per sequence), informative_length,
## and aln_length (alignment width).
## ---------------------------------------------------------------------------
.haplotype_groups <- function(alignment, method = c("strict", "overlap")) {
  method <- match.arg(method)
  seq_mat <- as.character(alignment)
  lineages <- rownames(alignment)

  code <- matrix(0L, nrow(seq_mat), ncol(seq_mat))
  code[seq_mat == "a"] <- 1L
  code[seq_mat == "c"] <- 2L
  code[seq_mat == "g"] <- 3L
  code[seq_mat == "t"] <- 4L
  informative_length <- rowSums(code > 0L)
  names(informative_length) <- lineages

  ## strict groups: identical rows share an id (raw characters, gaps included)
  seq_key <- apply(seq_mat, 1, paste, collapse = "")
  group <- match(seq_key, unique(seq_key))

  if (method == "overlap") {
    group <- .merge_overlap(code, informative_length, group)
  }

  list(lineages = lineages, group = group,
       informative_length = informative_length, aln_length = ncol(seq_mat))
}

## Merge each strict group into a strictly more complete strict group that
## contains it over its informative positions. Operates on one representative
## row per strict group and returns a re-labelled group vector.
.merge_overlap <- function(code, informative_length, group) {
  rep_rows <- which(!duplicated(group))
  old_id   <- group[rep_rows]
  M        <- code[rep_rows, , drop = FALSE]
  infl     <- informative_length[rep_rows]
  nrep     <- nrow(M)
  max_infl <- max(infl)

  ## Blocking index to avoid the O(nrep^2) all-pairs comparison. A container
  ## must agree with the short sequence at *every* one of the short sequence's
  ## informative positions -- in particular at any single position we pick. So
  ## for each short sequence we first restrict candidates to the reps that match
  ## it at its most discriminating informative position (the (column, base) that
  ## is rarest among the reps), then do the full check only on that small set.
  ## This cannot miss a true container, so the resulting groups are identical to
  ## the naive all-pairs version, just far fewer comparisons.
  ## base_count[j, b] = number of reps carrying base b (1:4) at column j.
  base_count <- matrix(0L, ncol(M), 4L)
  for (b in 1:4) base_count[, b] <- colSums(M == b)

  merged_to <- seq_len(nrep)                 # rep index each rep merges into
  for (a in order(infl)) {                    # shortest first
    ## a zero-information sequence, or one already as complete as any other,
    ## cannot be contained in a *strictly* more complete sequence
    if (infl[a] == 0 || infl[a] == max_infl) next
    pos   <- which(M[a, ] > 0L)
    bases <- M[a, pos]
    ## pick a's most discriminating informative position (rarest base there),
    ## then keep only strictly-more-complete reps that match a at that position
    p0   <- pos[which.min(base_count[cbind(pos, bases)])]
    cand <- which(infl > infl[a] & M[, p0] == M[a, p0])
    if (length(cand) == 0) next
    eq  <- M[cand, pos, drop = FALSE] == matrix(bases, length(cand), length(pos),
                                                byrow = TRUE)
    contains <- cand[rowSums(eq) == length(pos)]
    if (length(contains) > 0) {
      ## merge into the most complete container; break ties alphabetically by
      ## lineage name so the result does not depend on alignment row order
      best <- contains[infl[contains] == max(infl[contains])]
      best <- best[order(names(infl)[best])][1]
      merged_to[a] <- best
    }
  }

  ## resolve chains to roots, then relabel groups
  root <- function(x) { while (merged_to[x] != x) x <- merged_to[x]; x }
  roots <- vapply(seq_len(nrep), root, integer(1))
  new_for_old <- old_id[roots]
  names(new_for_old) <- old_id
  unname(new_for_old[as.character(group)])
}

## ---------------------------------------------------------------------------
## Taxonomy matching helpers, used by match_taxonomy().
## ---------------------------------------------------------------------------

## Per-species MalAvi family and order from the Hosts and Sites table.
## Each host can appear on many rows; take the most frequent family/order label.
.host_family_order <- function(hosts) {
  sp  <- trimws(hosts$SPECIES_NAME)
  fam <- trimws(hosts$FAMILY_NAME)
  ord <- trimws(hosts$ORDER_NAME)
  ok  <- !is.na(sp) & sp != ""
  sp <- sp[ok]; fam <- fam[ok]; ord <- ord[ok]
  species <- sort(unique(sp))
  data.frame(
    species = species,
    family  = vapply(species, function(s) .modal(fam[sp == s]), character(1)),
    order   = vapply(species, function(s) .modal(ord[sp == s]), character(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

## Specific epithet of a binomial, and the same with a trailing Latin
## gender/declension ending removed so masculine/feminine/neuter forms collapse
## (e.g. aegyptiacus/aegyptiaca -> aegyptiac; kingi/kingii -> king).
.epithet      <- function(x) sub("^[A-Za-z]+ ", "", x)
.epithet_stem <- function(e) sub("(us|a|um|is|os|on|ii|i|e)$", "", e)

## Resolve one binomial to an eBird species name, trying exact, then IOC/
## BirdLife/Howard & Moore synonyms, then the family/order-constrained epithet
## match. Returns list(ebird, type); ebird is NA if nothing resolves.
.resolve_name <- function(name, family, order, ref) {
  if (name %in% ref$SCI_NAME) return(list(ebird = name, type = "exact"))
  syn <- c(IOC = "IOC_name", BirdLife = "Birdlife_name", HowardMoore = "H_M_name")
  for (label in names(syn)) {
    i <- match(name, ref[[syn[label]]])
    if (!is.na(i)) return(list(ebird = ref$SCI_NAME[i], type = paste0("synonym:", label)))
  }
  .epithet_reassign(name, family, order, ref)
}

## Recover a genus reassignment / gender change by matching the specific epithet
## within the host's MalAvi family (Latin family token used by clootl), falling
## back to its order when that family name is not one clootl uses. The match is
## accepted only when the epithet resolves to a single eBird species, so epithet
## collisions between unrelated birds are left unmatched.
.epithet_reassign <- function(name, family, order, ref) {
  e <- .epithet(name)
  s <- .epithet_stem(e)
  pool  <- ref[!is.na(family) & ref$latin_family == family, , drop = FALSE]
  level <- "family"
  if (nrow(pool) == 0) {
    pool  <- ref[!is.na(order) & ref$ORDER1 == order, , drop = FALSE]
    level <- "order"
  }
  if (nrow(pool) == 0) return(list(ebird = NA_character_, type = NA_character_))
  cand <- unique(pool$SCI_NAME[.epithet(pool$SCI_NAME) == e])            # exact epithet
  if (length(cand) != 1)
    cand <- unique(pool$SCI_NAME[.epithet_stem(.epithet(pool$SCI_NAME)) == s])  # gender-relaxed
  if (length(cand) == 1)
    return(list(ebird = cand, type = paste0("reassigned:", level)))
  list(ebird = NA_character_, type = NA_character_)
}

## Build the synonymy table (groups with >1 lineage) and pick which to keep.
## select = "complete" (default) keeps the most complete sequence per group
## (ties: alphabetical); select = "random" keeps one member at random.
.build_synonymies <- function(lineages, group, informative_length,
                              select = "complete", keep = NULL) {
  sizes <- table(group)
  multi <- as.integer(names(sizes))[sizes > 1]

  rows <- list()
  kept <- character(0)
  for (i in seq_along(multi)) {
    members <- lineages[group == multi[i]]
    members <- members[order(-informative_length[members], members)]
    ## members are sorted most-complete first; random picks any member instead
    chosen <- if (select == "random") sample(members, 1L) else members[1]
    user_choice <- intersect(keep, members)
    if (length(user_choice) > 1) {
      stop("More than one 'keep' lineage falls in the same haplotype group: ",
           paste(user_choice, collapse = ", "), call. = FALSE)
    }
    if (length(user_choice) == 1) chosen <- user_choice
    kept <- c(kept, chosen)
    rows[[i]] <- data.frame(
      haplotype          = i,
      lineage            = members,
      informative_length = unname(informative_length[members]),
      status             = ifelse(members == chosen, "kept", "dropped"),
      stringsAsFactors   = FALSE
    )
  }
  synonymies <- if (length(rows)) do.call(rbind, rows) else
    data.frame(haplotype = integer(0), lineage = character(0),
               informative_length = integer(0), status = character(0))
  rownames(synonymies) <- NULL
  list(synonymies = synonymies, kept = kept)
}

## ---------------------------------------------------------------------------
## Sequence-QC helpers, shared by lineage_qc() and
## build_malavi_site_profile(). These are intentionally dependency-light (no
## Biostrings/DECIPHER): everything works with base R so the QC functions stay
## part of the package core. Not exported.
##
## Throughout, a "reference" is the curated, aligned MalAvi barcode (479 bp by
## default). Sequences are handled in two parallel representations:
##   - a character matrix of upper-case bases (rows = sequences, cols = sites)
##   - an integer matrix coding a/c/g/t -> 1:4 and everything else (gaps, Ns,
##     ambiguity codes) -> 0L, i.e. "uninformative". The integer coding makes
##     the hot paths (Hamming distance, nearest-lineage search, per-site base
##     counts) fast and vectorized, the same trick used by .haplotype_groups().
## ---------------------------------------------------------------------------

## Upper-case a sequence and strip whitespace. Leaves gaps and ambiguity codes
## in place so callers can decide how to treat them.
.qc_clean_seq <- function(x) toupper(gsub("\\s+", "", x))

## Resolve the `reference` argument to an upper-case character matrix (rows =
## sequences with names, cols = alignment positions). Accepts a DNAbin matrix,
## a named character vector of equal-length aligned strings, or NULL (in which
## case the bundled MalAvi alignment for `version` is used).
.qc_char_matrix <- function(reference = NULL, version = "latest") {
  if (is.null(reference)) {
    reference <- extract_alignment(version = version)
  }
  if (inherits(reference, "DNAbin")) {
    ## as.character() on a DNAbin matrix returns a lower-case character matrix
    m <- toupper(as.character(reference))
    if (is.null(rownames(m))) {
      stop("The reference alignment must have lineage names (row names).",
           call. = FALSE)
    }
    return(m)
  }
  if (is.character(reference)) {
    reference <- vapply(reference, .qc_clean_seq, character(1))
    widths <- nchar(reference)
    if (length(unique(widths)) != 1L) {
      stop("All reference sequences must be aligned to the same length.",
           call. = FALSE)
    }
    m <- do.call(rbind, strsplit(reference, "", fixed = TRUE))
    if (is.null(rownames(m))) rownames(m) <- names(reference)
    if (is.null(rownames(m))) {
      rownames(m) <- paste0("reference_", seq_len(nrow(m)))
    }
    return(m)
  }
  stop("`reference` must be a DNAbin alignment, a named character vector, or NULL.",
       call. = FALSE)
}

## Integer-code an upper-case character matrix: A/C/G/T -> 1:4, else 0L.
## Preserves dimensions and row names.
.qc_code_matrix <- function(charmat) {
  code <- matrix(0L, nrow(charmat), ncol(charmat),
                 dimnames = dimnames(charmat))
  code[charmat == "A"] <- 1L
  code[charmat == "C"] <- 2L
  code[charmat == "G"] <- 3L
  code[charmat == "T"] <- 4L
  code
}

## Integer-code a single upper-case query (character vector of bases).
.qc_code_vec <- function(qchars) {
  qcode <- integer(length(qchars))
  qcode[qchars == "A"] <- 1L
  qcode[qchars == "C"] <- 2L
  qcode[qchars == "G"] <- 3L
  qcode[qchars == "T"] <- 4L
  qcode
}

## Nearest reference lineages to a coded query, by Hamming distance computed
## only over positions where BOTH the query and the reference carry an
## unambiguous base (gaps/Ns are skipped). Vectorized over all references at
## once: no per-reference re-parsing. Returns a data.frame ordered by ascending
## distance with columns lineage, distance, n_comparable, and index (row in
## refcode).
##
## Tie-break note: distance uses pairwise deletion -- a position where the
## reference carries an N/gap is skipped, not counted as a mismatch. That is
## deliberate, so genuinely partial reference entries (e.g. lineages whose MalAvi
## record is terminally gapped) can still match a full-length query. But it also
## means a reference that is identical to the query EXCEPT for a single ambiguous
## base ties a true, fully-overlapping exact match at distance 0: the ambiguous
## position is simply dropped from the comparison. With only `order(distance)` the
## winner among such ties was then decided by arbitrary alignment order, which
## could report the ambiguous near-twin (e.g. P_CARCAR11, which carries one N) in
## place of the genuine exact match (P_SEIAUR01) -- silently mislabeling a perfect
## ASV. We therefore break ties by DESCENDING number of comparable positions, so
## the most-complete match wins: a real exact match is never displaced by an
## N-padded one, while a uniquely-matching partial entry (the only candidate at
## its distance) still wins because it has no competitor to lose the tie to.
.qc_nearest <- function(qcode, refcode, ref_names, top_n = 5L) {
  n <- nrow(refcode)
  L <- ncol(refcode)
  qmat <- matrix(qcode, n, L, byrow = TRUE)        # broadcast the query
  both_known <- (refcode > 0L) & (qmat > 0L)       # positions comparable in both
  mism <- both_known & (refcode != qmat)           # disagreements among those
  distance <- rowSums(mism)
  n_comparable <- rowSums(both_known)              # positions actually compared
  ord <- order(distance, -n_comparable)            # ties -> most-complete match wins
  ord <- ord[seq_len(min(top_n, n))]
  data.frame(lineage = ref_names[ord], distance = distance[ord],
             n_comparable = n_comparable[ord], index = ord,
             stringsAsFactors = FALSE)
}

## NCBI genetic code 4 (mold, protozoan, and coelenterate mitochondrial code):
## TGA = W (not stop), ATA = M, and TAA/TAG remain stops. This is the code that
## fits the avian haemosporidian (apicomplexan) cytochrome b barcode: across the
## bundled alignment, frame 1 (positions 1,4,7,...) is essentially stop-free
## under this code. Hardcoded to avoid a Biostrings dependency.
.qc_genetic_code_4 <- function() {
  c(
    TTT = "F", TTC = "F", TTA = "L", TTG = "L",
    TCT = "S", TCC = "S", TCA = "S", TCG = "S",
    TAT = "Y", TAC = "Y", TAA = "*", TAG = "*",
    TGT = "C", TGC = "C", TGA = "W", TGG = "W",
    CTT = "L", CTC = "L", CTA = "L", CTG = "L",
    CCT = "P", CCC = "P", CCA = "P", CCG = "P",
    CAT = "H", CAC = "H", CAA = "Q", CAG = "Q",
    CGT = "R", CGC = "R", CGA = "R", CGG = "R",
    ATT = "I", ATC = "I", ATA = "M", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "S", AGG = "S",
    GTT = "V", GTC = "V", GTA = "V", GTG = "V",
    GCT = "A", GCC = "A", GCA = "A", GCG = "A",
    GAT = "D", GAC = "D", GAA = "E", GAG = "E",
    GGT = "G", GGC = "G", GGA = "G", GGG = "G"
  )
}

## Codon position (1, 2, or 3) of a 1-based alignment position, assuming the
## barcode is in reading frame 1 (translation starts at position 1).
.qc_codon_position <- function(position) ((position - 1L) %% 3L) + 1L

## TRUE if a base change is a transition (purine<->purine or pyrimidine<->
## pyrimidine); FALSE for a transversion. Non-A/C/G/T inputs give FALSE.
.qc_is_transition <- function(from, to) {
  paste0(from, to) %in% c("AG", "GA", "CT", "TC")
}

## Translate an upper-case base vector in frame 1 under a genetic-code table.
## Any codon containing a gap/N/ambiguity (so not in the table) becomes "X".
.qc_translate <- function(qchars, code = .qc_genetic_code_4()) {
  n_codons <- length(qchars) %/% 3L
  if (n_codons == 0L) return(character(0))
  idx <- seq_len(n_codons * 3L)
  codons <- apply(matrix(qchars[idx], nrow = 3L), 2L, paste0, collapse = "")
  aa <- unname(code[codons])
  aa[is.na(aa)] <- "X"
  aa
}

## Per-site base profile of a reference alignment, given its character matrix.
## Smoothed frequencies use a pseudocount. Returns a data.frame with one row per
## alignment position (see build_malavi_site_profile() for the column meanings).
.qc_site_profile <- function(charmat, pseudocount = 0.01) {
  codemat <- .qc_code_matrix(charmat)
  L <- ncol(codemat)
  bases <- c("A", "C", "G", "T")

  ## per-position counts of each base (columns A,C,G,T)
  counts <- vapply(1:4, function(b) colSums(codemat == b), numeric(L))
  colnames(counts) <- bases
  n_nonmissing <- rowSums(counts)

  ## smoothed frequencies (pseudocount on each of the four bases)
  freqs <- (counts + pseudocount) / (n_nonmissing + length(bases) * pseudocount)

  ## major base = most common observed base; undefined (NA) for all-gap columns
  major_idx <- max.col(counts, ties.method = "first")
  major_base <- bases[major_idx]
  major_base[n_nonmissing == 0] <- NA_character_

  n_observed <- rowSums(counts > 0)
  observed_alleles <- apply(counts > 0, 1L,
                            function(keep) paste(bases[keep], collapse = ""))

  ## Shannon entropy from the empirical (unsmoothed) frequencies, 0*log0 := 0
  emp <- counts / n_nonmissing
  ent_terms <- emp * log2(emp)
  ent_terms[!is.finite(ent_terms)] <- 0
  entropy <- -rowSums(ent_terms)

  data.frame(
    position          = seq_len(L),
    codon_position    = .qc_codon_position(seq_len(L)),
    n_seqs            = n_nonmissing,
    count_A           = counts[, "A"], count_C = counts[, "C"],
    count_G           = counts[, "G"], count_T = counts[, "T"],
    freq_A            = freqs[, "A"], freq_C = freqs[, "C"],
    freq_G            = freqs[, "G"], freq_T = freqs[, "T"],
    major_base        = major_base,
    n_observed_alleles = n_observed,
    observed_alleles  = observed_alleles,
    invariant         = n_observed == 1,
    entropy           = entropy,
    stringsAsFactors  = FALSE
  )
}

## Score a query (upper-case base vector) against a site profile: per-site
## smoothed log-probability of the query base, plus a per-site flag
## ("ok", "rare_base_at_site", "invariant_site_change",
## "base_never_observed_at_site", or "ambiguous_or_invalid").
.qc_score_site <- function(qchars, profile, rare_freq = 0.01) {
  L <- nrow(profile)
  bidx <- match(qchars, c("A", "C", "G", "T"))     # NA for gaps/Ns/ambiguities
  valid <- !is.na(bidx)

  freq_mat <- as.matrix(profile[, c("freq_A", "freq_C", "freq_G", "freq_T")])
  count_mat <- as.matrix(profile[, c("count_A", "count_C", "count_G", "count_T")])

  p <- rep(NA_real_, L)
  observed <- rep(NA, L)
  p[valid] <- freq_mat[cbind(which(valid), bidx[valid])]
  observed[valid] <- count_mat[cbind(which(valid), bidx[valid])] > 0

  flags <- rep("ok", L)
  flags[!valid] <- "ambiguous_or_invalid"
  flags[valid & !observed] <- "base_never_observed_at_site"
  invariant_change <- valid & observed & profile$invariant &
    qchars != profile$major_base
  flags[invariant_change] <- "invariant_site_change"
  rare <- valid & observed & !profile$invariant & p < rare_freq
  flags[rare] <- "rare_base_at_site"

  list(log_likelihood = sum(log(p), na.rm = TRUE),
       mean_log_probability = mean(log(p), na.rm = TRUE),
       site_flags = flags)
}

## Annotate the differences between a query and its nearest reference sequence
## (both upper-case base vectors). Only positions where both carry an
## unambiguous base are considered, so gaps in a partial reference do not appear
## as spurious mutations. Returns a per-mutation data.frame.
.qc_annotate_mutations <- function(qchars, rchars, profile,
                                   code = .qc_genetic_code_4()) {
  empty <- data.frame(
    position = integer(0), nearest_base = character(0), query_base = character(0),
    codon_position = integer(0), site_entropy = numeric(0),
    site_invariant = logical(0), query_base_observed_at_site = logical(0),
    transition = logical(0), transversion = logical(0),
    nearest_codon = character(0), query_codon = character(0),
    nearest_aa = character(0), query_aa = character(0),
    synonymous = logical(0), nonsynonymous = logical(0),
    warning = character(0), stringsAsFactors = FALSE
  )

  both <- qchars %in% c("A", "C", "G", "T") & rchars %in% c("A", "C", "G", "T")
  diff_positions <- which(both & qchars != rchars)
  if (length(diff_positions) == 0L) return(empty)

  rows <- lapply(diff_positions, function(pos) {
    ## the codon (frame 1) containing this position
    codon_start <- pos - ((pos - 1L) %% 3L)
    codon_idx <- codon_start:(codon_start + 2L)
    nearest_codon <- paste0(rchars[codon_idx], collapse = "")
    query_codon   <- paste0(qchars[codon_idx], collapse = "")
    nearest_aa <- unname(code[nearest_codon]); if (is.na(nearest_aa)) nearest_aa <- "X"
    query_aa   <- unname(code[query_codon]);   if (is.na(query_aa))   query_aa   <- "X"

    observed_alleles <- strsplit(profile$observed_alleles[pos], "")[[1]]
    query_observed <- qchars[pos] %in% observed_alleles
    is_trans <- .qc_is_transition(rchars[pos], qchars[pos])

    bits <- character(0)
    if (isTRUE(profile$invariant[pos])) bits <- c(bits, "invariant_site_change")
    if (!query_observed)                bits <- c(bits, "query_base_never_observed_at_site")
    if (query_aa == "*")                bits <- c(bits, "stop_codon")
    if (query_aa != nearest_aa)         bits <- c(bits, "nonsynonymous_change")
    if (.qc_codon_position(pos) == 2L)  bits <- c(bits, "second_codon_position_change")
    if (!is_trans)                      bits <- c(bits, "transversion")

    data.frame(
      position = pos, nearest_base = rchars[pos], query_base = qchars[pos],
      codon_position = .qc_codon_position(pos),
      site_entropy = profile$entropy[pos], site_invariant = profile$invariant[pos],
      query_base_observed_at_site = query_observed,
      transition = is_trans, transversion = !is_trans,
      nearest_codon = nearest_codon, query_codon = query_codon,
      nearest_aa = nearest_aa, query_aa = query_aa,
      synonymous = nearest_aa == query_aa, nonsynonymous = nearest_aa != query_aa,
      warning = paste(bits, collapse = ";"), stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

## Fixed internal weights and cutoffs for the lineage QC score. These are
## deliberately NOT user-facing arguments: they are the heuristic guts of the
## plausibility score, tuned once, and exposing them as function arguments only
## added clutter. The two knobs users actually want to change (expected_length
## and rare_base_frequency) are plain arguments of lineage_qc(); everything else
## lives here. .lineage_qc_settings() merges the two knobs with these weights
## into the single list the QC core consumes.
.lineage_qc_weights <- function() {
  list(
    ## penalty added to the running penalty per offending site/mutation
    invariant_site_penalty  = 4,    # a change at a never-varying site
    unobserved_base_penalty = 3,    # a base never seen at that site in MalAvi
    rare_base_penalty       = 1.5,  # a base seen but rare at that site
    nonsynonymous_penalty   = 1,    # an amino-acid-changing difference
    second_position_penalty = 1.5,  # a 2nd-codon-position difference (often nonsyn)
    transversion_penalty    = 0.75, # a transversion (rarer than a transition)

    ## Hamming-distance bins to the nearest known lineage (flag wording only)
    near_known_distance = 2,
    divergent_distance  = 5,

    ## sliding-window chimera screen
    chimera_window = 120, chimera_step = 20,
    chimera_delta_threshold = 3, chimera_min_parent_switches = 2,

    ## final-score cutoffs that map the score to a call
    pass_score = 0.85, review_score = 0.60, strong_warning_score = 0.35
  )
}

## Assemble the full settings list the QC core uses from the two user-facing
## knobs plus the fixed internal weights.
.lineage_qc_settings <- function(expected_length = 479, rare_base_frequency = 0.01) {
  c(list(expected_length = expected_length,
         rare_base_frequency = rare_base_frequency),
    .lineage_qc_weights())
}

## Most frequent non-empty label in a character vector (ties broken by sort
## order). Used to assign one parasite genus / family / order to a lineage that
## appears on many host rows.
.modal <- function(x) {
  x <- x[!is.na(x) & x != "" & x != "N/A"]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

## ---------------------------------------------------------------------------
## Helpers for the database-wide screen (lineage_studies(), lineage_screen()).
## ---------------------------------------------------------------------------

## Parse the genus-prefixed alignment sequence names into a bare lineage name
## and a parasite genus. MalAvi alignment names look like "H_ABSUP01" or
## "P_SGS1_Plasmodium_relictum": a one-letter genus prefix, the lineage name,
## then an optional species suffix. The lineage name is the token between the
## first and the next underscore; the genus comes from the prefix letter.
.lineage_parse_names <- function(x) {
  genus_code <- sub("^([A-Za-z]+)_.*$", "\\1", x)
  lineage    <- sub("^[A-Za-z]+_([^_]+).*$", "\\1", x)
  no_us <- !grepl("_", x)                 # names with no prefix at all
  lineage[no_us]    <- x[no_us]
  genus_code[no_us] <- NA_character_
  map <- c(H = "Haemoproteus", L = "Leucocytozoon", P = "Plasmodium")
  data.frame(seq_name = x, genus_code = genus_code,
             genus = unname(map[genus_code]), lineage = lineage,
             stringsAsFactors = FALSE)
}

## Count, per alignment row (lineage), the "singleton" substitutions: a base that
## the lineage ALONE carries at a well-covered site (a singleton minority base,
## differing from the site consensus). Singleton substitutions are the signature
## of sequencing error -- a real, shared variant is carried by more than one
## lineage. Each singleton substitution is classified, against the consensus codon
## (frame 1, genetic code 4) with the singleton base swapped in, as synonymous,
## non-synonymous, or stop-creating (a stop is also counted as non-synonymous).
##
## A site counts only if at least `min_cov_count` of the sequences carry an
## unambiguous base there, so consensus is well defined and sparse alignment
## columns do not masquerade as singleton substitutions. Operates on the coded
## matrix (.qc_code_matrix) and the site profile (.qc_site_profile). Returns a
## data.frame with one row per sequence, aligned with the rows of `codemat`.
.qc_singleton_substitutions <- function(codemat, site_profile, min_cov_count,
                                      code = .qc_genetic_code_4()) {
  n_seq <- nrow(codemat)
  L     <- ncol(codemat)
  bases <- c("A", "C", "G", "T")
  count_mat <- as.matrix(site_profile[, c("count_A", "count_C", "count_G", "count_T")])
  major     <- site_profile$major_base          # consensus base per site (NA if all gaps)
  major_idx <- match(major, bases)

  n_subst <- integer(n_seq); n_non <- integer(n_seq)
  n_syn   <- integer(n_seq); n_stop <- integer(n_seq)

  for (b in 1:4) {
    ## sites where base b is a singleton, the site is well covered, and b is not
    ## the consensus base (so it is a genuine singleton deviation)
    cols <- which(count_mat[, b] == 1L &
                    site_profile$n_seqs >= min_cov_count &
                    major_idx != b)
    for (j in cols) {
      i <- which(codemat[, j] == b)[1]          # the unique carrier of base b
      n_subst[i] <- n_subst[i] + 1L

      ## classify against the consensus codon with the singleton base swapped in
      cs <- j - ((j - 1L) %% 3L)                # codon start (frame 1)
      if (cs + 2L > L) next                     # incomplete trailing codon -> unclassified
      cons <- major[cs:(cs + 2L)]
      if (anyNA(cons)) next                     # consensus codon undefined -> unclassified
      mut <- cons; mut[(j - cs) + 1L] <- bases[b]
      cons_aa <- code[paste0(cons, collapse = "")]
      mut_aa  <- code[paste0(mut,  collapse = "")]
      if (is.na(cons_aa) || is.na(mut_aa)) next
      if (mut_aa == "*") n_stop[i] <- n_stop[i] + 1L
      if (mut_aa == cons_aa) n_syn[i] <- n_syn[i] + 1L else n_non[i] <- n_non[i] + 1L
    }
  }

  data.frame(n_singleton_substitutions = n_subst,
             n_singleton_nonsynonymous = n_non,
             n_singleton_synonymous    = n_syn,
             n_singleton_stop          = n_stop,
             stringsAsFactors = FALSE)
}

## Crude sliding-window chimera screen. For each window along the barcode, find
## the nearest reference lineage; count how often that nearest lineage switches
## from window to window, and compare the best single full-length parent with a
## rough two-parent mosaic distance. This is a heuristic flag for manual review,
## NOT a formal recombination test.
.qc_detect_chimera <- function(qcode, refcode, ref_names,
                               window = 120L, step = 20L, top_n = 1L) {
  L <- length(qcode)
  starts <- seq(1L, L - window + 1L, by = step)
  if (utils::tail(starts, 1L) + window - 1L < L) {
    starts <- c(starts, L - window + 1L)        # make sure the tail is covered
  }

  win <- lapply(starts, function(s) {
    cols <- s:(s + window - 1L)
    nearest <- .qc_nearest(qcode[cols], refcode[, cols, drop = FALSE],
                           ref_names, top_n = top_n)
    data.frame(window_start = s, window_end = s + window - 1L,
               nearest_lineage = nearest$lineage[1],
               nearest_distance = nearest$distance[1],
               stringsAsFactors = FALSE)
  })
  windows <- do.call(rbind, win)

  parent_switches <- sum(windows$nearest_lineage[-1] !=
                           windows$nearest_lineage[-nrow(windows)])

  full_nearest <- .qc_nearest(qcode, refcode, ref_names, top_n = 5L)
  best_single_distance <- full_nearest$distance[1]

  ## rough two-parent mosaic: sum the best per-window distances, rescaled to the
  ## full barcode length (windows overlap, hence the scaling). Intentionally
  ## approximate -- see the "future refinements" note in the QC roadmap.
  best_window_distance_sum <- sum(windows$nearest_distance)
  approx_two_parent <- best_window_distance_sum *
    (L / (nrow(windows) * window))
  chimera_delta <- best_single_distance - approx_two_parent

  list(windows = windows, best_single_lineage = full_nearest$lineage[1],
       best_single_distance = best_single_distance,
       parent_switches = parent_switches,
       approximate_two_parent_distance = approx_two_parent,
       chimera_delta = chimera_delta)
}
