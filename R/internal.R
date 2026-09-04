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

## Resolve a requested version ("latest" or an explicit date tag) to the concrete
## bundled date string, for stamping onto outputs as provenance. Returns
## NA_character_ if nothing is bundled or the tag is unknown (the loaders raise
## the actual error; this is only for the metadata stamp, so it must not itself
## stop). Kept separate from .malavi_file(), which returns a path.
.malavi_resolve_version <- function(version = "latest", prefix = "malavi_db_") {
  vers <- .malavi_versions(prefix)
  if (length(vers) == 0) return(NA_character_)
  if (identical(version, "latest")) return(vers[1])
  if (!version %in% vers) return(NA_character_)
  version
}

## Attach a compact provenance list to an object as attr(x, "malavi_meta"),
## dropping any NULL entries. Used to stamp outputs with the MalAvi version,
## method, genetic code, etc. without changing their return structure (an extra
## attribute leaves names(), print.data.frame(), and element access untouched).
## The S3 print methods and interested callers read it back with
## attr(x, "malavi_meta").
.malavi_attach_meta <- function(x, ...) {
  meta <- list(...)
  meta <- meta[!vapply(meta, is.null, logical(1))]
  attr(x, "malavi_meta") <- meta
  x
}

## Pretty one-line rendering of a malavi_meta list for print methods, e.g.
## "MalAvi 2026-03-23 | method: overlap". NULL/empty -> NULL (print nothing).
.malavi_meta_line <- function(x) {
  meta <- attr(x, "malavi_meta")
  if (is.null(meta) || length(meta) == 0) return(NULL)
  bits <- character(0)
  if (!is.null(meta$malavi_version) && !is.na(meta$malavi_version))
    bits <- c(bits, paste0("MalAvi ", meta$malavi_version))
  if (!is.null(meta$method))        bits <- c(bits, paste0("method: ", meta$method))
  if (!is.null(meta$select))        bits <- c(bits, paste0("select: ", meta$select))
  if (!is.null(meta$min_comparable)) bits <- c(bits, paste0("min_comparable: ", meta$min_comparable))
  if (!is.null(meta$clootl_version)) bits <- c(bits, paste0("clootl ", meta$clootl_version))
  if (length(bits) == 0) return(NULL)
  paste(bits, collapse = "  |  ")
}

## Tidy the whitespace of every character column of a MalAvi table: collapse any
## run of whitespace -- including the stray line breaks and tabs the source
## spreadsheets embed in free-text fields (e.g. SPECIES_NAME "Setophaga citrina\n",
## COMMENT, SITE_NAME) -- to a single space, and trim the ends. Left uncleaned,
## such embedded newlines (a) break exact joins on host names (e.g. against the
## `taxonomy` dataset) even though the species is present, and (b) split one value
## into several distinct strings, which silently inflates counts of distinct hosts.
## Only formatting is altered; a value's identity is never changed and NA is
## preserved (gsub()/trimws() return NA for NA input). Non-character columns are
## left untouched.
.clean_table_ws <- function(df) {
  char_cols <- which(vapply(df, is.character, logical(1)))
  for (j in char_cols) {
    ## PCRE rather than the POSIX class: [[:space:]] leaves a non-breaking space
    ## (U+00A0) in place, and the release has them -- 22 GENBANK_ACC values in the
    ## 2026-03-23 release begin with one. They are invisible, survive trimws(), and
    ## break an accession join exactly the way a stray newline breaks a name join.
    ## \h covers the horizontal spaces including U+00A0, \v the vertical ones.
    df[[j]] <- trimws(gsub("[\\s\\h\\v]+", " ", df[[j]], perl = TRUE))
  }
  df
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
## row per strict group and returns a re-labeled group vector.
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

## Genus (first token) and specific epithet (everything after it) of a binomial,
## plus the epithet with a trailing Latin gender/declension ending removed so
## masculine/feminine/neuter forms collapse (e.g. aegyptiacus/aegyptiaca ->
## aegyptiac; kingi/kingii -> king).
.genus        <- function(x) sub(" .*$", "", x)
.epithet      <- function(x) sub("^[A-Za-z]+ ", "", x)
.epithet_stem <- function(e) sub("(us|a|um|is|os|on|ii|i|e)$", "", e)

## Explode clootl's three alternate-authority synonym columns into one long
## lookup: one row per (synonym name, eBird species, authority).
##
## The columns do not hold one name per cell. When another authority SPLITS an
## eBird species, clootl records every one of that authority's names in the cell,
## joined by ";" -- the bundled snapshot has 115 such IOC cells, 302 BirdLife and
## 86 Howard & Moore. Comparing the cell to a name with == therefore never matched
## anything inside a joined cell, and 153 MalAvi host names sit in one. Most were
## caught by the exact step anyway, but nine had to be rescued by a hand-written
## override that the synonym step should have resolved on its own, and one --
## Phaethornis baroni, listed by BirdLife under Phaethornis longirostris -- fell
## through to the family-pool epithet step and came back as Metallura baroni.
##
## `lump` records that the name came out of a joined cell, i.e. that the authority
## recognizes it as a species while eBird folds it into a broader one. That is
## worth telling the user, because the MalAvi host concept is then narrower than
## the eBird species it maps to.
.syn_table <- function(ref) {
  sources <- c(IOC = "IOC_name", BirdLife = "Birdlife_name", HowardMoore = "H_M_name")
  parts <- lapply(names(sources), function(label) {
    values <- ref[[sources[[label]]]]
    pieces <- strsplit(values, ";", fixed = TRUE)
    data.frame(
      synonym  = trimws(unlist(pieces, use.names = FALSE)),
      ebird    = rep(ref$SCI_NAME, lengths(pieces)),
      source   = label,
      ## a name from a cell that listed more than one is a split in that authority
      lump     = rep(lengths(pieces) > 1L, lengths(pieces)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, parts)
  out[!is.na(out$synonym) & nzchar(out$synonym) & out$synonym != "NA", , drop = FALSE]
}

## Match a host name against clootl's alternate-authority synonyms (IOC, then
## BirdLife, then Howard & Moore), via the exploded table above. The same synonym
## string can sit on more than one eBird species (e.g. Howard & Moore
## "Trochalopteron cachinnans" is carried by both Montecincla cachinnans and
## M. jerdoni), and splitting the joined cells creates more such cases, not fewer
## -- BirdLife's "Saxicola torquatus" covers four eBird species. A plain match()
## would silently take whichever row comes first, so when a synonym is ambiguous
## we keep only the candidate whose own epithet agrees (allowing Latin gender)
## with the host's epithet; if that still does not single one out, we decline the
## match rather than guess. Returns list(ebird, type); ebird is NA if nothing
## resolves. `syn` is passed in by match_taxonomy() so the table is built once
## per call rather than once per host name.
.syn_resolve <- function(name, ref, syn = .syn_table(ref)) {
  s <- .epithet_stem(.epithet(name))
  for (label in c("IOC", "BirdLife", "HowardMoore")) {
    rows <- which(syn$source == label & syn$synonym == name)
    if (length(rows) == 0) next
    cand <- unique(syn$ebird[rows])
    if (length(cand) > 1) {                       # ambiguous synonym: prefer epithet agreement
      cand <- unique(cand[.epithet_stem(.epithet(cand)) == s])
      rows <- rows[syn$ebird[rows] %in% cand]
    }
    if (length(cand) == 1) {
      ## "-lump" says eBird treats this authority's species as part of a broader
      ## one, so the MalAvi host concept is narrower than the name it maps to
      suffix <- if (any(syn$lump[rows])) "-lump" else ""
      return(list(ebird = cand, type = paste0("synonym:", label, suffix)))
    }
  }
  list(ebird = NA_character_, type = NA_character_)
}

## Recover a species whose genus name is unchanged but whose epithet has shifted
## (usually Latin gender agreement, e.g. Saxicola maura -> Saxicola maurus). We
## search the whole taxonomy for the same genus, taking an exact-epithet hit
## first and a gender-relaxed one otherwise, and accept it only when it resolves
## to a single eBird species. Keeping the genus fixed makes this the strongest
## identity signal available, so it is tried before the family/order-constrained
## search, which can otherwise be misled by a wrong MalAvi family label into a
## same-epithet match in an unrelated genus. Returns list(ebird, type).
.same_genus_reassign <- function(name, ref) {
  g    <- .genus(name)
  pool <- ref[.genus(ref$SCI_NAME) == g, , drop = FALSE]
  if (nrow(pool) == 0) return(list(ebird = NA_character_, type = NA_character_))
  cand <- unique(pool$SCI_NAME[.epithet(pool$SCI_NAME) == .epithet(name)])        # exact epithet
  if (length(cand) != 1)
    cand <- unique(pool$SCI_NAME[.epithet_stem(.epithet(pool$SCI_NAME)) ==
                                   .epithet_stem(.epithet(name))])                # gender-relaxed
  if (length(cand) == 1)
    return(list(ebird = cand, type = "reassigned:genus"))
  list(ebird = NA_character_, type = NA_character_)
}

## Resolve one binomial to an eBird species name, trying exact, then the
## (ambiguity-aware) IOC/BirdLife/Howard & Moore synonyms, then a same-genus
## epithet shift, then the family/order-constrained epithet match. Returns
## list(ebird, type); ebird is NA if nothing resolves.
.resolve_name <- function(name, family, order, ref, syn_table = NULL) {
  if (name %in% ref$SCI_NAME) return(list(ebird = name, type = "exact"))
  syn <- if (is.null(syn_table)) .syn_resolve(name, ref) else .syn_resolve(name, ref, syn_table)
  if (!is.na(syn$ebird)) return(syn)
  sg <- .same_genus_reassign(name, ref)
  if (!is.na(sg$ebird)) return(sg)
  .epithet_reassign(name, family, order, ref)
}

## Recover a genus reassignment / gender change by matching the specific epithet
## within the host's MalAvi family (Latin family token used by clootl), falling
## back to its order when that family name is not one clootl uses. The match is
## accepted only when the epithet resolves to a single eBird species, so epithet
## collisions between unrelated birds are left unmatched.
##
## This is the weakest step in the resolver, and the reason is that the pool is
## built from MalAvi's FAMILY_NAME -- the least maintained field in the release.
## When MalAvi files a genus under an old or wrong family, the pool becomes the
## CURRENT clootl membership of that family name, the true species is not in it
## (clootl has moved it elsewhere), and any lone same-epithet bird in the pool
## wins. That produced five false species in the shipped crosswalk, among them
## Tiaris obscura -> Akialoa obscura, an extinct Hawaiian honeycreeper standing in
## for a Peruvian grassquit.
##
## So the epithet alone is not enough to move a host into a different genus. Where
## clootl still uses the MalAvi genus, it already says which families that genus
## belongs to, and a candidate outside them is rejected. Note this is deliberately
## about the genus's home in clootl, not about MalAvi's family label -- the label
## is what is untrustworthy here. When clootl has retired the MalAvi genus
## altogether (Hemispingus, say) there is nothing to check against and the step is
## as blind as before; those cases need a row in data-raw/manual_taxonomy.csv.
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
  if (length(cand) != 1) return(list(ebird = NA_character_, type = NA_character_))

  ## the genus guard described above: only applies when the match moves the host
  ## into a different genus AND clootl still recognizes the MalAvi genus
  if (.genus(cand) != .genus(name)) {
    home_families <- unique(ref$latin_family[.genus(ref$SCI_NAME) == .genus(name)])
    cand_family   <- unique(ref$latin_family[ref$SCI_NAME == cand])
    if (length(home_families) > 0 && !any(cand_family %in% home_families))
      return(list(ebird = NA_character_, type = NA_character_))
  }
  list(ebird = cand, type = paste0("reassigned:", level))
}

## Flag the crosswalk rows that rest on the weakest evidence, so a maintainer can
## eyeball them after every rebuild (called from build_taxonomy.R). These are the
## matches with the least independent support, and the only place the Oriolus
## brachyrhynchus -> Corvus brachyrhynchos style collision can hide:
##   - "weak_reassignment": a family/order-pool epithet match that *changed the
##     genus* and only agreed after Latin gender relaxation (not an exact
##     epithet). A wrong MalAvi family label can drive a coincidental stem into
##     an unrelated genus here. (Same-genus shifts and synonym matches carry
##     stronger support and are not flagged.)
##   - "retired_genus": a family/order-pool epithet match that changed the genus,
##     where clootl no longer uses the MalAvi genus at all. This is the exact blind
##     spot of the genus guard in .epithet_reassign(): with no clootl entry for the
##     MalAvi genus there is no home family to check the candidate against, so the
##     match rests on the epithet and MalAvi's family label alone. Most of these are
##     the standard transfers (Dendroica -> Setophaga, Megalaima -> Psilopogon), but
##     Hemispingus frontalis -> Crithagra frontalis, an African finch standing in for
##     a Peruvian tanager, hid here through several releases. Flagged whether or not
##     the epithet matched exactly, because an exact epithet is no evidence at all
##     when the pool is the wrong family.
##   - "legacy": resolved only through the decades-old original malaviR hand key,
##     which predates the current eBird taxonomy and is worth re-checking.
## Manual overrides are excluded -- a human already decided those. NOTE this is a
## review aid, not a correctness test: the great majority of flagged rows are
## legitimate (e.g. Carduelis -> Spinus gender agreement); the point is that the
## genuinely wrong ones, when they occur, will be somewhere in this short list.
## Returns a data frame of the rows to review; empty if none.
.audit_taxonomy <- function(key) {
  genus_changed <- .genus(key$malavi_species) != .genus(key$ebird_species)
  epithet_exact <- .epithet(key$malavi_species) == .epithet(key$ebird_species)
  ## genera clootl still recognizes; a MalAvi genus absent from this is one the
  ## genus guard could not check
  known_genus <- .genus(key$malavi_species) %in% unique(.genus(clootl_ref$SCI_NAME))

  pooled <- key$match_type %in% c("reassigned:family", "reassigned:order") &
    genus_changed & !is.na(key$ebird_species)
  weak    <- pooled & !epithet_exact
  retired <- pooled & !known_genus & !weak     # weak already lists the relaxed ones
  legacy  <- key$match_type == "legacy"

  flag   <- weak | retired | legacy
  reason <- ifelse(legacy[flag], "legacy",
                   ifelse(retired[flag], "retired_genus", "weak_reassignment"))
  out <- data.frame(
    malavi_species = key$malavi_species[flag],
    ebird_species  = key$ebird_species[flag],
    match_type     = key$match_type[flag],
    reason         = reason,
    stringsAsFactors = FALSE
  )
  out[order(out$reason, out$malavi_species), , drop = FALSE]
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

## Minimum overlap, as a fraction of the query's determined positions, before a
## reference's comparison with the query is trusted.
##
## The gate is RELATIVE to the query, not a fixed number of positions. Until
## version 1.1.2 it was a fixed 300, which made it a no-op for exactly the queries
## that need it: a 180 bp partial query can never reach 300 comparable positions
## against ANY reference, so every reference counted as thin and the guard stopped
## discriminating. Expressed as a fraction, the same rule applies at every query
## length. For a full-length 479 bp query the floor is 288 positions -- essentially
## the fixed 300 used before -- and 81 of the 5,365 sequences in the bundled
## alignment fall below it (1.5%); the median reference overlap is 478.
.QC_MIN_COMPARABLE_FRACTION <- 0.6

## Positions a reference must share with the query before its comparison counts.
## `q_known` is the number of positions the QUERY has an unambiguous base at. The
## floor is never 0: a reference that shares nothing with the query must always
## fail, including when the query itself is entirely unknown (an all-N query, where
## `q_known` is 0 and every reference would otherwise pass a floor of 0).
.qc_min_comparable <- function(q_known) {
  max(1L, as.integer(ceiling(.QC_MIN_COMPARABLE_FRACTION * q_known)))
}

## Nearest reference lineages to a coded query, by Hamming distance computed
## only over positions where BOTH the query and the reference carry an
## unambiguous base (gaps/Ns are skipped). Vectorized over all references at
## once: no per-reference re-parsing. Returns a data.frame ordered as described
## below with columns lineage, distance, n_comparable, and index (row in refcode).
##
## ORDERING: order(exact, thin, rate, -n_comparable), i.e.
##   1. exact matches backed by enough overlap, then
##   2. references that share enough of the query, by ascending mismatch RATE, then
##   3. references that share too little of the query, again by rate,
##   with ties at every level broken by descending overlap.
##
## Three separate problems shaped this, and each of the keys answers one of them.
##
## RATE, not count (`rate`). Ordering on the raw mismatch count lets a reference that
## overlaps the query in few positions win simply by having less to disagree over.
## Measured on a real submission (2026-08-20): a candidate lineage was reported nearest
## to a reference with 23 mismatches over only 133 comparable positions -- 82.7%
## identity, and that reference is the least-covered sequence in the bundled alignment.
## Its true nearest relative had 38 mismatches over 477 positions: 92.0% identity, and
## the clade the lineage actually belongs to. Reporting the first as "nearest" pointed a
## curator at the wrong genus.
##
## ENOUGH OVERLAP BEFORE A MATCH IS EVIDENCE (`exact`, `thin`). Distance uses pairwise
## deletion, so a position where either sequence is N or a gap is skipped rather than
## counted as a mismatch. That is deliberate -- genuinely partial reference entries have
## to be able to match a full-length query -- but it means a reference that shares NO
## determined position with the query also has distance 0. Before 1.1.2 the exact bucket
## was `distance == 0` alone, so such a reference sorted ahead of every reference that
## actually disagrees somewhere, and lineage_qc() called the query a known lineage. Two
## confirmed cases: a 180 bp partial query with three real substitutions was reported
## exactly matching a reference sharing 6 positions with it, and an all-N query matched
## alignment row 1 over 0 positions. Since partial queries began to be placed and screened
## (2026-08-20) this is ordinary use, not a corner case: 3,338 of the 5,365 bundled
## sequences are partial. So a distance-0 agreement only enters the exact bucket when it
## rests on at least `.qc_min_comparable(q_known)` shared positions, and references below
## that floor sort last whatever their distance. They are ordered, not dropped, and
## `n_comparable` is returned so the caller can see how much a match rests on.
##
## TIES AMONG GENUINE EXACT MATCHES (`-n_comparable`). A reference identical to the query
## except for a single ambiguous base ties a true, fully-overlapping exact match at
## distance 0, because the ambiguous position is dropped from the comparison. With only
## `order(distance)` the winner was then decided by arbitrary alignment order, which could
## report the ambiguous near-twin (e.g. P_CARCAR11, which carries one N) in place of the
## genuine exact match (P_SEIAUR01) -- silently mislabeling a perfect ASV. Breaking ties by
## DESCENDING overlap makes the most-complete match win, while a uniquely-matching partial
## entry still wins because it has no competitor to lose the tie to.
.qc_nearest <- function(qcode, refcode, ref_names, top_n = 5L) {
  n <- nrow(refcode)
  L <- ncol(refcode)
  qmat <- matrix(qcode, n, L, byrow = TRUE)        # broadcast the query
  both_known <- (refcode > 0L) & (qmat > 0L)       # positions comparable in both
  mism <- both_known & (refcode != qmat)           # disagreements among those
  distance <- rowSums(mism)
  n_comparable <- rowSums(both_known)              # positions actually compared

  ## how much of the query a reference must cover before its comparison is trusted
  min_comparable <- .qc_min_comparable(sum(qcode > 0L))
  enough <- n_comparable >= min_comparable

  exact <- as.integer(!(distance == 0L & enough))  # 0 sorts first, so 0 == "exact"
  thin  <- as.integer(!enough)
  rate  <- ifelse(n_comparable > 0L, distance / pmax(n_comparable, 1L), Inf)
  ord <- order(exact, thin, rate, -n_comparable)
  ord <- ord[seq_len(min(top_n, n))]
  data.frame(lineage = ref_names[ord], distance = distance[ord],
             n_comparable = n_comparable[ord], index = ord,
             stringsAsFactors = FALSE)
}

## NCBI genetic code 4 (mold, protozoan, and coelenterate mitochondrial code).
## Code 4 is identical to the standard code except that TGA is tryptophan
## rather than a stop; TAA and TAG remain the only stops, and ATA/AGA/AGG keep
## their standard meanings (I/R/R). Do not copy the invertebrate mitochondrial
## code (table 5) here: it also sets ATA = M and AGA/AGG = S, which is wrong for
## haemosporidians. This is the code that fits the avian haemosporidian
## (apicomplexan) cytochrome b barcode: across the bundled alignment, frame 1
## (positions 1,4,7,...) is essentially stop-free under this code. Hardcoded to
## avoid a Biostrings dependency; test-lineage_qc.R checks it against the
## standard code plus the single TGA difference.
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
    ATT = "I", ATC = "I", ATA = "I", ATG = "M",
    ACT = "T", ACC = "T", ACA = "T", ACG = "T",
    AAT = "N", AAC = "N", AAA = "K", AAG = "K",
    AGT = "S", AGC = "S", AGA = "R", AGG = "R",
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

## Number of stop codons obtained by translating an upper-case base vector in
## each of the three forward reading frames (frame f starts at position f).
## Used only to diagnose a query that is the right length but appears to have
## been padded on the wrong end, which shifts it out of frame 1. Returns an
## integer vector of length 3, named "1", "2", "3".
.qc_frame_stop_counts <- function(qchars, code = .qc_genetic_code_4()) {
  counts <- vapply(1:3, function(frame) {
    if (length(qchars) < frame + 2L) return(0L)
    ## drop the leading frame - 1 bases, then translate what remains in frame 1
    shifted <- qchars[seq.int(frame, length(qchars))]
    sum(.qc_translate(shifted, code) == "*")
  }, integer(1))
  names(counts) <- as.character(1:3)
  counts
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
  ## Deliberately NOT conditioned on `observed`. At an invariant site only one
  ## base has been seen, so any query base that differs from the major base is by
  ## definition unobserved there; requiring `observed` made this test always
  ## FALSE, and the weight-4 invariant-site penalty and its flag were unreachable
  ## from version 1.1.0 to 1.1.1 while the mutations table reported the change
  ## anyway. The assignment order is what separates the two categories: an
  ## invariant-site change is the more specific and the more suspicious of the
  ## two statements, so it overwrites the never-observed flag set just above.
  invariant_change <- valid & profile$invariant & qchars != profile$major_base
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

## Fixed internal cutoffs for the lineage QC screen. These are deliberately NOT
## user-facing arguments: exposing them only added clutter. The two knobs users
## actually want to change (expected_length and rare_base_frequency) are plain
## arguments of lineage_qc(); everything else lives here.
## .lineage_qc_settings() merges the two knobs with these into the single list
## the QC core consumes.
##
## Version 1.2.0 removed the weighted penalty score this list used to carry
## (invariant_site_penalty 4, unobserved_base_penalty 3, rare_base_penalty 1.5,
## nonsynonymous_penalty 1, second_position_penalty 1.5, transversion_penalty
## 0.75, mapped through exp(-penalty/10) to pass/review/strong_warning cutoffs of
## 0.85/0.60/0.35). Leave-one-out on 60 curated bundled lineages put 26% of them
## in strong_warning or possible_error, and every lineage more than 20 bp from its
## nearest neighbor came out possible_error -- including L_PIRIE01, a described
## Leucocytozoon species, at a score of 0. No term was normalized by distance, so
## the composite was in effect a divergence measure wearing a plausibility label,
## and divergence is the one property a genuinely new lineage has. The individual
## counts it was built from are kept and reported; they are checkable facts.
.lineage_qc_weights <- function() {
  list(
    ## Hamming-distance bins to the nearest known lineage (flag wording only)
    near_known_distance = 2,
    divergent_distance  = 5,

    ## sliding-window chimera screen
    chimera_window = 120, chimera_step = 20,
    ## Raised from 3 to 8 on 2026-08-20, and joined by a parent-distance guard. Calibrated
    ## on 200 real lineages (leave-one-out) and 200 synthetic chimeras: 91.5% of true
    ## chimeras caught at 2.0% false positives, against 23.2% false positives before.
    ## See the note above .qc_detect_chimera. chimera_min_parent_switches is retained and
    ## reported but no longer takes part in the call: 98.4% of ordinary lineages met it.
    chimera_delta_threshold = 8, chimera_min_parent_switches = 2,
    chimera_min_parent_distance = 5, chimera_min_segment = 60
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

## Two-parent breakpoint chimera screen.
##
## For every breakpoint b, the best two-parent explanation of the query is
##   min_A mismatches(query[1..b], A) + min_B mismatches(query[b+1..L], B),
## computed for all references at once from row-wise cumulative sums. `chimera_delta` is
## how much better that is than the best single parent: the evidence that no one lineage
## explains the query but two do.
##
## WHAT THIS REPLACED, AND WHY (2026-08-20). The previous screen counted how often the
## nearest lineage changed between overlapping 120 bp windows, and compared the best single
## parent against a rescaled sum of per-window distances. Measured by leave-one-out over the
## release, it called **23.2% of ordinary lineages** possible chimeras. Two reasons:
##
##   * `parent_switches >= 2` was met by 98.4% of ordinary lineages. With a median of 1 bp
##     to the nearest relative, which lineage wins any given window is close to arbitrary,
##     so the term excluded almost nothing.
##   * A 120 bp window nearly always finds some lineage matching it exactly, so the
##     approximate two-parent distance sat near zero and `chimera_delta` collapsed to the
##     distance to the nearest single lineage. The call was a **divergence** measure wearing
##     a chimera label -- and divergence is the one property a genuinely new lineage has.
##
## Calibrated against labeled data rather than intuition (`data-raw/chimera_v2_eval.R`):
## 200 real lineages leave-one-out as negatives, 200 synthetic chimeras spliced from real
## lineages with both parents left in the reference as positives. Negatives have median
## delta 1 (95th percentile 4); real chimeras median 18. At the shipped threshold of 8:
## **91.5% of true chimeras caught, 2.0% false positives** -- against 23.2% before.
##
## Three guards the old screen lacked, each removing a way to manufacture evidence:
##   * the two parents must be different lineages;
##   * each segment must be at least `min_segment` long, so a breakpoint a few bases from an
##     end cannot be built out of trimming noise;
##   * the parents must differ from EACH OTHER by `min_parent_distance` over the positions
##     the query defines -- alternating between two lineages 1 bp apart is not evidence;
##   * references must cover most of the query, so a thinly covered one cannot win by having
##     little to disagree over (the same failure fixed in .qc_nearest).
##
## Still a heuristic flag for manual review, NOT a formal recombination test.
.qc_detect_chimera <- function(qcode, refcode, ref_names,
                               window = 120L, step = 20L, top_n = 1L,
                               min_segment = 60L, min_cover = 0.9) {
  L <- length(qcode)
  n <- nrow(refcode)

  ## per-window nearest lineage, kept for the detail output and for parent_switches.
  ## It no longer takes part in the call -- see the note above.
  starts <- seq(1L, L - window + 1L, by = step)
  if (utils::tail(starts, 1L) + window - 1L < L) starts <- c(starts, L - window + 1L)
  win <- lapply(starts, function(s) {
    cols <- s:(s + window - 1L)
    nearest <- .qc_nearest(qcode[cols], refcode[, cols, drop = FALSE], ref_names,
                           top_n = top_n)
    data.frame(window_start = s, window_end = s + window - 1L,
               nearest_lineage = nearest$lineage[1],
               nearest_distance = nearest$distance[1], stringsAsFactors = FALSE)
  })
  windows <- do.call(rbind, win)
  parent_switches <- sum(windows$nearest_lineage[-1] !=
                           windows$nearest_lineage[-nrow(windows)])

  qmat  <- matrix(qcode, n, L, byrow = TRUE)
  known <- (refcode > 0L) & (qmat > 0L)
  mism  <- known & (refcode != qmat)

  q_known <- sum(qcode > 0L)
  keep <- rowSums(known) >= min_cover * q_known
  empty <- list(windows = windows, best_single_lineage = NA_character_,
                best_single_distance = NA_integer_, parent_switches = parent_switches,
                approximate_two_parent_distance = NA_real_, chimera_delta = 0,
                breakpoint = NA_integer_, parent_a = NA_character_,
                parent_b = NA_character_, parent_distance = NA_integer_)
  if (sum(keep) < 2L || L < 2L * min_segment) return(empty)

  mism <- mism[keep, , drop = FALSE]
  refc <- refcode[keep, , drop = FALSE]
  rn   <- ref_names[keep]

  Pm  <- t(apply(mism, 1L, cumsum))          # mismatches in 1..b, per reference
  tot <- Pm[, L]
  d1  <- min(tot)
  a1  <- rn[which.min(tot)]

  bs   <- seq.int(min_segment, L - min_segment)
  pre  <- Pm[, bs, drop = FALSE]
  post <- tot - pre
  i_pre  <- apply(pre,  2L, which.min)
  i_post <- apply(post, 2L, which.min)
  d2 <- pre[cbind(i_pre, seq_along(bs))] + post[cbind(i_post, seq_along(bs))]
  d2[i_pre == i_post] <- Inf                 # one parent is not a mosaic
  if (!any(is.finite(d2))) {
    empty$best_single_lineage <- a1; empty$best_single_distance <- d1
    return(empty)
  }

  k <- which.min(d2)
  A <- i_pre[k]; B <- i_post[k]
  both <- refc[A, ] > 0L & refc[B, ] > 0L & qcode > 0L
  parent_distance <- sum(refc[A, both] != refc[B, both])

  list(windows = windows,
       best_single_lineage = a1, best_single_distance = d1,
       parent_switches = parent_switches,
       approximate_two_parent_distance = d2[k],
       chimera_delta = d1 - d2[k],
       breakpoint = bs[k], parent_a = rn[A], parent_b = rn[B],
       parent_distance = parent_distance)
}


## Consensus code per alignment column: the most frequent unambiguous base, 0 if a
## column has none. Used only for registration, where comparing against one consensus
## is both faster and steadier than comparing against every reference in turn.
.qc_consensus <- function(refcode) {
  apply(refcode, 2L, function(col) {
    col <- col[col > 0L]
    if (!length(col)) return(0L)
    tab <- tabulate(col, nbins = 4L)
    as.integer(which.max(tab))
  })
}

## Place a query that is SHORTER than the reference frame into that frame.
##
## Returns list(offset, rate, n_comparable) for the best placement, or NULL when nothing
## lands convincingly. `offset` is 0-based: the query's first base sits at frame position
## offset + 1.
##
## Why this exists: 3,340 of MalAvi's 5,368 lineages cover only part of the 479 bp barcode,
## so a partial query is the ordinary case rather than an error. Before this, lineage_qc()
## returned `invalid_sequence` for any length but 479 and skipped every sequence metric --
## including the reading-frame diagnosis that exists to explain exactly that situation,
## which sat below the gate and was unreachable. A primer-trimmed amplicon is 478 or 476 bp
## and a one-primer read is shorter still; none of them could be screened at all.
##
## Registration slides the query and scores it against the reference consensus, which is
## what the malavi_rebuild screen does. It refuses rather than guessing: a query that does
## not land clearly better than random is not silently placed somewhere arbitrary.
.qc_register <- function(qcode, refcode, max_rate = 0.35, min_comparable = 30L) {
  L <- ncol(refcode)
  qlen <- length(qcode)
  if (qlen >= L) return(NULL)

  consensus <- .qc_consensus(refcode)
  best <- NULL
  for (off in 0:(L - qlen)) {
    idx <- seq.int(off + 1L, off + qlen)
    ref <- consensus[idx]
    known <- ref > 0L & qcode > 0L
    n_known <- sum(known)
    if (n_known < min_comparable) next
    rate <- sum(ref[known] != qcode[known]) / n_known
    if (is.null(best) || rate < best$rate) {
      best <- list(offset = off, rate = rate, n_comparable = n_known)
    }
  }
  if (is.null(best) || best$rate > max_rate) return(NULL)
  best
}
