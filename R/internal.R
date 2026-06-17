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
  modal <- function(x) {
    x <- x[!is.na(x) & x != "" & x != "N/A"]
    if (!length(x)) return(NA_character_)
    names(sort(table(x), decreasing = TRUE))[1]
  }
  species <- sort(unique(sp))
  data.frame(
    species = species,
    family  = vapply(species, function(s) modal(fam[sp == s]), character(1)),
    order   = vapply(species, function(s) modal(ord[sp == s]), character(1)),
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
