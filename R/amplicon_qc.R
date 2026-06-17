## Abundance-aware quality-control screen for denoised amplicon sequence
## variants (ASVs). Builds on lineage_qc() per variant and adds abundance flags.
## Fixed/derived abundance settings live in .amplicon_qc_settings() (internal.R).

#' Quality-control screen for denoised amplicon sequence variants
#'
#' \strong{Experimental, use with caution.} A QC screen for
#' amplicon sequence variants (ASVs) from amplifications of the 479bp barcoding MalAvi region. My intention is that you would use this after
#' sequencing with short-read deep sequencing technology (2 x 300bp), merging the pairs, and denoising with software like \pkg{dada2},
#' vsearch, or unoise. Basically, it's a second check, because my experience is that those tools will still leave you with zillions of ASVs.
#' The way it works is that you feed the function your ASVs (the actual sequence...needs to be MalAvi aligned 479bp) and the number of final reads it has (so a two column data frame).
#' The function first evaluates the ASVs \code{\link{lineage_qc}}, then it tries to flag possible erroneous ASVs by computing some statistics related
#' to their relative abundance in the pool and how genetically similar they are to abundant variants. You can use this to, for example, remove ASVs
#' that appear in less than 5% of your reads, or rare ASVs that are very close (a bp off) from very abundant ASVs. I hope it will be useful as we start to do more deep sequencing.
#'
#' The reference is coded and its site profile built once, then reused across all variants.
#'
#' @param variants A two-column \code{data.frame}: the \strong{first} column holds
#'   the ASV sequences (each MalAvi-aligned, 479 bp) and the \strong{second}
#'   column holds their final read counts. The function errors if \code{variants}
#'   is not a data.frame, has fewer than two columns, has the columns in the wrong
#'   order (first not sequences, second not numeric counts), or contains any
#'   sequence that is not the reference length (479 bp for the bundled data) --
#'   ASVs must already be aligned to the MalAvi barcode. Only the first two
#'   columns are used.
#' @param reference Reference alignment of curated lineages: a \code{DNAbin}
#'   alignment you provide, a named character vector you provide, or \code{NULL} (default) to use the
#'   bundled MalAvi alignment.
#' @param site_profile Optional precomputed per-site base profile from
#'   \code{\link{build_malavi_site_profile}} -- a table that records, for each
#'   alignment position, which bases occur across the reference lineages and how
#'   often (the consensus base, how variable the site is, and so on). It is what
#'   lets the screen judge whether a variant's base at a site looks typical. Built
#'   automatically from \code{reference} when \code{NULL}; supply it only to avoid
#'   rebuilding it when you call the function repeatedly.
#' @param version Which bundled MalAvi release to screen against. Defaults to
#'   \code{"latest"} (the newest bundled release); pass a date string such as
#'   \code{"2026-03-23"} only if you want an older one. Consulted only when
#'   \code{reference = NULL}.
#' @param min_freq Within-pool relative frequency (read fraction) below which a
#'   variant is flagged as low-frequency (default 0.01, i.e. 1\% of the reads).
#' @param nearest_neighbor_diff How many times more abundant a near (one base off)
#'   neighbor must be before a rarer variant is flagged as a likely sequencing
#'   error derived from that abundant variant (default 10).
#' @param allow_ambiguity,chimera_check Passed through to \code{\link{lineage_qc}}
#'   for each variant.
#' @return An object of class \code{malavi_amplicon_qc}: the input
#'   \code{data.frame} with added columns, the most useful being
#'   \describe{
#'     \item{\code{relative_frequency}}{the variant's read fraction within the pool.}
#'     \item{\code{distance_to_nearest_more_abundant_neighbor}, \code{fold_less_abundant_than_nearest_neighbor}}{
#'       distance to, and fold-abundance of, the nearest more-abundant neighbor sequence.}
#'     \item{\code{lineage_call}, \code{lineage_score}}{the \code{\link{lineage_qc}}
#'       verdict and score for the variant's sequence.}
#'     \item{\code{nearest_malavi_lineage}, \code{distance_to_nearest_malavi},
#'       \code{n_nonsynonymous}}{nearest known MalAvi lineage, distance to it, and
#'       non-synonymous changes versus it.}
#'     \item{\code{lineage_flags}, \code{amplicon_flags}}{the warning
#'       tags from the lineage screen and the abundance screen.}
#'     \item{\code{amplicon_call}}{the overall per-variant verdict: \code{passes},
#'       \code{review}, \code{low_frequency_review}, \code{possible_amplicon_artifact},
#'       or \code{strong_warning}.}
#'   }
#' @seealso \code{\link{lineage_qc}}, \code{\link{lineage_screen}}
#' @examples
#' ## two known lineages plus a rare one-base error derivative of the first
#' aln <- extract_alignment()
#' s1 <- paste(as.character(aln[1, ]), collapse = "")
#' s2 <- paste(as.character(aln[2, ]), collapse = "")
#' err <- s1; substr(err, 10, 10) <- if (substr(s1, 10, 10) == "a") "c" else "a"
#' variants <- data.frame(sequence = c(s1, s2, err), count = c(10000, 4000, 5),
#'                        stringsAsFactors = FALSE)
#' aqc <- amplicon_qc(variants)
#' aqc
#' @export
amplicon_qc <- function(variants, reference = NULL, site_profile = NULL,
                        version = "latest", min_freq = 0.01,
                        nearest_neighbor_diff = 10,
                        allow_ambiguity = FALSE, chimera_check = TRUE) {
  ## ---- validate the two-column input: sequences then counts ----
  if (!is.data.frame(variants) || ncol(variants) < 2) {
    stop("`variants` must be a data.frame whose first column is the ASV ",
         "sequences and whose second column is their read counts.", call. = FALSE)
  }
  if (is.numeric(variants[[1]]) || !is.numeric(variants[[2]])) {
    stop("`variants` columns look wrong: the first column must be the ",
         "sequences (character) and the second the read counts (numeric).",
         call. = FALSE)
  }

  ## standardize to sequence/count columns (only the first two are used)
  dat <- data.frame(sequence = as.character(variants[[1]]),
                    count     = as.numeric(variants[[2]]),
                    stringsAsFactors = FALSE)
  dat$sequence <- vapply(dat$sequence, .qc_clean_seq, character(1))

  ## settings: fixed lineage weights (479 bp barcode) + the two abundance knobs
  lineage_settings  <- .lineage_qc_settings()
  amplicon_settings <- .amplicon_qc_settings(min_freq, nearest_neighbor_diff)

  ## build the coded reference + site profile once, shared across all variants
  charmat   <- .qc_char_matrix(reference, version)

  ## every ASV must already be aligned to the reference barcode, i.e. have the
  ## same length as the reference alignment (479 bp for the bundled MalAvi data)
  ref_len  <- ncol(charmat)
  seq_lens <- nchar(dat$sequence)
  if (any(seq_lens != ref_len)) {
    bad_lens <- sort(unique(seq_lens[seq_lens != ref_len]))
    stop(sum(seq_lens != ref_len), " of ", nrow(dat), " sequence(s) are not ",
         ref_len, " bp, the reference alignment length (observed length(s): ",
         paste(utils::head(bad_lens, 5), collapse = ", "),
         if (length(bad_lens) > 5) ", ..." else "",
         "). ASVs must already be aligned to the MalAvi barcode before screening; ",
         "see blast_malavi() to place an unaligned sequence.", call. = FALSE)
  }

  refcode   <- .qc_code_matrix(charmat)
  ref_names <- rownames(charmat)
  if (is.null(site_profile)) site_profile <- .qc_site_profile(charmat)
  code <- .qc_genetic_code_4()

  total_reads <- sum(dat$count, na.rm = TRUE)
  dat$relative_frequency <- dat$count / total_reads
  ## sort most abundant first so each variant's neighbors are the rows above it
  dat <- dat[order(dat$count, decreasing = TRUE), , drop = FALSE]

  ## integer-code the variants once for fast pairwise distances
  seq_codes <- lapply(strsplit(dat$sequence, "", fixed = TRUE), .qc_code_vec)

  ## nearest more-abundant variant for each variant (NA for the most abundant)
  neighbor <- .amplicon_nearest_neighbor(seq_codes, dat$count)
  dat <- cbind(dat, neighbor)

  ## biological plausibility of each variant via lineage_qc()
  lineage_results <- lapply(dat$sequence, function(q) {
    .lineage_qc_core(q, charmat, refcode, ref_names, site_profile, code,
                     allow_ambiguity, chimera_check, lineage_settings,
                     details = FALSE)
  })
  dat$lineage_call  <- vapply(lineage_results, function(x) x$call, character(1))
  dat$lineage_score <- vapply(lineage_results, function(x) x$score, numeric(1))
  dat$nearest_malavi_lineage <- vapply(lineage_results,
                                       function(x) x$nearest$lineage[1], character(1))
  dat$distance_to_nearest_malavi <- vapply(lineage_results,
                                           function(x) x$nearest$distance[1], numeric(1))
  dat$n_nonsynonymous <- vapply(lineage_results,
                                function(x) as.integer(x$summary$n_nonsynonymous), integer(1))
  dat$lineage_flags <- vapply(lineage_results,
                              function(x) paste(x$flags, collapse = ";"), character(1))

  ## abundance-aware amplicon flags + an overall call, per variant
  dat$amplicon_flags <- vapply(seq_len(nrow(dat)),
    function(i) .amplicon_flags(dat[i, ], amplicon_settings),
    character(1))
  dat$amplicon_call <- vapply(seq_len(nrow(dat)),
    function(i) .amplicon_call(dat$amplicon_flags[i], dat$lineage_call[i]),
    character(1))

  dat$pool_total_reads <- total_reads
  rownames(dat) <- NULL

  attr(dat, "min_freq") <- min_freq
  class(dat) <- c("malavi_amplicon_qc", class(dat))
  dat
}

## For abundance-sorted variants (most abundant first), find each variant's
## nearest more-abundant neighbor by Hamming distance and the fold-abundance of
## that neighbor. The first (most abundant) variant has no more-abundant
## neighbor. `seq_codes` is a list of integer-coded variants aligned with
## `counts`.
.amplicon_nearest_neighbor <- function(seq_codes, counts) {
  n <- length(seq_codes)
  rows <- lapply(seq_len(n), function(i) {
    if (i == 1L) {
      return(data.frame(nearest_more_abundant_neighbor_count = NA_real_,
                        distance_to_nearest_more_abundant_neighbor = NA_real_,
                        fold_less_abundant_than_nearest_neighbor = NA_real_,
                        stringsAsFactors = FALSE))
    }
    qi <- seq_codes[[i]]
    ## distance to each more-abundant variant, over positions known in both
    dists <- vapply(seq_len(i - 1L), function(j) {
      pj <- seq_codes[[j]]
      both <- qi > 0L & pj > 0L
      sum(both & qi != pj)
    }, numeric(1))
    best <- which.min(dists)
    data.frame(nearest_more_abundant_neighbor_count = counts[best],
               distance_to_nearest_more_abundant_neighbor = dists[best],
               fold_less_abundant_than_nearest_neighbor = counts[best] / counts[i],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## Build the abundance-aware flag string for a single variant row. Two simple
## checks: (1) the variant sits below the minimum read fraction, and (2) it is
## one base off a much more abundant variant (a likely error derivative). The
## lineage_qc() verdict for bad sequences is also surfaced as a flag.
.amplicon_flags <- function(row, th) {
  flags <- character(0)

  if (row$relative_frequency < th$min_freq) {
    flags <- c(flags, paste0("below_min_freq_", th$min_freq))
  }

  dist_nb <- row$distance_to_nearest_more_abundant_neighbor
  fold_nb <- row$fold_less_abundant_than_nearest_neighbor
  if (!is.na(dist_nb) && !is.na(fold_nb) &&
      dist_nb <= th$oneoff_distance && fold_nb >= th$nearest_neighbor_diff) {
    flags <- c(flags, "oneoff_from_much_more_abundant_variant")
  }

  if (row$lineage_call %in% c("possible_error", "strong_warning",
                              "possible_chimera", "invalid_or_strong_warning",
                              "invalid_sequence")) {
    flags <- c(flags, paste0("lineage_qc_", row$lineage_call))
  }

  paste(unique(flags), collapse = ";")
}

## Reduce a variant's amplicon flags + lineage call to a single overall call.
.amplicon_call <- function(flags, lineage_call) {
  if (grepl("lineage_qc_invalid|lineage_qc_possible_error", flags)) {
    return("strong_warning")
  }
  if (grepl("oneoff_from_much_more_abundant", flags)) {
    return("possible_amplicon_artifact")
  }
  if (grepl("below_min_freq", flags)) {
    return("low_frequency_review")
  }
  if (lineage_call %in% c("known_lineage", "plausible_new_lineage")) {
    return("passes")
  }
  "review"
}

#' @export
print.malavi_amplicon_qc <- function(x, ...) {
  cat("MalAvi amplicon QC\n")
  cat("  variants: ", nrow(x), "\n", sep = "")
  cat("  minimum relative frequency: ", attr(x, "min_freq"), "\n\n", sep = "")

  cols <- intersect(
    c("sequence", "count", "relative_frequency", "amplicon_call", "lineage_call",
      "lineage_score", "nearest_malavi_lineage", "distance_to_nearest_malavi",
      "amplicon_flags"),
    names(x)
  )
  show <- as.data.frame(x)[, cols, drop = FALSE]
  ## abbreviate the (up to 479 bp) sequence so the print stays readable; the full
  ## sequences are still in the returned object
  if ("sequence" %in% names(show)) {
    s <- show$sequence
    long <- nchar(s) > 15
    show$sequence <- ifelse(long, paste0(substr(s, 1, 12), "..."), s)
  }
  print(show)
  invisible(x)
}
