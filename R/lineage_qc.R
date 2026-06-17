## Quality-control screen for a single candidate MalAvi cytochrome b lineage.
## See internal.R for the shared, dependency-light QC helpers and for
## the fixed scoring weights (.lineage_qc_weights / .lineage_qc_settings).

#' Build a per-site base profile of a MalAvi alignment
#'
#' \strong{Experimental.} A helper used by \code{\link{lineage_qc}} and
#' \code{\link{amplicon_qc}}. In an aligned set of sequences, every column (a
#' "site", i.e. one position along the 479 bp barcode) lines up the same base
#' across all lineages. This function walks the alignment column by column and
#' records, for each site, which DNA bases (A/C/G/T) occur there and how often --
#' the consensus (most common) base, how variable the site is, and so on. The QC
#' functions use that summary to judge whether a query's base at a given site
#' looks typical of the lineages currently in MalAvi (e.g. a base never seen at
#' that site, or a change at a site that never varies, is suspicious).
#'
#' Building the profile from the bundled alignment is very fast (seconds). But, if you
#' are screening many sequences, consider building it once and passing it to \code{lineage_qc} or
#' \code{amplicon_qc} through their \code{site_profile} argument.
#'
#' @param reference A reference alignment: a \code{DNAbin} alignment (e.g. from
#'   \code{\link{extract_alignment}}) or a named character vector of equal-length
#'   aligned sequences. Leave it \code{NULL} (the default) to profile the bundled
#'   MalAvi alignment, which is the usual case; only supply \code{reference} if
#'   you want to profile your own alignment instead.
#' @param version Which bundled MalAvi release to profile. Only consulted when
#'   \code{reference = NULL} (it selects which bundled alignment to load, and is
#'   ignored if you pass your own \code{reference}). A date string such as
#'   \code{"2026-03-23"}, or \code{"latest"} (default) for the newest bundled
#'   release.
#' @param pseudocount Smoothing pseudocount added to each base count when
#'   computing per-site frequencies (default 0.01). Without it, a base that never
#'   occurs at a site would have frequency exactly 0, and the QC score (which
#'   takes logs of these frequencies) would hit \code{log(0) = -Inf}. Adding a
#'   small amount to every base lifts those zeros just off the floor. A small
#'   value (0.01) barely perturbs the observed frequencies, so genuinely
#'   unobserved bases stay very rare (and thus stay flag-worthy); a larger value
#'   would make never-seen bases look more ordinary.
#' @return A \code{data.frame} with one row per alignment position (one row per
#'   column of the alignment):
#'   \describe{
#'     \item{\code{position}}{the site's index along the barcode, i.e. its 1-based
#'       column number in the alignment (1 = first aligned base).}
#'     \item{\code{codon_position}}{where the site falls within its codon under
#'       reading frame 1: 1, 2, or 3 (position 1 -> 1, 2 -> 2, 3 -> 3, 4 -> 1, ...).
#'       Third-position sites tolerate synonymous change, so this flags how
#'       constrained the site is expected to be.}
#'     \item{\code{n_seqs}}{number of sequences with a called base (A/C/G/T)
#'       at the site (i.e. not a gap or ambiguity).}
#'     \item{\code{count_A}, \code{count_C}, \code{count_G}, \code{count_T}}{number
#'       of sequences carrying each base at the site.}
#'     \item{\code{freq_A}, \code{freq_C}, \code{freq_G}, \code{freq_T}}{the
#'       smoothed \emph{proportion} of sequences carrying each base at the site,
#'       computed as \code{(count + pseudocount) / (n_seqs + 4 * pseudocount)};
#'       the four base frequencies at a site therefore sum to 1.}
#'     \item{\code{major_base}}{the most common base at the site (the consensus).}
#'     \item{\code{n_observed_alleles}}{how many of the four bases are seen at all.}
#'     \item{\code{observed_alleles}}{those bases as a string, e.g. \code{"AG"}.}
#'     \item{\code{invariant}}{\code{TRUE} if only one base is ever seen at the site.}
#'     \item{\code{entropy}}{Shannon entropy of the observed bases; 0 at an
#'       invariant site, higher at variable ones.}
#'   }
#' @seealso \code{\link{lineage_qc}}, \code{\link{amplicon_qc}}
#' @examples
#' profile <- build_malavi_site_profile()
#' head(profile)
#' @export
build_malavi_site_profile <- function(reference = NULL, version = "latest",
                                      pseudocount = 0.01) {
  charmat <- .qc_char_matrix(reference, version)
  .qc_site_profile(charmat, pseudocount = pseudocount)
}

#' Quality-control screen for a candidate MalAvi lineage
#'
#' \strong{Experimental, use with caution.} A lightweight plausibility screen for
#' a single, already-aligned cytochrome \emph{b} barcode (479 bp by default). It
#' does \strong{not} decide whether a lineage is real; it flags features that are
#' surprising relative to curated MalAvi diversity so you can review them by hand.
#' Read warnings as "worth a look", not "this is definitely an error".
#'
#' The query must already be aligned to the MalAvi barcode (same length and
#' reading frame as the reference). For an unaligned sequence, find its closest
#' lineages first with \code{\link{blast_malavi}}. The screen rolls a handful of
#' checks into one \code{score} in \code{[0, 1]} (1 = typical of known MalAvi
#' diversity, 0 = highly suspicious or invalid). The score is a transparent,
#' rule-based heuristic, \strong{not} a calibrated probability.
#'
#' @section Flags:
#' \code{flags} is a character vector of short tags you can filter on. The ones
#' you may see:
#' \describe{
#'   \item{\code{wrong_length_*}}{the query is not the expected barcode length.}
#'   \item{\code{contains_gaps}, \code{contains_N}, \code{contains_ambiguity_codes},
#'     \code{invalid_or_disallowed_characters}}{the query is not clean A/C/G/T.}
#'   \item{\code{contains_stop_codon}}{translation (frame 1, genetic code 4)
#'     contains a stop codon -- a strong sign of an error or wrong frame.}
#'   \item{\code{exact_match_to_known_lineage}, \code{near_known_lineage},
#'     \code{moderately_divergent_from_known_lineages},
#'     \code{highly_divergent_from_known_lineages}}{how far the query sits from the
#'     nearest known lineage.}
#'   \item{\code{N_changes_at_invariant_sites}}{changes at sites that never vary
#'     in MalAvi (very suspicious).}
#'   \item{\code{N_bases_never_observed_at_their_sites}}{bases never seen at that
#'     site in MalAvi.}
#'   \item{\code{N_rare_bases_at_sites}}{bases that are real but rare at that site.}
#'   \item{\code{N_nonsynonymous_changes_vs_nearest_lineage},
#'     \code{N_second_codon_position_changes_vs_nearest_lineage},
#'     \code{N_transversions_vs_nearest_lineage}}{unusual mutation types relative
#'     to the nearest lineage.}
#'   \item{\code{possible_chimera_or_mixed_template_pattern}}{the sliding-window
#'     screen suggests a mosaic of two parents.}
#' }
#'
#' @param query A single candidate barcode as a character string. Whitespace is
#'   stripped and the sequence is upper-cased.
#' @param reference Reference alignment of curated lineages: a \code{DNAbin}
#'   alignment, a named character vector, or \code{NULL} (default) to use the
#'   bundled MalAvi alignment for \code{version}.
#' @param site_profile Optional precomputed profile from
#'   \code{\link{build_malavi_site_profile}}; built from \code{reference} when
#'   \code{NULL}. Supply it to avoid rebuilding when screening many sequences.
#' @param version MalAvi release to use when \code{reference} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @param expected_length Expected barcode length in bp (479 for MalAvi). A query
#'   of a different length is reported as \code{invalid_sequence}.
#' @param rare_base_frequency Smoothed frequency below which a base counts as
#'   "rare" at its site (default 0.01).
#' @param allow_ambiguity If \code{FALSE} (default), gaps and IUPAC ambiguity
#'   codes in the query are flagged.
#' @param chimera_check Whether to run the sliding-window chimera screen
#'   (default \code{TRUE}).
#' @param pseudocount Smoothing pseudocount passed to
#'   \code{\link{build_malavi_site_profile}} when \code{site_profile} is
#'   \code{NULL}.
#' @param details If \code{TRUE}, the result also carries the full amino-acid
#'   translation, the per-site profile score, and the chimera result. Default
#'   \code{FALSE} keeps the object small.
#' @param genetic_code Genetic code for translation. Only NCBI code \code{4}
#'   (protozoan mitochondrial) is implemented; it is the correct code for avian
#'   haemosporidians.
#' @return An object of class \code{malavi_lineage_qc}: a list with
#'   \describe{
#'     \item{\code{call}}{the overall verdict: \code{known_lineage} (exact match),
#'       \code{plausible_new_lineage}, \code{review}, \code{strong_warning},
#'       \code{possible_error}, \code{possible_chimera}, or
#'       \code{invalid_or_strong_warning} / \code{invalid_sequence}.}
#'     \item{\code{score}}{the plausibility score in \code{[0, 1]} (1 = typical,
#'       0 = suspicious). Artifact risk is simply \code{1 - score}.}
#'     \item{\code{summary}}{a one-row data frame with \code{call}, \code{score},
#'       the \code{nearest_lineage} and \code{nearest_distance} (Hamming distance
#'       to it), \code{n_mutations} versus that lineage, of which
#'       \code{n_nonsynonymous} change the protein, and \code{n_stop_codons}.}
#'     \item{\code{flags}}{the character vector of warnings (see \emph{Flags}).}
#'     \item{\code{counts}}{a named integer vector with the full per-category
#'       counts behind the score (invariant-site changes, never-observed bases,
#'       rare bases, second-position changes, transversions).}
#'     \item{\code{nearest}}{a data frame of the five nearest known lineages and
#'       their distances.}
#'     \item{\code{mutations}}{one row per difference from the nearest lineage,
#'       annotated with codon position, amino-acid change, transition/transversion,
#'       and whether the query base is seen at that site in MalAvi.}
#'   }
#'   With \code{details = TRUE} the list also holds \code{translation},
#'   \code{site_profile_score}, and \code{chimera}.
#' @seealso \code{\link{amplicon_qc}}, \code{\link{lineage_screen}},
#'   \code{\link{build_malavi_site_profile}}, \code{\link{blast_malavi}}
#' @examples
#' ## screen a known lineage against the bundled alignment (expect a clean pass)
#' aln <- extract_alignment()
#' seq <- paste(as.character(aln[1, ]), collapse = "")
#' qc <- lineage_qc(seq)
#' qc
#' qc$summary
#' @export
lineage_qc <- function(query, reference = NULL, site_profile = NULL,
                       version = "latest",
                       expected_length = 479, rare_base_frequency = 0.01,
                       allow_ambiguity = FALSE, chimera_check = TRUE,
                       pseudocount = 0.01, details = FALSE, genetic_code = 4) {
  if (genetic_code != 4) {
    stop("lineage_qc() currently implements only genetic_code = 4 ",
         "(protozoan mitochondrial code), which is correct for haemosporidians.",
         call. = FALSE)
  }

  settings <- .lineage_qc_settings(expected_length, rare_base_frequency)

  ## build (or reuse) the coded reference and its site profile once
  charmat   <- .qc_char_matrix(reference, version)
  refcode   <- .qc_code_matrix(charmat)
  ref_names <- rownames(charmat)
  if (is.null(site_profile)) {
    site_profile <- .qc_site_profile(charmat, pseudocount = pseudocount)
  }
  code <- .qc_genetic_code_4()

  .lineage_qc_core(query, charmat, refcode, ref_names, site_profile, code,
                   allow_ambiguity, chimera_check, settings, details)
}

## Core lineage QC, working from a pre-built coded reference and site profile so
## amplicon_qc() can call it per variant without re-coding the reference.
## `settings` is the merged list from .lineage_qc_settings(). All arguments are
## already validated/prepared by the callers.
.lineage_qc_core <- function(query, charmat, refcode, ref_names, site_profile,
                             code, allow_ambiguity, chimera_check, settings,
                             details) {
  query  <- .qc_clean_seq(query)
  qchars <- strsplit(query, "", fixed = TRUE)[[1]]
  expected_length <- settings$expected_length

  flags <- character(0)
  query_length <- nchar(query)

  ## ---- basic sequence validation ----
  if (query_length != expected_length) {
    flags <- c(flags, paste0("wrong_length_expected_", expected_length,
                             "_observed_", query_length))
  }
  valid_chars <- if (allow_ambiguity) "^[ACGTRYSWKMBDHVN-]+$" else "^[ACGT]+$"
  if (!grepl(valid_chars, query)) flags <- c(flags, "invalid_or_disallowed_characters")
  if (grepl("-", query, fixed = TRUE)) flags <- c(flags, "contains_gaps")
  if (grepl("N", query, fixed = TRUE)) flags <- c(flags, "contains_N")
  if (grepl("[RYSWKMBDHV]", query))    flags <- c(flags, "contains_ambiguity_codes")

  ## If the length does not match the reference, the per-site and per-codon
  ## metrics cannot be computed. Return a useful invalid-sequence result.
  if (query_length != nrow(site_profile)) {
    out <- list(
      summary = data.frame(call = "invalid_sequence", score = 0,
                           nearest_lineage = NA_character_, nearest_distance = NA_real_,
                           n_mutations = NA_integer_, n_nonsynonymous = NA_integer_,
                           n_stop_codons = NA_integer_, stringsAsFactors = FALSE),
      call = "invalid_sequence", score = 0, flags = unique(flags),
      query_length = query_length,
      message = paste("Query length does not match the expected MalAvi barcode",
                      "length; downstream QC was skipped. Align the query to the",
                      "reference first (see blast_malavi).")
    )
    class(out) <- c("malavi_lineage_qc", class(out))
    return(out)
  }

  ## ---- translation / stop codons ----
  aa <- .qc_translate(qchars, code)
  n_stop_codons <- sum(aa == "*")
  has_stop_codon <- n_stop_codons > 0
  has_unknown_aa <- any(aa == "X")
  if (has_stop_codon) flags <- c(flags, "contains_stop_codon")
  if (has_unknown_aa) flags <- c(flags, "contains_unknown_amino_acid_after_translation")

  ## ---- nearest known lineage ----
  qcode   <- .qc_code_vec(qchars)
  nearest <- .qc_nearest(qcode, refcode, ref_names, top_n = 5L)
  min_distance <- nearest$distance[1]
  if (min_distance == 0) {
    flags <- c(flags, "exact_match_to_known_lineage")
  } else if (min_distance <= settings$near_known_distance) {
    flags <- c(flags, "near_known_lineage")
  } else if (min_distance <= settings$divergent_distance) {
    flags <- c(flags, "moderately_divergent_from_known_lineages")
  } else {
    flags <- c(flags, "highly_divergent_from_known_lineages")
  }

  ## ---- site-profile surprise score ----
  profile_score <- .qc_score_site(qchars, site_profile,
                                  rare_freq = settings$rare_base_frequency)
  site_flags <- profile_score$site_flags
  n_invariant_changes     <- sum(site_flags == "invariant_site_change")
  n_unobserved_site_bases <- sum(site_flags == "base_never_observed_at_site")
  n_rare_site_bases       <- sum(site_flags == "rare_base_at_site")
  if (n_invariant_changes > 0)
    flags <- c(flags, paste0(n_invariant_changes, "_changes_at_invariant_sites"))
  if (n_unobserved_site_bases > 0)
    flags <- c(flags, paste0(n_unobserved_site_bases, "_bases_never_observed_at_their_sites"))
  if (n_rare_site_bases > 0)
    flags <- c(flags, paste0(n_rare_site_bases, "_rare_bases_at_sites"))

  ## ---- mutations relative to the nearest lineage ----
  nearest_chars  <- charmat[nearest$index[1], ]
  mutations      <- .qc_annotate_mutations(qchars, nearest_chars, site_profile, code)
  n_mutations       <- nrow(mutations)
  n_nonsynonymous   <- if (n_mutations) sum(mutations$nonsynonymous) else 0L
  n_second_position <- if (n_mutations) sum(mutations$codon_position == 2) else 0L
  n_transversions   <- if (n_mutations) sum(mutations$transversion) else 0L
  if (n_nonsynonymous > 0)
    flags <- c(flags, paste0(n_nonsynonymous, "_nonsynonymous_changes_vs_nearest_lineage"))
  if (n_second_position > 0)
    flags <- c(flags, paste0(n_second_position, "_second_codon_position_changes_vs_nearest_lineage"))
  if (n_transversions > 0)
    flags <- c(flags, paste0(n_transversions, "_transversions_vs_nearest_lineage"))

  ## ---- chimera screen ----
  chimera <- NULL
  chimera_flagged <- FALSE
  if (chimera_check) {
    chimera <- .qc_detect_chimera(qcode, refcode, ref_names,
                                  window = settings$chimera_window,
                                  step = settings$chimera_step)
    chimera_flagged <- chimera$chimera_delta >= settings$chimera_delta_threshold &&
      chimera$parent_switches >= settings$chimera_min_parent_switches
    if (chimera_flagged) flags <- c(flags, "possible_chimera_or_mixed_template_pattern")
  }

  ## ---- transparent penalty -> bounded score ----
  ## Start at penalty 0 and add interpretable amounts. This is an empirical
  ## plausibility score, not a statistical probability.
  penalty <- settings$invariant_site_penalty * n_invariant_changes +
    settings$unobserved_base_penalty * n_unobserved_site_bases +
    settings$rare_base_penalty * n_rare_site_bases +
    settings$nonsynonymous_penalty * n_nonsynonymous +
    settings$second_position_penalty * n_second_position +
    settings$transversion_penalty * n_transversions
  if (has_stop_codon) penalty <- penalty + 10
  if (any(c("contains_gaps", "contains_N", "contains_ambiguity_codes",
            "invalid_or_disallowed_characters") %in% flags)) penalty <- penalty + 2
  if (chimera_flagged) penalty <- penalty + 4

  score <- max(0, min(1, exp(-penalty / 10)))

  ## ---- final call ----
  if (has_stop_codon) {
    call <- "invalid_or_strong_warning"
  } else if ("exact_match_to_known_lineage" %in% flags) {
    call <- "known_lineage"
  } else if (chimera_flagged) {
    call <- "possible_chimera"
  } else if (score >= settings$pass_score) {
    call <- "plausible_new_lineage"
  } else if (score >= settings$review_score) {
    call <- "review"
  } else if (score >= settings$strong_warning_score) {
    call <- "strong_warning"
  } else {
    call <- "possible_error"
  }

  summary <- data.frame(
    call = call, score = score,
    nearest_lineage = nearest$lineage[1], nearest_distance = min_distance,
    n_mutations = n_mutations, n_nonsynonymous = n_nonsynonymous,
    n_stop_codons = n_stop_codons, stringsAsFactors = FALSE
  )

  ## full per-category counts kept out of the printed summary to reduce clutter
  counts <- c(n_invariant_site_changes = n_invariant_changes,
              n_bases_never_observed   = n_unobserved_site_bases,
              n_rare_site_bases        = n_rare_site_bases,
              n_second_position_changes = n_second_position,
              n_transversions          = n_transversions)

  result <- list(summary = summary, call = call, score = score,
                 flags = unique(flags), counts = counts,
                 nearest = nearest, mutations = mutations)
  if (details) {
    result$translation <- list(amino_acid_sequence = paste(aa, collapse = ""),
                               n_stop_codons = n_stop_codons,
                               has_stop_codon = has_stop_codon,
                               has_unknown_aa = has_unknown_aa)
    result$site_profile_score <- profile_score
    result$chimera <- chimera
  }
  class(result) <- c("malavi_lineage_qc", class(result))
  result
}

#' @export
print.malavi_lineage_qc <- function(x, ...) {
  cat("MalAvi lineage QC\n")
  cat("  call:   ", x$call, "\n", sep = "")
  cat("  score:  ", formatC(x$score, format = "f", digits = 2),
      "   (0 = suspicious, 1 = typical of known MalAvi diversity)\n", sep = "")
  if (identical(x$call, "invalid_sequence")) {
    cat("\n", x$message, "\n", sep = "")
  } else {
    s <- x$summary
    cat("  nearest lineage: ", s$nearest_lineage,
        "  (distance ", s$nearest_distance, ")\n", sep = "")
    cat("  mutations vs nearest: ", s$n_mutations,
        "  (", s$n_nonsynonymous, " nonsynonymous, ",
        s$n_stop_codons, " stop codons)\n", sep = "")
  }
  cat("\nFlags:\n")
  if (length(x$flags) == 0) {
    cat("  none\n")
  } else {
    cat(paste0("  - ", x$flags, collapse = "\n"), "\n")
  }
  invisible(x)
}
