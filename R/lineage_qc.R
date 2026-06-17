## Quality-control screen for a single candidate MalAvi cytochrome b lineage.
## See internal.R for the shared, dependency-light QC helpers used here.

#' Default thresholds for \code{lineage_qc}
#'
#' These are new experimental functions, that we're still testing. Please use
#' with caution.
#'
#' Returns the list of tunable thresholds and penalties used by
#' \code{\link{lineage_qc}}. Edit the returned list and pass it back via the
#' \code{thresholds} argument to customize the screen (for example, raise or
#' lower \code{rare_base_frequency} to change what counts as a rare base at a
#' site).
#'
#' @return A named list of thresholds. The elements are:
#'   \describe{
#'     \item{\code{expected_length}}{expected barcode length in bp (479).}
#'     \item{\code{invariant_site_penalty}, \code{unobserved_base_penalty},
#'       \code{rare_base_penalty}}{penalty weights for changes at invariant
#'       sites, bases never observed at a site, and rare bases at a site.}
#'     \item{\code{rare_base_frequency}}{smoothed frequency below which a base is
#'       "rare" at a site (default 0.01).}
#'     \item{\code{nonsynonymous_penalty}, \code{second_position_penalty},
#'       \code{transversion_penalty}}{penalty weights for those mutation types
#'       relative to the nearest lineage.}
#'     \item{\code{near_known_distance}, \code{divergent_distance}}{Hamming-
#'       distance cutoffs binning the query as near / moderately divergent /
#'       highly divergent from known lineages.}
#'     \item{\code{chimera_window}, \code{chimera_step},
#'       \code{chimera_delta_threshold}, \code{chimera_min_parent_switches}}{the
#'       sliding-window chimera screen settings.}
#'     \item{\code{pass_score}, \code{review_score},
#'       \code{strong_warning_score}}{score cutoffs for the final call.}
#'   }
#' @seealso \code{\link{lineage_qc}}
#' @examples
#' th <- default_lineage_qc_thresholds()
#' th$rare_base_frequency <- 0.02   # treat more bases as "rare"
#' @export
default_lineage_qc_thresholds <- function() {
  list(
    expected_length = 479,

    ## site-profile penalties
    invariant_site_penalty = 4,
    unobserved_base_penalty = 3,
    rare_base_penalty = 1.5,
    rare_base_frequency = 0.01,

    ## mutation-pattern penalties (relative to the nearest lineage)
    nonsynonymous_penalty = 1,
    second_position_penalty = 1.5,
    transversion_penalty = 0.75,

    ## distance bins
    near_known_distance = 2,
    divergent_distance = 5,

    ## sliding-window chimera screen
    chimera_window = 120,
    chimera_step = 20,
    chimera_delta_threshold = 3,
    chimera_min_parent_switches = 2,

    ## final-score cutoffs for the call
    pass_score = 0.85,
    review_score = 0.60,
    strong_warning_score = 0.35
  )
}

#' Build a per-site base profile of a MalAvi alignment
#'
#' These are new experimental functions, that we're still testing. Please use
#' with caution.
#'
#' Summarizes a curated MalAvi alignment one column (site) at a time: per-base
#' counts, smoothed frequencies, the major base, which bases are observed,
#' whether the site is invariant, its Shannon entropy, and its codon position
#' (assuming reading frame 1). This profile is what \code{\link{lineage_qc}} and
#' \code{\link{amplicon_qc}} use to decide whether a query's bases are typical of
#' known MalAvi diversity.
#'
#' Building the profile from the bundled alignment takes about a second. If you
#' will screen many sequences, build it once and pass it to \code{lineage_qc} /
#' \code{amplicon_qc} via their \code{site_profile} argument to avoid rebuilding
#' it on every call.
#'
#' @param reference A reference alignment: a \code{DNAbin} alignment (e.g. from
#'   \code{\link{extract_alignment}}), a named character vector of equal-length
#'   aligned sequences, or \code{NULL} (default) to use the bundled MalAvi
#'   alignment for \code{version}.
#' @param version MalAvi release to use when \code{reference} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @param pseudocount Smoothing pseudocount added to each base count when
#'   computing per-site frequencies (default 0.5).
#' @return A \code{data.frame} with one row per alignment position and columns
#'   \code{position}, \code{codon_position}, \code{n_nonmissing},
#'   \code{count_A/C/G/T}, \code{freq_A/C/G/T}, \code{major_base},
#'   \code{n_observed_alleles}, \code{observed_alleles}, \code{invariant}, and
#'   \code{entropy}.
#' @seealso \code{\link{lineage_qc}}, \code{\link{amplicon_qc}}
#' @examples
#' profile <- build_malavi_site_profile()
#' head(profile)
#' @export
build_malavi_site_profile <- function(reference = NULL, version = "latest",
                                      pseudocount = 0.5) {
  charmat <- .qc_char_matrix(reference, version)
  .qc_site_profile(charmat, pseudocount = pseudocount)
}

#' Quality-control screen for a candidate MalAvi lineage
#'
#' These are new experimental functions, that we're still testing. Please use
#' with caution.
#'
#' A lightweight plausibility screen for a single, aligned 479 bp MalAvi
#' cytochrome \emph{b} barcode. It does \strong{not} decide whether a lineage is
#' real; it flags features that are surprising relative to curated MalAvi
#' diversity so that they can be reviewed by hand. Warnings should be read as
#' "manual review recommended", not "this is definitely an error".
#'
#' The query must already be aligned to the MalAvi barcode (same length and frame
#' as the reference, 479 bp by default). For an unaligned sequence, find its
#' closest lineages first with \code{\link{blast_malavi}}. The screen checks:
#' sequence length, gaps / ambiguity codes, stop codons (translated in frame 1
#' under the protozoan mitochondrial genetic code, NCBI code 4), Hamming distance
#' to the nearest known lineage, changes at invariant or rarely varying sites,
#' bases never seen at a site, unusual nonsynonymous / second-codon-position /
#' transversion changes relative to the nearest lineage, and a crude
#' sliding-window chimera pattern.
#'
#' The penalties are summed and mapped to an \code{overall_score} in
#' \code{[0, 1]} (1 = very typical of known MalAvi diversity, 0 = highly
#' suspicious or invalid), with \code{artifact_risk = 1 - overall_score}. The
#' score is a transparent, rule-based heuristic, \strong{not} a calibrated
#' probability; treat it as a triage aid and adjust the weights and cutoffs via
#' \code{thresholds} (see \code{\link{default_lineage_qc_thresholds}}).
#'
#' @param query A single candidate barcode as a character string. Whitespace is
#'   removed and the sequence is upper-cased.
#' @param reference Reference alignment of curated lineages: a \code{DNAbin}
#'   alignment, a named character vector, or \code{NULL} (default) to use the
#'   bundled MalAvi alignment for \code{version}.
#' @param site_profile Optional precomputed profile from
#'   \code{\link{build_malavi_site_profile}}. If \code{NULL} (default) it is built
#'   from \code{reference}. Supply it to avoid rebuilding when screening many
#'   sequences.
#' @param version MalAvi release to use when \code{reference} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @param genetic_code Genetic code for translation. Only NCBI code \code{4}
#'   (protozoan mitochondrial) is implemented and is correct for avian
#'   haemosporidians.
#' @param allow_ambiguity If \code{FALSE} (default), gaps and IUPAC ambiguity
#'   codes in the query are flagged.
#' @param pseudocount Smoothing pseudocount passed to
#'   \code{build_malavi_site_profile} when \code{site_profile} is \code{NULL}.
#' @param chimera_check Whether to run the sliding-window chimera screen
#'   (default \code{TRUE}).
#' @param thresholds A list of thresholds and penalties; see
#'   \code{\link{default_lineage_qc_thresholds}}.
#' @param return_details If \code{TRUE} (default), the result also carries the
#'   full translation, site-profile score, chimera result, and thresholds used.
#' @return An object of class \code{malavi_lineage_qc}: a list with a one-row
#'   \code{summary} data frame, the \code{call}, \code{overall_score},
#'   \code{artifact_risk}, character vector of \code{flags}, the \code{nearest}
#'   lineages table, and the \code{mutation_table} of differences from the
#'   nearest lineage (plus \code{translation}, \code{site_profile_score},
#'   \code{chimera}, and \code{thresholds} when \code{return_details = TRUE}).
#' @seealso \code{\link{amplicon_qc}}, \code{\link{build_malavi_site_profile}},
#'   \code{\link{blast_malavi}}, \code{\link{default_lineage_qc_thresholds}}
#' @examples
#' ## screen a known lineage against the bundled alignment (expect a clean pass)
#' aln <- extract_alignment()
#' seq <- paste(as.character(aln[1, ]), collapse = "")
#' qc <- lineage_qc(seq)
#' qc
#' qc$summary$call
#' @export
lineage_qc <- function(query, reference = NULL, site_profile = NULL,
                       version = "latest", genetic_code = 4,
                       allow_ambiguity = FALSE, pseudocount = 0.5,
                       chimera_check = TRUE,
                       thresholds = default_lineage_qc_thresholds(),
                       return_details = TRUE) {
  if (genetic_code != 4) {
    stop("lineage_qc() currently implements only genetic_code = 4 ",
         "(protozoan mitochondrial code), which is correct for haemosporidians.",
         call. = FALSE)
  }

  ## build (or reuse) the coded reference and its site profile once
  charmat <- .qc_char_matrix(reference, version)
  refcode <- .qc_code_matrix(charmat)
  ref_names <- rownames(charmat)
  if (is.null(site_profile)) {
    site_profile <- .qc_site_profile(charmat, pseudocount = pseudocount)
  }
  code <- .qc_genetic_code_4()

  res <- .lineage_qc_core(query, charmat, refcode, ref_names, site_profile,
                          code, allow_ambiguity, chimera_check, thresholds,
                          return_details)
  res
}

## Core lineage QC, working from a pre-built coded reference and site profile so
## amplicon_qc() can call it per variant without re-coding the reference. All
## arguments are already validated/prepared by the callers.
.lineage_qc_core <- function(query, charmat, refcode, ref_names, site_profile,
                             code, allow_ambiguity, chimera_check, thresholds,
                             return_details) {
  query <- .qc_clean_seq(query)
  qchars <- strsplit(query, "", fixed = TRUE)[[1]]
  expected_length <- thresholds$expected_length

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
      summary = data.frame(call = "invalid_sequence", overall_score = 0,
                           artifact_risk = 1, query_length = query_length,
                           stringsAsFactors = FALSE),
      call = "invalid_sequence", overall_score = 0, artifact_risk = 1,
      flags = unique(flags), query_length = query_length,
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
  qcode <- .qc_code_vec(qchars)
  nearest <- .qc_nearest(qcode, refcode, ref_names, top_n = 5L)
  min_distance <- nearest$distance[1]
  if (min_distance == 0) {
    flags <- c(flags, "exact_match_to_known_lineage")
  } else if (min_distance <= thresholds$near_known_distance) {
    flags <- c(flags, "near_known_lineage")
  } else if (min_distance <= thresholds$divergent_distance) {
    flags <- c(flags, "moderately_divergent_from_known_lineages")
  } else {
    flags <- c(flags, "highly_divergent_from_known_lineages")
  }

  ## ---- site-profile surprise score ----
  profile_score <- .qc_score_site(qchars, site_profile,
                                  rare_freq = thresholds$rare_base_frequency)
  site_flags <- profile_score$site_flags
  n_invariant_changes      <- sum(site_flags == "invariant_site_change")
  n_unobserved_site_bases  <- sum(site_flags == "base_never_observed_at_site")
  n_rare_site_bases        <- sum(site_flags == "rare_base_at_site")
  if (n_invariant_changes > 0)
    flags <- c(flags, paste0(n_invariant_changes, "_changes_at_invariant_sites"))
  if (n_unobserved_site_bases > 0)
    flags <- c(flags, paste0(n_unobserved_site_bases, "_bases_never_observed_at_their_sites"))
  if (n_rare_site_bases > 0)
    flags <- c(flags, paste0(n_rare_site_bases, "_rare_bases_at_sites"))

  ## ---- mutations relative to the nearest lineage ----
  nearest_chars <- charmat[nearest$index[1], ]
  mutation_table <- .qc_annotate_mutations(qchars, nearest_chars, site_profile, code)
  n_mutations       <- nrow(mutation_table)
  n_nonsynonymous   <- if (n_mutations) sum(mutation_table$nonsynonymous) else 0
  n_second_position <- if (n_mutations) sum(mutation_table$codon_position == 2) else 0
  n_transversions   <- if (n_mutations) sum(mutation_table$transversion) else 0
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
                                  window = thresholds$chimera_window,
                                  step = thresholds$chimera_step)
    chimera_flagged <- chimera$chimera_delta >= thresholds$chimera_delta_threshold &&
      chimera$parent_switches >= thresholds$chimera_min_parent_switches
    if (chimera_flagged) flags <- c(flags, "possible_chimera_or_mixed_template_pattern")
  }

  ## ---- transparent penalty -> bounded score ----
  ## Start at penalty 0 and add interpretable amounts. This is an empirical
  ## plausibility score, not a statistical probability.
  penalty <- thresholds$invariant_site_penalty * n_invariant_changes +
    thresholds$unobserved_base_penalty * n_unobserved_site_bases +
    thresholds$rare_base_penalty * n_rare_site_bases +
    thresholds$nonsynonymous_penalty * n_nonsynonymous +
    thresholds$second_position_penalty * n_second_position +
    thresholds$transversion_penalty * n_transversions
  if (has_stop_codon) penalty <- penalty + 10
  if (any(c("contains_gaps", "contains_N", "contains_ambiguity_codes",
            "invalid_or_disallowed_characters") %in% flags)) penalty <- penalty + 2
  if (chimera_flagged) penalty <- penalty + 4

  overall_score <- max(0, min(1, exp(-penalty / 10)))
  artifact_risk <- 1 - overall_score

  ## ---- final call ----
  if (has_stop_codon) {
    call <- "invalid_or_strong_warning"
  } else if ("exact_match_to_known_lineage" %in% flags) {
    call <- "known_lineage"
  } else if (chimera_flagged) {
    call <- "possible_chimera"
  } else if (overall_score >= thresholds$pass_score) {
    call <- "plausible_new_lineage"
  } else if (overall_score >= thresholds$review_score) {
    call <- "review"
  } else if (overall_score >= thresholds$strong_warning_score) {
    call <- "strong_warning"
  } else {
    call <- "possible_error"
  }

  summary <- data.frame(
    call = call, overall_score = overall_score, artifact_risk = artifact_risk,
    nearest_lineage = nearest$lineage[1], nearest_distance = min_distance,
    n_mutations_vs_nearest = n_mutations,
    n_invariant_site_changes = n_invariant_changes,
    n_unobserved_site_bases = n_unobserved_site_bases,
    n_rare_site_bases = n_rare_site_bases, n_nonsynonymous = n_nonsynonymous,
    n_second_position_changes = n_second_position, n_transversions = n_transversions,
    site_log_likelihood = profile_score$log_likelihood,
    site_mean_log_probability = profile_score$mean_log_probability,
    n_stop_codons = n_stop_codons, stringsAsFactors = FALSE
  )

  result <- list(summary = summary, call = call, overall_score = overall_score,
                 artifact_risk = artifact_risk, flags = unique(flags),
                 nearest = nearest, mutation_table = mutation_table)
  if (return_details) {
    result$translation <- list(amino_acid_sequence = paste(aa, collapse = ""),
                               n_stop_codons = n_stop_codons,
                               has_stop_codon = has_stop_codon,
                               has_unknown_aa = has_unknown_aa)
    result$site_profile_score <- profile_score
    result$chimera <- chimera
    result$thresholds <- thresholds
  }
  class(result) <- c("malavi_lineage_qc", class(result))
  result
}

#' @export
print.malavi_lineage_qc <- function(x, ...) {
  cat("MalAvi lineage QC\n")
  cat("-----------------\n")
  print(x$summary)
  cat("\nFlags:\n")
  if (length(x$flags) == 0) {
    cat("  none\n")
  } else {
    cat(paste0("  - ", x$flags, collapse = "\n"), "\n")
  }
  invisible(x)
}
