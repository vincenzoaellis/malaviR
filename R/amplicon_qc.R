## Abundance-aware quality-control screen for denoised amplicon sequence
## variants (ASVs). Builds on lineage_qc() per variant and adds abundance flags.

#' Default thresholds for \code{amplicon_qc}
#'
#' These are new experimental functions, that we're still testing. Please use
#' with caution.
#'
#' Returns the list of abundance-related thresholds used by
#' \code{\link{amplicon_qc}}. Edit the returned list and pass it back via the
#' \code{amplicon_thresholds} argument to customize the screen.
#'
#' @return A named list with:
#'   \describe{
#'     \item{\code{min_relative_frequency}}{variants below this within-sample
#'       relative frequency are flagged as below the minimum.}
#'     \item{\code{low_frequency_warning}, \code{very_low_frequency_warning}}{
#'       relative-frequency cutoffs for the low / very-low frequency flags.}
#'     \item{\code{parent_fold_abundance}}{fold-abundance of a more common
#'       variant above which a near, rarer variant looks like an error
#'       derivative.}
#'     \item{\code{oneoff_distance}, \code{twooff_distance}}{base distances from
#'       a much more abundant variant used by the error-derivative flags.}
#'     \item{\code{index_hop_relative_frequency}}{relative frequency below which
#'       a variant may reflect index hopping or low-level contamination.}
#'   }
#' @seealso \code{\link{amplicon_qc}}
#' @examples
#' default_amplicon_qc_thresholds()
#' @export
default_amplicon_qc_thresholds <- function() {
  list(
    min_relative_frequency = 0.01,
    low_frequency_warning = 0.01,
    very_low_frequency_warning = 0.001,

    ## a low-frequency variant a base or two from a much more abundant variant
    ## may be an error derivative of it
    parent_fold_abundance = 10,
    oneoff_distance = 1,
    twooff_distance = 2,

    ## possible low-level contamination / index hopping
    index_hop_relative_frequency = 0.001
  )
}

#' Quality-control screen for denoised amplicon sequence variants
#'
#' These are new experimental functions, that we're still testing. Please use
#' with caution.
#'
#' An abundance-aware QC screen for amplicon sequence variants (ASVs) produced by
#' denoising tools such as \pkg{dada2}, vsearch, or unoise. Each variant is run
#' through \code{\link{lineage_qc}} for biological plausibility, and additional
#' flags use the read counts to catch likely technical artifacts: variants below
#' a minimum relative abundance, variants a base or two away from a much more
#' abundant variant (possible error derivatives), and very rare variants that
#' may reflect index hopping or low-level contamination.
#'
#' As with \code{\link{lineage_qc}}, the variant sequences must already be
#' aligned to the MalAvi barcode (same length and frame as the reference). The
#' reference is coded and its site profile built once, then reused across all
#' variants.
#'
#' @param variants A \code{data.frame} of variants with at least a sequence
#'   column and a count column (named by \code{sequence_col} / \code{count_col}).
#'   An optional \code{sample_col} runs the screen separately per sample;
#'   relative frequencies and parent comparisons are then computed within each
#'   sample.
#' @param reference Reference alignment of curated lineages: a \code{DNAbin}
#'   alignment, a named character vector, or \code{NULL} (default) to use the
#'   bundled MalAvi alignment for \code{version}.
#' @param site_profile Optional precomputed profile from
#'   \code{\link{build_malavi_site_profile}}; built from \code{reference} if
#'   \code{NULL}.
#' @param version MalAvi release to use when \code{reference} is \code{NULL}.
#' @param sequence_col,count_col Column names in \code{variants} holding the
#'   sequence and the read count.
#' @param sample_col Optional column name grouping variants into samples. If
#'   \code{NULL} (default), all variants are treated as one pool.
#' @param lineage_qc_thresholds,amplicon_thresholds Threshold lists; see
#'   \code{\link{default_lineage_qc_thresholds}} and
#'   \code{\link{default_amplicon_qc_thresholds}}.
#' @param allow_ambiguity,chimera_check Passed through to \code{\link{lineage_qc}}
#'   for each variant.
#' @return An object of class \code{malavi_amplicon_qc}: the input
#'   \code{data.frame} with added columns including \code{relative_frequency},
#'   the nearest more-abundant variant and distance to it,
#'   \code{lineage_call} / \code{lineage_score} / \code{artifact_risk} and
#'   nearest MalAvi lineage from \code{\link{lineage_qc}}, the per-variant
#'   \code{amplicon_flags}, and an overall \code{amplicon_call}.
#' @seealso \code{\link{lineage_qc}}, \code{\link{default_amplicon_qc_thresholds}}
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
                        version = "latest", sequence_col = "sequence",
                        count_col = "count", sample_col = NULL,
                        lineage_qc_thresholds = default_lineage_qc_thresholds(),
                        amplicon_thresholds = default_amplicon_qc_thresholds(),
                        allow_ambiguity = FALSE, chimera_check = TRUE) {
  if (!sequence_col %in% names(variants)) {
    stop("`sequence_col` ('", sequence_col, "') not found in variants.", call. = FALSE)
  }
  if (!count_col %in% names(variants)) {
    stop("`count_col` ('", count_col, "') not found in variants.", call. = FALSE)
  }

  variants[[sequence_col]] <- vapply(variants[[sequence_col]], .qc_clean_seq, character(1))

  ## build the coded reference + site profile once, shared across all variants
  charmat <- .qc_char_matrix(reference, version)
  refcode <- .qc_code_matrix(charmat)
  ref_names <- rownames(charmat)
  if (is.null(site_profile)) site_profile <- .qc_site_profile(charmat)
  code <- .qc_genetic_code_4()
  min_rf <- amplicon_thresholds$min_relative_frequency

  ## with no sample column, treat all variants as one pseudo-sample
  if (is.null(sample_col)) {
    variants$.malavi_qc_sample <- "all"
    sample_col <- ".malavi_qc_sample"
  }

  out_list <- lapply(unique(variants[[sample_col]]), function(s) {
    dat <- variants[variants[[sample_col]] == s, , drop = FALSE]

    total_reads <- sum(dat[[count_col]], na.rm = TRUE)
    dat$relative_frequency <- dat[[count_col]] / total_reads
    ## sort most abundant first so each variant's "parents" are the rows above it
    dat <- dat[order(dat[[count_col]], decreasing = TRUE), , drop = FALSE]

    seqs   <- dat[[sequence_col]]
    counts <- dat[[count_col]]
    ## integer-code the variants once for fast pairwise distances
    seq_codes <- lapply(strsplit(seqs, "", fixed = TRUE), .qc_code_vec)

    ## nearest more-abundant variant for each variant (NA for the most abundant)
    parent <- .amplicon_nearest_parent(seq_codes, counts)
    dat <- cbind(dat, parent)

    ## biological plausibility of each variant via lineage_qc()
    lineage_results <- lapply(seqs, function(q) {
      .lineage_qc_core(q, charmat, refcode, ref_names, site_profile, code,
                       allow_ambiguity, chimera_check, lineage_qc_thresholds,
                       return_details = FALSE)
    })
    dat$lineage_call  <- vapply(lineage_results, function(x) x$call, character(1))
    dat$lineage_score <- vapply(lineage_results, function(x) x$overall_score, numeric(1))
    dat$artifact_risk <- vapply(lineage_results, function(x) x$artifact_risk, numeric(1))
    dat$nearest_malavi_lineage <- vapply(lineage_results,
                                         function(x) x$nearest$lineage[1], character(1))
    dat$distance_to_nearest_malavi <- vapply(lineage_results,
                                             function(x) x$nearest$distance[1], numeric(1))
    dat$lineage_flags <- vapply(lineage_results,
                                function(x) paste(x$flags, collapse = ";"), character(1))

    ## abundance-aware amplicon flags + an overall call, per variant
    dat$amplicon_flags <- vapply(seq_len(nrow(dat)),
      function(i) .amplicon_flags(dat[i, ], amplicon_thresholds, min_rf),
      character(1))
    dat$amplicon_call <- vapply(seq_len(nrow(dat)),
      function(i) .amplicon_call(dat$amplicon_flags[i], dat$lineage_call[i]),
      character(1))

    dat$sample_total_reads <- total_reads
    dat
  })

  out <- do.call(rbind, out_list)
  out$.malavi_qc_sample <- NULL
  rownames(out) <- NULL

  attr(out, "amplicon_thresholds") <- amplicon_thresholds
  attr(out, "lineage_qc_thresholds") <- lineage_qc_thresholds
  class(out) <- c("malavi_amplicon_qc", class(out))
  out
}

## For abundance-sorted variants (most abundant first), find each variant's
## nearest more-abundant variant by Hamming distance and the fold-abundance of
## that parent. The first (most abundant) variant has no parent. `seq_codes` is
## a list of integer-coded variants aligned with `counts`.
.amplicon_nearest_parent <- function(seq_codes, counts) {
  n <- length(seq_codes)
  rows <- lapply(seq_len(n), function(i) {
    if (i == 1L) {
      return(data.frame(nearest_more_abundant_count = NA_real_,
                        distance_to_nearest_more_abundant = NA_real_,
                        fold_less_abundant_than_parent = NA_real_,
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
    data.frame(nearest_more_abundant_count = counts[best],
               distance_to_nearest_more_abundant = dists[best],
               fold_less_abundant_than_parent = counts[best] / counts[i],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

## Build the abundance-aware flag string for a single variant row.
.amplicon_flags <- function(row, th, min_rf) {
  flags <- character(0)
  rf <- row$relative_frequency

  if (rf < min_rf) flags <- c(flags, paste0("below_min_relative_frequency_", min_rf))
  if (rf < th$very_low_frequency_warning) {
    flags <- c(flags, "very_low_frequency_variant")
  } else if (rf < th$low_frequency_warning) {
    flags <- c(flags, "low_frequency_variant")
  }

  dist_parent <- row$distance_to_nearest_more_abundant
  fold_parent <- row$fold_less_abundant_than_parent
  if (!is.na(dist_parent) && !is.na(fold_parent)) {
    if (dist_parent <= th$oneoff_distance && fold_parent >= th$parent_fold_abundance) {
      flags <- c(flags, "oneoff_from_much_more_abundant_variant")
    }
    if (dist_parent <= th$twooff_distance && fold_parent >= th$parent_fold_abundance &&
        rf < min_rf) {
      flags <- c(flags, "low_frequency_near_abundant_variant_possible_error_derivative")
    }
  }

  if (rf < th$index_hop_relative_frequency) {
    flags <- c(flags, "possible_index_hopping_or_low_level_contamination")
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
  if (grepl("possible_index_hopping|oneoff_from_much_more_abundant|possible_error_derivative",
            flags)) {
    return("possible_amplicon_artifact")
  }
  if (grepl("below_min_relative_frequency|low_frequency_variant", flags)) {
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
  cat("------------------\n")
  cat("Variants:", nrow(x), "\n")
  th <- attr(x, "amplicon_thresholds")
  cat("Minimum relative frequency:", th$min_relative_frequency, "\n\n")

  cols <- intersect(
    c("sequence", "count", "relative_frequency", "amplicon_call", "lineage_call",
      "lineage_score", "nearest_malavi_lineage", "distance_to_nearest_malavi",
      "amplicon_flags"),
    names(x)
  )
  print(as.data.frame(x)[, cols, drop = FALSE])
  invisible(x)
}
