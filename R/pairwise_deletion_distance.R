## Public wrapper around the shared pairwise-deletion Hamming machinery in
## internal.R (.qc_nearest). Exposes per-reference distance AND the number of
## comparable positions, so callers (e.g. within-sample ASV reconciliation) can
## measure sequence distances without letting undetermined positions inflate them.

#' Pairwise-deletion Hamming distance from a query to reference sequences
#'
#' \strong{Experimental.} Counts base differences between one aligned query
#' sequence and each of a set of aligned reference sequences using
#' \strong{pairwise deletion}: a position is compared only where \emph{both} the
#' query and the reference carry a definite A/C/G/T base. Positions that are an
#' ambiguity code (e.g. \code{N}) or a gap in \emph{either} sequence are skipped --
#' they are not counted as differences and not counted as comparable. Each
#' reference therefore reports both its \code{distance} (mismatches among the
#' positions actually compared) and \code{n_comparable} (how many positions that
#' was), so a small distance backed by few comparable positions is not mistaken
#' for strong agreement.
#'
#' This is the same distance used internally by \code{\link{lineage_qc}} and the
#' near-neighbor search, exposed for reuse. Typical uses: how far a denoised ASV
#' sits from the dominant ASV in the same sample (pass that sample's ASVs as
#' \code{reference}), or from candidate MalAvi lineages (the default reference).
#'
#' Results are ordered by ascending \code{distance}, then by \strong{descending}
#' \code{n_comparable}, so among equally distant references the most completely
#' overlapping one comes first. This is deliberate: a fully-overlapping exact match
#' is never displaced by a reference that merely ties it because an ambiguous
#' position was dropped from the comparison (the same tie-break
#' \code{\link{lineage_qc}} relies on).
#'
#' @param query A single aligned sequence (character scalar), the same length as
#'   the reference alignment (479 bp for the bundled MalAvi data). Case and
#'   whitespace are normalized; ambiguity codes and gaps are treated as missing.
#' @param reference A \code{DNAbin} alignment, a named character vector of
#'   equal-length aligned sequences, or \code{NULL} (default) to compare against
#'   the bundled MalAvi alignment for \code{version}.
#' @param version MalAvi release to use when \code{reference} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @param top_n Return only the nearest \code{top_n} references. \code{NULL}
#'   (default) returns all references, ordered nearest first.
#' @return A \code{data.frame} with one row per reference (or per nearest
#'   \code{top_n}), with columns \code{lineage}, \code{distance} (mismatches over
#'   comparable positions), and \code{n_comparable} (positions actually compared),
#'   ordered nearest first.
#' @seealso \code{\link{lineage_qc}}, \code{\link{ambiguous_pairs}},
#'   \code{\link{blast_malavi}}
#' @examples
#' aln <- extract_alignment()
#' q <- paste(as.character(aln[1, ]), collapse = "")
#' ## nearest MalAvi lineages to this sequence (it is itself, so distance 0 first)
#' head(pairwise_deletion_distance(q, top_n = 5))
#' @export
pairwise_deletion_distance <- function(query, reference = NULL,
                                       version = "latest", top_n = NULL) {
  if (!is.character(query) || length(query) != 1L) {
    stop("`query` must be a single aligned sequence (a length-1 character vector).",
         call. = FALSE)
  }
  query <- .qc_clean_seq(query)

  ## coerce the reference to a coded matrix (rows = named sequences); the helper
  ## accepts a DNAbin alignment, a named character vector, or NULL (bundled data)
  charmat   <- .qc_char_matrix(reference, version)
  ref_len   <- ncol(charmat)
  ref_names <- rownames(charmat)

  ## the query must already be aligned to the same frame as the reference
  if (nchar(query) != ref_len) {
    stop("`query` is ", nchar(query), " bp but the reference alignment is ",
         ref_len, " bp; align the query to the reference frame first ",
         "(see frame_to_malavi() or blast_malavi()).", call. = FALSE)
  }

  refcode <- .qc_code_matrix(charmat)
  qcode   <- .qc_code_vec(strsplit(query, "", fixed = TRUE)[[1]])

  n  <- nrow(refcode)
  tn <- if (is.null(top_n)) n else min(as.integer(top_n), n)

  ## .qc_nearest does the vectorized pairwise-deletion comparison and the
  ## distance / -n_comparable ordering; drop its internal row-index column
  res <- .qc_nearest(qcode, refcode, ref_names, top_n = tn)
  res$index <- NULL
  rownames(res) <- NULL
  res
}
