## Detect "ambiguous pairs" of MalAvi reference lineages: two lineages that agree
## wherever both are determined but are separated only by ambiguity/gap sites,
## with neither cleanly contained in the other. See internal.R for the shared,
## dependency-light sequence-coding helpers reused here.

#' Find ambiguous reference pairs in a MalAvi alignment
#'
#' \strong{Experimental.} Reports pairs of MalAvi reference lineages that are an
#' \emph{ambiguous pair}: under pairwise deletion they have \strong{distance 0}
#' (they disagree at no position where \emph{both} carry a definite A/C/G/T base),
#' yet each one is determined at \strong{at least one} position where the other is
#' an ambiguity code or gap. Because each lineage carries a base the other lacks,
#' neither is cleanly contained in the other, and we cannot tell from the
#' reference sequences alone whether the two names denote the same haplotype or
#' two genuinely one-base-different haplotypes -- the deciding site is undetermined
#' in at least one of them.
#'
#' This is the companion to, and deliberately distinct from,
#' \code{\link{synonymy_report}}:
#' \describe{
#'   \item{synonymy (\code{synonymy_report})}{one sequence is \emph{contained} in
#'     another (identical wherever the shorter one is determined), so the more
#'     complete sequence can represent both without discarding any observed base.}
#'   \item{ambiguous pair (this function)}{the two are \emph{mutually partial} --
#'     each is determined where the other is not -- so neither contains the other
#'     and there is no single sequence that can stand in for both without
#'     discarding an observed base. These are exactly the pairs
#'     \code{synonymy_report} does \strong{not} flag.}
#' }
#' The function identifies these pairs.
#'
#' A true one-base \emph{neighbor} (a real base-versus-base difference at a site
#' determined in both) has distance \eqn{\ge 1} and is therefore never reported
#' here; that is a genuine distinct lineage, to be resolved by other evidence
#' (e.g. within-sample read abundance), not by ambiguity.
#'
#' Only pairs in which \emph{both} members carry ambiguity can qualify: if either
#' sequence is fully determined it cannot have a position for the other to be
#' privately determined at, so any distance-0 relationship to it is a containment
#' (synonymy), not an ambiguous pair. The search is therefore restricted to the
#' partial (ambiguity-bearing) sequences, which is both correct and fast.
#'
#' @param alignment A \code{DNAbin} alignment, a named character vector of
#'   equal-length aligned sequences, or \code{NULL} (default) to use the bundled
#'   MalAvi alignment for \code{version}.
#' @param version MalAvi release to use when \code{alignment} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @param min_comparable Minimum number of positions determined in \emph{both}
#'   members for a pair to be reported (default \code{1}, which only excludes the
#'   degenerate case of two sequences that never overlap on a determined base).
#'   A distance of 0 over very few comparable positions is weak evidence of
#'   ambiguity; raise this to focus on well-overlapping pairs.
#' @return A list with:
#'   \describe{
#'     \item{\code{summary}}{a one-row \code{data.frame}: \code{n_sequences},
#'       \code{n_partial} (sequences carrying any ambiguity/gap),
#'       \code{n_ambiguous_pairs}, \code{n_lineages_in_pairs} (distinct lineages
#'       appearing in any pair), \code{n_same_genus}, \code{n_cross_genus}, and the
#'       \code{min_comparable} used.}
#'     \item{\code{by_genus}}{counts of pairs split by the genus combination.}
#'     \item{\code{pairs}}{a \code{data.frame}, one row per ambiguous pair, with
#'       \code{lineage_a}, \code{lineage_b}, \code{genus_a}, \code{genus_b},
#'       \code{same_genus}, \code{n_comparable} (positions determined in both, all
#'       agreeing), \code{a_private} / \code{b_private} (positions determined in
#'       one member but ambiguous in the other), ordered strongest first (most
#'       comparable positions, then fewest private positions).}
#'   }
#' @seealso \code{\link{synonymy_report}}, \code{\link{lineage_qc}},
#'   \code{\link{clean_alignment}}
#' @examples
#' ap <- ambiguous_pairs()
#' ap$summary
#' head(ap$pairs)
#' @export
ambiguous_pairs <- function(alignment = NULL, version = "latest",
                            min_comparable = 1L) {
  ## resolve the input to an upper-case character matrix (rows = lineages with
  ## names, cols = alignment positions); reuses the shared QC coercion so DNAbin,
  ## named character vectors, and the bundled alignment are all accepted
  charmat <- .qc_char_matrix(alignment, version)
  if (is.null(rownames(charmat))) {
    stop("The alignment must have lineage names (row names).", call. = FALSE)
  }
  min_comparable <- as.integer(min_comparable)

  ## integer-code: A/C/G/T -> 1:4, everything else (N, other IUPAC codes, gaps)
  ## -> 0, i.e. "not a determined base"
  code <- .qc_code_matrix(charmat)
  ref_names <- rownames(code)
  n_sequences <- nrow(code)
  aln_length  <- ncol(code)

  ## a position is "determined" where it carries a definite base; a sequence is
  ## "partial" if it has any undetermined position. Only partial sequences can be
  ## members of an ambiguous pair (see Details), so restrict to them.
  determined <- code > 0L
  is_partial <- rowSums(determined) < aln_length
  partial_index <- which(is_partial)
  n_partial <- length(partial_index)

  empty_pairs <- data.frame(
    lineage_a = character(0), lineage_b = character(0),
    genus_a = character(0), genus_b = character(0), same_genus = logical(0),
    n_comparable = integer(0), a_private = integer(0), b_private = integer(0),
    stringsAsFactors = FALSE)

  ## fewer than two partial sequences => no pair is possible
  if (n_partial < 2L) {
    return(.ambiguous_pairs_result(empty_pairs, n_sequences, n_partial,
                                   min_comparable))
  }

  ## ---- vectorized pairwise computation, partial sequences only ----
  ## Cp: coded bases; Dp: determined indicator (1/0); Ap: ambiguous indicator.
  ## All matrix products below are over the n_partial x n_partial pair space.
  Cp <- code[partial_index, , drop = FALSE]
  Dp <- (Cp > 0L) + 0                      # numeric so BLAS handles the products
  Ap <- 1 - Dp

  ## comparable[i, j] = positions determined in BOTH i and j
  comparable <- tcrossprod(Dp)

  ## agree[i, j] = positions where BOTH are determined and carry the SAME base,
  ## summed over the four bases (one indicator matrix per base)
  agree <- matrix(0, n_partial, n_partial)
  for (base_code in 1:4) {
    onehot <- (Cp == base_code) + 0
    agree <- agree + tcrossprod(onehot)
  }

  ## pairwise-deletion distance over comparable positions; 0 means "no conflict
  ## anywhere both are determined"
  mismatch <- comparable - agree

  ## a_private[i, j] = positions determined in i but ambiguous in j; the matching
  ## "b_private" (determined in j, ambiguous in i) is simply its transpose
  a_private_mat <- tcrossprod(Dp, Ap)

  ## an ambiguous pair: distance 0, enough overlap, and BOTH sides private
  ## (mutually partial => neither contained). Take the upper triangle so each
  ## unordered pair is reported once.
  keep <- (mismatch == 0) & (comparable >= min_comparable) &
    (a_private_mat > 0) & (t(a_private_mat) > 0) & upper.tri(mismatch)
  hits <- which(keep)

  ## free the large intermediates before assembling the (small) result
  b_private_mat <- t(a_private_mat)
  rm(agree, mismatch, Dp, Ap, Cp); gc(verbose = FALSE)

  if (length(hits) == 0L) {
    return(.ambiguous_pairs_result(empty_pairs, n_sequences, n_partial,
                                   min_comparable))
  }

  rc <- arrayInd(hits, c(n_partial, n_partial))   # row/col index of each hit
  a_idx <- partial_index[rc[, 1]]                 # back to alignment row numbers
  b_idx <- partial_index[rc[, 2]]

  pairs <- data.frame(
    lineage_a    = ref_names[a_idx],
    lineage_b    = ref_names[b_idx],
    genus_a      = .ambiguous_genus(ref_names[a_idx]),
    genus_b      = .ambiguous_genus(ref_names[b_idx]),
    n_comparable = as.integer(comparable[hits]),
    a_private    = as.integer(a_private_mat[hits]),
    b_private    = as.integer(b_private_mat[hits]),
    stringsAsFactors = FALSE)
  pairs$same_genus <- pairs$genus_a == pairs$genus_b
  pairs <- pairs[c("lineage_a", "lineage_b", "genus_a", "genus_b", "same_genus",
                   "n_comparable", "a_private", "b_private")]

  ## strongest first: most comparable positions, then fewest private positions
  ord <- order(-pairs$n_comparable, pmax(pairs$a_private, pairs$b_private))
  pairs <- pairs[ord, , drop = FALSE]
  rownames(pairs) <- NULL

  .ambiguous_pairs_result(pairs, n_sequences, n_partial, min_comparable)
}

## Assemble the summary/by_genus/pairs list from a (possibly empty) pairs frame.
.ambiguous_pairs_result <- function(pairs, n_sequences, n_partial,
                                    min_comparable) {
  lineages_in_pairs <- unique(c(pairs$lineage_a, pairs$lineage_b))

  summary <- data.frame(
    n_sequences          = n_sequences,
    n_partial            = n_partial,
    n_ambiguous_pairs    = nrow(pairs),
    n_lineages_in_pairs  = length(lineages_in_pairs),
    n_same_genus         = sum(pairs$same_genus),
    n_cross_genus        = sum(!pairs$same_genus),
    min_comparable       = min_comparable,
    stringsAsFactors = FALSE)

  ## label each pair by its (sorted) genus combination, e.g. "Haemoproteus"
  ## within-genus or "Haemoproteus|Plasmodium" cross-genus
  combo <- if (nrow(pairs)) {
    apply(cbind(pairs$genus_a, pairs$genus_b), 1L,
          function(g) paste(sort(g), collapse = "|"))
  } else character(0)
  by_genus <- as.data.frame(table(genus_combination = combo),
                            responseName = "n_pairs", stringsAsFactors = FALSE)

  list(summary = summary, by_genus = by_genus, pairs = pairs)
}

## Map a MalAvi lineage name to its parasite genus from the leading prefix, the
## same convention used by synonymy_report(); anything else is "other".
.ambiguous_genus <- function(x) {
  g <- c("P_" = "Plasmodium", "H_" = "Haemoproteus",
         "L_" = "Leucocytozoon")[substr(x, 1, 2)]
  g[is.na(g)] <- "other"
  unname(g)
}
