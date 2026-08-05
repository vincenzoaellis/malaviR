#' Identify and collapse repeated haplotypes in a MalAvi alignment
#'
#' Shorter MalAvi lineages (i.e., < 479 bp) sometimes match perfectly to longer
#' sequences that have different lineage names ("synonymies"), and it has been
#' pointed out in the literature that this inflates estimates of parasite
#' diversity (Tamayo-Quintero et al. 2025). This function finds groups of
#' lineages that share a haplotype, returns a table of those synonymies, and
#' produces a de-duplicated alignment that keeps one lineage per group.
#'
#' By default this function is deterministic: the most complete (i.e., longest)
#' sequence in each group is kept (ties broken alphabetically). Set
#' \code{select = "random"} for the quick random selection of the earlier
#' \code{malaviR} version, which keeps one lineage per group at random (call
#' \code{\link{set.seed}} first for reproducibility). In either case, supply
#' \code{keep} to override the choice for specific groups (i.e., if you want to
#' choose particular lineages to represent a haplotype group); any group without
#' a supplied choice falls back to the \code{select} rule.
#'
#' Two definitions of "same haplotype" are available via \code{method}:
#' \describe{
#'   \item{\code{"overlap"} (default)}{collapses a partial sequence into any
#'     strictly more complete sequence that is identical to it over the partial's
#'     informative (non-gap/non-\code{N}) positions, in addition to collapsing
#'     fully identical sequences. This catches the partial-sequence synonymies
#'     highlighted by Tamayo-Quintero et al. (2025), but is slower on large
#'     alignments.}
#'   \item{\code{"strict"}}{collapses only sequences that are identical across the
#'     whole alignment, including gaps -- the behavior of the original function
#'     from the earlier \code{malaviR} version.}
#' }
#' The \code{informative_length} column (count of A/C/G/T bases) helps flag the
#' short, partial sequences at the heart of the problem.
#'
#' Note on \code{method = "overlap"}: a partial sequence is collapsed into a more
#' complete one when they are \emph{compatible} at every position where both carry
#' a definite A/C/G/T base. Positions missing in the partial sequence are treated
#' as \strong{unknown}, not as evidence of identity, so an overlap group is a group
#' of sequences \emph{compatible over their observed sites} rather than sequences
#' proven identical across the unobserved ones. Keeping the most complete member
#' discards no observed base, but it does assume the shorter sequence would have
#' matched at the positions it never determined.
#'
#' @param alignment A \code{DNAbin} alignment, or \code{NULL} (default) to use the
#'   bundled MalAvi alignment for \code{version} -- the same default as
#'   \code{\link{synonymy_report}}.
#' @param version MalAvi release to use when \code{alignment} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @param method How to define a repeated haplotype: \code{"overlap"} (default)
#'   or \code{"strict"} (see Details).
#' @param select How to pick the lineage kept from each synonymy group when it is
#'   not named in \code{keep}: \code{"complete"} (default) keeps the most complete
#'   sequence (ties broken alphabetically); \code{"random"} keeps one at random
#'   (set a seed first for reproducibility).
#' @param keep Optional character vector of lineage names to keep. For each
#'   synonymy group containing one of these names, that name is kept; an error is
#'   raised if a single group contains more than one supplied name.
#' @return A list with elements:
#'   \describe{
#'     \item{\code{synonymies}}{a \code{data.frame}, one row per lineage in a
#'       repeated-haplotype group, with columns \code{haplotype} (group id),
#'       \code{lineage}, \code{informative_length}, and \code{status}
#'       (\code{"kept"} or \code{"dropped"}).}
#'     \item{\code{kept}}{character vector of lineages kept.}
#'     \item{\code{dropped}}{character vector of lineages dropped.}
#'     \item{\code{alignment_clean}}{the \code{DNAbin} alignment with dropped
#'       lineages removed.}
#'   }
#' @references
#' Tamayo-Quintero J, Martinez-de la Puente J, Matta NE, Pacheco MA,
#' Rivera-Gutierrez HF (2025). Imprudent use of MalAvi names biases the
#' estimation of parasite diversity of avian haemosporidians. PLoS Pathogens
#' 21(2): e1012911. \doi{10.1371/journal.ppat.1012911}
#' @seealso \code{\link{synonymy_report}}, \code{\link{extract_alignment}}
#' @examples
#' ## defaults to the bundled alignment
#' res <- clean_alignment()
#' head(res$synonymies)
#'
#' ## or pass your own
#' aln <- extract_alignment(genus = "Plasmodium")
#' res_plas <- clean_alignment(aln)
#'
#' ## quick random pick (reproducible with a seed)
#' set.seed(1)
#' res_rand <- clean_alignment(aln, select = "random")
#' @export
clean_alignment <- function(alignment = NULL, method = c("overlap", "strict"),
                            select = c("complete", "random"), keep = NULL,
                            version = "latest") {
  method <- match.arg(method)
  select <- match.arg(select)
  ## record whether the bundled alignment was used *before* reassigning
  ## `alignment`, so the version stamp below is correct
  used_bundled <- is.null(alignment)
  if (is.null(alignment)) alignment <- extract_alignment(version = version)
  if (!inherits(alignment, "DNAbin")) {
    stop("The alignment should be of class 'DNAbin'.", call. = FALSE)
  }

  g <- .haplotype_groups(alignment, method)
  if (!any(table(g$group) > 1)) {
    stop("The alignment has no repeated haplotypes.", call. = FALSE)
  }

  syn <- .build_synonymies(g$lineages, g$group, g$informative_length,
                           select = select, keep = keep)
  dropped <- syn$synonymies$lineage[syn$synonymies$status == "dropped"]
  alignment_clean <- alignment[!g$lineages %in% dropped, ]

  out <- list(synonymies = syn$synonymies, kept = syn$kept, dropped = dropped,
              alignment_clean = alignment_clean)
  ## MalAvi version is only meaningful when the bundled alignment was cleaned
  mv <- if (used_bundled) .malavi_resolve_version(version) else NA_character_
  out <- .malavi_attach_meta(out, malavi_version = mv, method = method,
                             select = select)
  class(out) <- c("malavi_alignment_clean", class(out))
  out
}

#' @export
print.malavi_alignment_clean <- function(x, ...) {
  cat("MalAvi cleaned alignment\n")
  meta <- .malavi_meta_line(x)
  if (!is.null(meta)) cat("  ", meta, "\n", sep = "")
  cat("  synonymy groups collapsed: ", length(x$kept), "\n", sep = "")
  cat("  lineages dropped:          ", length(x$dropped), "\n", sep = "")
  cat("  sequences in $alignment_clean: ", nrow(x$alignment_clean), "\n", sep = "")
  if (identical(attr(x, "malavi_meta")$method, "overlap")) {
    cat("\nNote: overlap groups are compatible over observed A/C/G/T sites;\n",
        "one representative per group is kept in $alignment_clean.\n", sep = "")
  }
  cat("\nSee $synonymies for the per-lineage kept/dropped table.\n")
  invisible(x)
}
