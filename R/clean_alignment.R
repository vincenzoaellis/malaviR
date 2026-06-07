#' Identify and collapse repeated haplotypes in a MalAvi alignment
#'
#' Different MalAvi lineage names are sometimes assigned to the same sequence
#' ("synonymies"), which inflates estimates of parasite diversity, especially
#' when short, partial sequences are involved (Tamayo-Quintero et al. 2025). This
#' function finds groups of lineages that share a haplotype, returns a table of
#' those synonymies, and produces a de-duplicated alignment that keeps one
#' lineage per group.
#'
#' By default this function is deterministic: the most complete sequence in each
#' group is kept (ties broken alphabetically). Set \code{select = "random"} for
#' the quick, non-deterministic behaviour of earlier \code{malaviR} versions,
#' which keep one lineage per group at random (call \code{\link{set.seed}} first
#' for reproducibility). In either case, supply \code{keep} to override the choice
#' for specific groups; any group without a supplied choice falls back to the
#' \code{select} rule.
#'
#' Two definitions of "same haplotype" are available via \code{method}:
#' \describe{
#'   \item{\code{"strict"} (default)}{sequences identical across the whole
#'     alignment, including gaps -- the behaviour of the original function.}
#'   \item{\code{"overlap"}}{additionally collapses a partial sequence into any
#'     strictly more complete sequence that is identical to it over the partial's
#'     informative (non-gap/non-\code{N}) positions. This catches the
#'     partial-sequence synonymies highlighted by Tamayo-Quintero et al. (2025),
#'     but is slower on large alignments.}
#' }
#' The \code{informative_length} column (count of A/C/G/T bases) helps flag the
#' short, partial sequences at the heart of the problem.
#'
#' @param alignment A \code{DNAbin} alignment (e.g. from
#'   \code{\link{extract_alignment}}).
#' @param method How to define a repeated haplotype: \code{"strict"} (default) or
#'   \code{"overlap"} (see Details).
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
#' aln <- extract_alignment()
#' res <- clean_alignment(aln)
#' head(res$synonymies)
#'
#' ## quick random pick (reproducible with a seed)
#' set.seed(1)
#' res_rand <- clean_alignment(aln, select = "random")
#' @export
clean_alignment <- function(alignment, method = c("strict", "overlap"),
                            select = c("complete", "random"), keep = NULL) {
  method <- match.arg(method)
  select <- match.arg(select)
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

  list(synonymies = syn$synonymies, kept = syn$kept, dropped = dropped,
       alignment_clean = alignment_clean)
}
