#' Quantify MalAvi haplotype synonymies for investigation
#'
#' Summarizes how many lineage names share a haplotype with another name and
#' returns the lineage names in groups so they can be examined. By default it
#' reports on the bundled MalAvi alignment using the overlap-aware definition of
#' a haplotype (which catches short, partial sequences identical to a longer
#' one), but any alignment and either method may be used.
#'
#' This is a reporting companion to \code{\link{clean_alignment}}: use this to see
#' the size of the problem and which lineages to check, and \code{clean_alignment}
#' to actually produce a de-duplicated alignment.
#'
#' @param alignment A \code{DNAbin} alignment. If \code{NULL} (default), the
#'   bundled MalAvi alignment for \code{version} is used.
#' @param method How to define a repeated haplotype: \code{"overlap"} (default)
#'   or \code{"strict"}. See \code{\link{clean_alignment}}.
#' @param version MalAvi release to use when \code{alignment} is \code{NULL}; a
#'   date string or \code{"latest"} (default).
#' @return A list with:
#'   \describe{
#'     \item{\code{summary}}{a one-row \code{data.frame} of counts:
#'       \code{n_sequences}, \code{n_haplotypes} (distinct haplotypes),
#'       \code{n_synonymous_haplotypes} (haplotypes carrying >1 lineage name),
#'       \code{n_lineages_in_synonymies}, \code{n_redundant_names}
#'       (\code{n_sequences - n_haplotypes}, the diversity inflation),
#'       \code{pct_diversity_inflation}, and \code{n_partial_sequences}.}
#'     \item{\code{by_genus}}{redundant-name counts split by parasite genus.}
#'     \item{\code{synonymies}}{a \code{data.frame} of the synonymy groups, one
#'       row per lineage, with \code{haplotype}, \code{lineage}, \code{genus},
#'       \code{informative_length}, \code{is_partial}, and \code{status} -- the
#'       list of names to investigate.}
#'   }
#' @references
#' Tamayo-Quintero J, Martinez-de la Puente J, Matta NE, Pacheco MA,
#' Rivera-Gutierrez HF (2025). Imprudent use of MalAvi names biases the
#' estimation of parasite diversity of avian haemosporidians. PLoS Pathogens
#' 21(2): e1012911. \doi{10.1371/journal.ppat.1012911}
#' @seealso \code{\link{clean_alignment}}, \code{\link{extract_alignment}}
#' @examples
#' rep <- synonymy_report(method = "strict")
#' rep$summary
#' head(rep$synonymies)
#' @export
synonymy_report <- function(alignment = NULL, method = c("overlap", "strict"),
                            version = "latest") {
  method <- match.arg(method)
  if (is.null(alignment)) alignment <- extract_alignment(version = version)
  if (!inherits(alignment, "DNAbin")) {
    stop("The alignment should be of class 'DNAbin'.", call. = FALSE)
  }

  g <- .haplotype_groups(alignment, method)
  syn <- .build_synonymies(g$lineages, g$group, g$informative_length)
  synonymies <- syn$synonymies

  ## label each synonymy lineage with parasite genus and partial-sequence flag
  genus_of <- function(x) {
    p <- c("P_" = "Plasmodium", "H_" = "Haemoproteus",
           "L_" = "Leucocytozoon")[substr(x, 1, 2)]
    p[is.na(p)] <- "other"
    unname(p)
  }
  ## (works on a zero-row frame too, so a synonymy-free alignment is handled)
  synonymies$genus <- genus_of(synonymies$lineage)
  synonymies$is_partial <- synonymies$informative_length < g$aln_length
  synonymies <- synonymies[c("haplotype", "lineage", "genus",
                             "informative_length", "is_partial", "status")]

  n_sequences  <- length(g$lineages)
  n_haplotypes <- length(unique(g$group))
  n_redundant  <- n_sequences - n_haplotypes

  summary <- data.frame(
    n_sequences              = n_sequences,
    n_haplotypes             = n_haplotypes,
    n_synonymous_haplotypes  = length(unique(synonymies$haplotype)),
    n_lineages_in_synonymies = nrow(synonymies),
    n_redundant_names        = n_redundant,
    pct_diversity_inflation  = round(100 * n_redundant / n_haplotypes, 2),
    n_partial_sequences      = sum(g$informative_length < g$aln_length)
  )

  dropped <- synonymies[synonymies$status == "dropped", , drop = FALSE]
  by_genus <- as.data.frame(table(genus = dropped$genus),
                            responseName = "n_redundant_names")

  list(summary = summary, by_genus = by_genus, synonymies = synonymies)
}
