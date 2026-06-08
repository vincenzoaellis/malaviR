#' Get the MalAvi sequence alignment
#'
#' Returns the aligned MalAvi cytochrome \emph{b} sequences from the database
#' bundled in the package, as a \code{DNAbin} object. MalAvi is no longer
#' downloaded from the web; the alignment comes from the release shipped with
#' \code{malaviR} (see \code{\link{malavi_version}}).
#'
#' Lineage names are prefixed by parasite genus: \code{P_} (\emph{Plasmodium}),
#' \code{H_} (\emph{Haemoproteus}), \code{L_} (\emph{Leucocytozoon}); any other
#' prefix is treated as \code{"other"}. Use \code{genus} to subset the alignment
#' to one or more genera. Note that some tip labels also carry a morphological
#' species name appended after the lineage code (e.g.
#' \code{"H_COLL2_Haemoproteus_pallidus"}).
#'
#' @param version MalAvi release to read, as a date string (e.g.
#'   \code{"2026-03-23"}) or \code{"latest"} (default).
#' @param genus Parasite genus/genera to return. Either \code{"all"} (default,
#'   the whole alignment) or one or more of \code{"Plasmodium"},
#'   \code{"Haemoproteus"}, \code{"Leucocytozoon"}, and \code{"other"}.
#' @return A \code{DNAbin} alignment, optionally subset by genus.
#' @seealso \code{\link{extract_table}}, \code{\link{clean_alignment}}
#' @examples
#' aln <- extract_alignment()
#' dim(aln)
#' plas <- extract_alignment(genus = "Plasmodium")
#' @export
extract_alignment <- function(version = "latest",
                              genus = c("all", "Plasmodium", "Haemoproteus",
                                        "Leucocytozoon", "other")) {
  genus <- match.arg(genus, several.ok = TRUE)
  alignment <- .malavi_load(version, "malavi_db_")$alignment

  if ("all" %in% genus) return(alignment)

  ## map each lineage's name prefix to a genus
  prefix <- substr(rownames(alignment), 1, 2)
  lineage_genus <- c("P_" = "Plasmodium", "H_" = "Haemoproteus",
                     "L_" = "Leucocytozoon")[prefix]
  lineage_genus[is.na(lineage_genus)] <- "other"

  keep <- lineage_genus %in% genus
  if (!any(keep)) stop("No lineages match genus: ", paste(genus, collapse = ", "),
                       call. = FALSE)
  alignment[keep, ]
}
