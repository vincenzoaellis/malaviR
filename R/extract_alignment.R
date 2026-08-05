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
#' The alignment is a \strong{matrix} \code{DNAbin} (sequences in rows, alignment
#' columns in columns), not a list of sequences. So use \code{nrow()} and
#' \code{rownames()} to count and name the sequences: \code{length()} returns the
#' number of matrix cells (sequences times sites) and \code{names()} returns
#' \code{NULL}. Subset sequences with \code{aln[i, ]} and sites with
#' \code{aln[, j]}.
#'
#' @param version MalAvi release to read, as a date string (e.g.
#'   \code{"2026-03-23"}) or \code{"latest"} (default).
#' @param genus Parasite genus/genera to return. Either \code{"all"} (default,
#'   the whole alignment) or one or more of \code{"Plasmodium"},
#'   \code{"Haemoproteus"}, \code{"Leucocytozoon"}, and \code{"other"}.
#' @return A matrix \code{DNAbin} alignment, optionally subset by genus.
#' @seealso \code{\link{extract_table}}, \code{\link{clean_alignment}}
#' @examples
#' aln <- extract_alignment()
#' dim(aln)                  # sequences x sites
#' nrow(aln)                 # number of sequences (not length(aln))
#' head(rownames(aln))       # sequence names (not names(aln))
#' plas <- extract_alignment(genus = "Plasmodium")
#' @export
extract_alignment <- function(version = "latest",
                              genus = c("all", "Plasmodium", "Haemoproteus",
                                        "Leucocytozoon", "other")) {
  genus <- match.arg(genus, several.ok = TRUE)
  alignment <- .malavi_load(version, "malavi_db_")$alignment
  mv <- .malavi_resolve_version(version)

  if ("all" %in% genus) {
    return(.malavi_attach_meta(alignment, malavi_version = mv, genus = "all"))
  }

  ## map each lineage's name prefix to a genus
  prefix <- substr(rownames(alignment), 1, 2)
  lineage_genus <- c("P_" = "Plasmodium", "H_" = "Haemoproteus",
                     "L_" = "Leucocytozoon")[prefix]
  lineage_genus[is.na(lineage_genus)] <- "other"

  keep <- lineage_genus %in% genus
  if (!any(keep)) stop("No lineages match genus: ", paste(genus, collapse = ", "),
                       call. = FALSE)
  .malavi_attach_meta(alignment[keep, ], malavi_version = mv,
                      genus = paste(genus, collapse = ","))
}
