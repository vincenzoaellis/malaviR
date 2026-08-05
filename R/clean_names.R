#' Clean MalAvi lineage names to match the tables
#'
#' MalAvi alignment tip labels carry a parasite-genus prefix (e.g.
#' \code{"H_COLL2"}), and often a trailing morphological-species name as well
#' (e.g. \code{"H_COLL2_Haemoproteus_pallidus"}), whereas the data tables store
#' the lineage name alone (e.g. \code{"COLL2"}). This helper strips the prefix
#' and any trailing morphological-species name so names from an alignment can be
#' matched to the tables, and can optionally return the parasite genus alongside
#' the cleaned name.
#'
#' @param lin.names Character vector of lineage names of the form
#'   \code{"<genus prefix>_<lineage>"}, optionally followed by a
#'   morphological-species name (e.g. from \code{rownames()} of an alignment).
#' @param keep.genus If \code{FALSE} (default), return just the cleaned lineage
#'   names as a character vector. If \code{TRUE}, return a \code{data.frame} with
#'   the parasite genus (\code{P}/\code{H}/\code{L} expanded to
#'   \emph{Plasmodium}/\emph{Haemoproteus}/\emph{Leucocytozoon}) and the cleaned
#'   \code{Lineage_Name}.
#' @return A character vector, or a \code{data.frame} when \code{keep.genus = TRUE}.
#'
#'   Only the \code{P}/\code{H}/\code{L} prefixes are expanded to genus names; any
#'   other prefix is returned as-is, because MalAvi does use others (lineages whose
#'   \code{GENUS_NAME} is recorded as \code{"N/A"} carry an \code{N_} prefix, which
#'   is a missing genus rather than a genus called \code{N} -- see
#'   \code{\link{malavi_issues}}).
#'
#'   A name with no underscore has no genus prefix to read, so with
#'   \code{keep.genus = TRUE} its genus is \code{NA} and a warning is raised. This
#'   matters because such names do occur in practice -- a bare GenBank accession
#'   mixed in with alignment labels, for instance -- and returning the accession
#'   itself as the parasite genus would be a plausible-looking wrong value. The
#'   lineage name is still returned unchanged, since a name with no prefix is
#'   already in the form the tables use.
#' @examples
#' clean_names(c("H_COLL2_Haemoproteus_pallidus", "P_GRW04_Plasmodium_relictum", "L_CIAE02"))
#' clean_names(c("H_COLL2_Haemoproteus_pallidus", "L_CIAE02"), keep.genus = TRUE)
#'
#' ## a name with no genus prefix gets an NA genus, not a made-up one
#' clean_names(c("H_COLL2", "GU085191"), keep.genus = TRUE)
#' @importFrom dplyr mutate
#' @importFrom dplyr recode
#' @importFrom magrittr %>%
#' @export

clean_names <- function(lin.names, keep.genus = FALSE){
  parasiteGenus = NULL
  lin.names <- as.character(lin.names)
  ## genus prefix = text before the first "_"; lineage = next "_"-delimited token
  prefix  <- sub("^([^_]+)_.*$", "\\1", lin.names)
  lineage <- sub("^[^_]+_([^_]+).*$", "\\1", lin.names)

  ## A name with no underscore carries no genus prefix, and both sub() calls
  ## above leave such a name untouched -- so the whole name would be reported as
  ## the parasite genus. That is a wrong value that looks plausible (a bare
  ## GenBank accession would come back as a genus), so record the genus as
  ## missing instead and say so. The lineage name is left alone: a name with no
  ## prefix is already in the form the data tables use.
  has_prefix <- grepl("_", lin.names, fixed = TRUE)
  prefix[!has_prefix] <- NA_character_

  if(keep.genus == FALSE){
    return(lineage)
  } else{
    if(any(!has_prefix)){
      warning(sum(!has_prefix), " name(s) have no '_' and so no parasite-genus ",
              "prefix; their genus is returned as NA (e.g. '",
              lin.names[!has_prefix][1], "').", call. = FALSE)
    }
    out.n <- data.frame(parasiteGenus = prefix, Lineage_Name = lineage,
                        stringsAsFactors = FALSE) %>%
      mutate(parasiteGenus = recode(parasiteGenus,
                                    P = "Plasmodium",
                                    H = "Haemoproteus",
                                    L = "Leucocytozoon"))
    return(out.n)
  }
}
