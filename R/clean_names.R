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
#' @examples
#' clean_names(c("H_COLL2_Haemoproteus_pallidus", "P_GRW04_Plasmodium_relictum", "L_CIAE02"))
#' clean_names(c("H_COLL2_Haemoproteus_pallidus", "L_CIAE02"), keep.genus = TRUE)
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
  if(keep.genus == FALSE){
    return(lineage)
  } else{
    out.n <- data.frame(parasiteGenus = prefix, Lineage_Name = lineage,
                        stringsAsFactors = FALSE) %>%
      mutate(parasiteGenus = recode(parasiteGenus,
                                    P = "Plasmodium",
                                    H = "Haemoproteus",
                                    L = "Leucocytozoon"))
    return(out.n)
  }
}
