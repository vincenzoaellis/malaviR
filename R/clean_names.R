#' Clean MalAvi lineage names to match the database
#'
#' MalAvi alignment tip labels carry a parasite-genus prefix (e.g.
#' \code{"H_COLL2"}), whereas the data tables store the bare lineage name (e.g.
#' \code{"COLL2"}). This helper strips the prefix so names from an alignment can
#' be matched to the tables, and can optionally return the parasite genus
#' alongside the cleaned name.
#'
#' @param lin.names Character vector of lineage names of the form
#'   \code{"<genus prefix>_<lineage>"} (e.g. from \code{rownames()} of an
#'   alignment).
#' @param keep.genus If \code{FALSE} (default), return just the cleaned lineage
#'   names as a character vector. If \code{TRUE}, return a \code{data.frame} with
#'   the parasite genus (\code{P}/\code{H}/\code{L} expanded to
#'   \emph{Plasmodium}/\emph{Haemoproteus}/\emph{Leucocytozoon}) and the cleaned
#'   \code{Lineage_Name}.
#' @return A character vector, or a \code{data.frame} when \code{keep.genus = TRUE}.
#' @examples
#' clean_names(c("H_COLL2", "P_GRW4", "L_CIAE02"))
#' clean_names(c("H_COLL2", "P_GRW4"), keep.genus = TRUE)
#' @importFrom dplyr mutate
#' @importFrom dplyr recode
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export

clean_names <- function(lin.names, keep.genus = FALSE){
  X1=X2=parasiteGenus=NULL
  out <- data.frame(do.call("rbind", strsplit(as.character(lin.names), "_")))
  if(keep.genus == FALSE){
    return(as.character(out[,2]))
  } else{
    out.n <- out %>%
      select(X1, X2) %>%
      rename(parasiteGenus = X1, Lineage_Name = X2) %>%
      mutate(parasiteGenus = recode(parasiteGenus,
                                    P = "Plasmodium",
                                    H = "Haemoproteus",
                                    L = "Leucocytozoon"))
    return(out.n)
  }
}
