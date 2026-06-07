#' Version of the clootl taxonomy bundled in the package
#'
#' \code{\link{match_taxonomy}} matches host names against a snapshot of the
#' \pkg{clootl} (eBird/Clements) avian taxonomy that is bundled with
#' \code{malaviR}. This returns the taxonomy year of that snapshot. The snapshot
#' is bundled (rather than fetched at run time) because older clootl taxonomy
#' years and phylogenies can be difficult to access later; matching to the
#' bundled current taxonomy keeps results reproducible.
#'
#' @return The bundled clootl taxonomy year (an integer).
#' @seealso \code{\link{match_taxonomy}}
#' @examples
#' clootl_taxonomy_version()
#' @export
clootl_taxonomy_version <- function() {
  clootl_year
}
