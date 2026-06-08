#' Version of the clootl taxonomy bundled in the package
#'
#' \code{\link{match_taxonomy}} matches host names against a snapshot of the
#' \pkg{clootl} (eBird/Clements) avian taxonomy that is bundled with
#' \code{malaviR}. This returns the taxonomy year of that snapshot.
#'
#' @return The bundled clootl taxonomy year (an integer).
#' @references
#' McTavish EJ, Gerbracht JA, Holder MT, Iliff MJ, Lepage D, Rasmussen PC,
#' Redelings BD, Sanchez Reyes LL, Miller ET (2025). A complete and dynamic tree
#' of birds. Proceedings of the National Academy of Sciences 122(18):
#' e2409658122. \doi{10.1073/pnas.2409658122}
#' @seealso \code{\link{match_taxonomy}}
#' @examples
#' clootl_taxonomy_version()
#' @export
clootl_taxonomy_version <- function() {
  clootl_year
}
