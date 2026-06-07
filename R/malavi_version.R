#' MalAvi database version bundled in the package
#'
#' Returns the version (release date) of the MalAvi database snapshot that
#' \code{malaviR} reads from. MalAvi is no longer queried online, so the
#' "version" is simply the date stamp of the bundled release (e.g.
#' \code{"2026-03-23"}). Use \code{which = "all"} to list every release bundled in
#' your installation.
#'
#' @param which Either \code{"latest"} (default) to return the most recent
#'   bundled release, or \code{"all"} to return all bundled releases.
#' @return A character vector of version (date) string(s).
#' @seealso \code{\link{extract_table}}, \code{\link{extract_alignment}}
#' @examples
#' malavi_version()
#' @export
malavi_version <- function(which = c("latest", "all")) {
  which <- match.arg(which)
  vers <- .malavi_versions("malavi_db_")
  if (length(vers) == 0) {
    stop("No bundled MalAvi data found in the package. Reinstall malaviR.", call. = FALSE)
  }
  if (which == "latest") vers[1] else vers
}
