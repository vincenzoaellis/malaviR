#' List the MalAvi database releases bundled in the package
#'
#' \code{malaviR} ships one or more MalAvi database snapshots inside the package.
#' Each is identified by its release date. This function lists the versions
#' available in your installation, newest first.
#'
#' @return A character vector of version (date) strings, newest first. The first
#'   element is what functions use when \code{version = "latest"}.
#' @seealso \code{\link{malavi_version}}, \code{\link{extract_table}},
#'   \code{\link{extract_alignment}}
#' @examples
#' malavi_versions()
#' @export
malavi_versions <- function() {
  vers <- .malavi_versions("malavi_db_")
  if (length(vers) == 0) {
    stop("No bundled MalAvi data found in the package. Reinstall malaviR.", call. = FALSE)
  }
  vers
}
