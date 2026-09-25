#' Follow old lineage names forward through the Lineage History
#'
#' MalAvi lineage names are cited in hundreds of papers. Since 2026-09-25 every
#' change to what a name means -- a sequence completed under the old name, a
#' rename, a merge into another lineage, a removal -- is recorded in the
#' \code{"Lineage History"} table (see \code{\link{extract_table}}). This
#' function walks that table so an analysis made on an older release can be
#' mapped onto the current one.
#'
#' @param names Character vector of lineage names as used in the older analysis.
#' @param version MalAvi release to read, as \code{\link{extract_table}} takes it.
#' @return A \code{data.frame} with one row per input name: \code{name} (as
#'   given), \code{current} (the name in force in the release, \code{NA} when
#'   the lineage was removed, or the input name unchanged when the history has
#'   no event for it), \code{status} (one of \code{"unchanged"},
#'   \code{"renamed"}, \code{"merged"}, \code{"removed"},
#'   \code{"sequence_updated"}) and \code{trail} (the events walked, as one
#'   string, empty when unchanged). A name that is neither in the history nor in
#'   the release's Grand Lineage Summary is reported \code{status = "unknown"}.
#' @examples
#' lineage_history_map(c("SGS1", "STRALU01"))
#' @export
lineage_history_map <- function(names, version = "latest") {
  history <- extract_table("Lineage History", version = version)
  gls     <- extract_table("Grand Lineage Summary", version = version)
  known   <- unique(gls$LINEAGE_NAME)

  ## forward links: old name -> (event, new name), for the events that move a name
  moving  <- history[history$EVENT %in% c("renamed", "merged_into", "removed"), , drop = FALSE]
  updated <- unique(history$LINEAGE_NAME[history$EVENT == "sequence_updated"])

  one <- function(name) {
    here  <- name
    trail <- character(0)
    seen  <- character(0)
    final_event <- NA_character_
    while (here %in% moving$OLD_NAME && !here %in% seen) {
      seen <- c(seen, here)
      row  <- moving[match(here, moving$OLD_NAME), ]
      trail <- c(trail, paste(here, row$EVENT, row$NEW_NAME))
      final_event <- row$EVENT
      if (row$EVENT == "removed") {
        return(list(current = NA_character_, status = "removed",
                    trail = paste(trail, collapse = "; ")))
      }
      here <- row$NEW_NAME
    }
    status <- if (!is.na(final_event)) {
      if (final_event == "merged_into") "merged" else "renamed"
    } else if (here %in% updated) {
      "sequence_updated"
    } else if (here %in% known) {
      "unchanged"
    } else {
      "unknown"
    }
    list(current = if (status == "unknown") NA_character_ else here,
         status = status, trail = paste(trail, collapse = "; "))
  }

  rows <- lapply(as.character(names), one)
  out <- data.frame(
    name    = as.character(names),
    current = vapply(rows, function(r) r$current, character(1)),
    status  = vapply(rows, function(r) r$status, character(1)),
    trail   = vapply(rows, function(r) r$trail, character(1)),
    stringsAsFactors = FALSE)
  .malavi_attach_meta(out, malavi_version = .malavi_resolve_version(version))
}
