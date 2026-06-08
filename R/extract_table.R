#' Get a MalAvi data table
#'
#' Returns one of the MalAvi data tables from the database bundled in
#' the package. MalAvi is no longer downloaded from the web; the tables come from
#' the release shipped with \code{malaviR} (see \code{\link{malavi_version}}).
#'
#' The bundled release provides five tables:
#' \describe{
#'   \item{\code{"Hosts and Sites Table"}}{host records, sites, and references (\code{hosts_and_sites}).}
#'   \item{\code{"Grand Lineage Summary"}}{per-lineage summary, including the sequence (\code{grand_lineage_summary}).}
#'   \item{\code{"Morpho Species Summary"}}{lineages linked to morphologically described species (\code{morpho_species}).}
#'   \item{\code{"Table of References"}}{reference list (\code{references}).}
#'   \item{\code{"Vector Data Table"}}{vector records (\code{vector_data}).}
#' }
#' Either the descriptive name above or its short \code{snake_case} key may be
#' supplied.
#'
#' @param table Name of the table to return (see Details), or \code{"all"} to
#'   return a named list of all five tables. Defaults to \code{"Hosts and Sites Table"}.
#' @param version MalAvi release to read, as a date string (e.g.
#'   \code{"2026-03-23"}) or \code{"latest"} (default).
#' @return A \code{data.frame}, or for \code{table = "all"} a named list of
#'   \code{data.frame}s.
#' @seealso \code{\link{extract_alignment}}, \code{\link{malavi_version}}
#' @examples
#' hosts <- extract_table("Hosts and Sites Table")
#' head(hosts)
#' @export
extract_table <- function(table = "Hosts and Sites Table", version = "latest") {
  if (length(table) != 1L) {
    stop('\'table\' must be a single table name (or "all").', call. = FALSE)
  }
  ## map descriptive names (and the snake_case keys) to the bundled list elements
  lookup <- c(
    "Hosts and Sites Table"  = "hosts_and_sites",
    "Grand Lineage Summary"  = "grand_lineage_summary",
    "Morpho Species Summary" = "morpho_species",
    "Table of References"    = "references",
    "Vector Data Table"      = "vector_data",
    "hosts_and_sites"        = "hosts_and_sites",
    "grand_lineage_summary"  = "grand_lineage_summary",
    "morpho_species"         = "morpho_species",
    "references"             = "references",
    "vector_data"            = "vector_data"
  )

  if (!(table %in% c(names(lookup), "all"))) {
    stop('Please choose one of the following tables:\n',
         '"Hosts and Sites Table", "Grand Lineage Summary", "Morpho Species Summary",\n',
         '"Table of References", "Vector Data Table", or "all".', call. = FALSE)
  }

  db <- .malavi_load(version, "malavi_db_")

  if (table == "all") {
    return(db[c("hosts_and_sites", "grand_lineage_summary", "morpho_species",
                "references", "vector_data")])
  }
  db[[lookup[[table]]]]
}
