#' Get a MalAvi data table
#'
#' Returns one of the MalAvi data tables from the database bundled in
#' the package. MalAvi is no longer downloaded from the web; the tables come from
#' the release shipped with \code{malaviR} (see \code{\link{malavi_version}}).
#'
#' The bundled release provides seven tables:
#' \describe{
#'   \item{\code{"Hosts and Sites Table"}}{host records, sites, and references (\code{hosts_and_sites}).}
#'   \item{\code{"Grand Lineage Summary"}}{per-lineage summary, including the sequence (\code{grand_lineage_summary}).}
#'   \item{\code{"Morpho Species Summary"}}{lineages linked to morphologically described species (\code{morpho_species}).}
#'   \item{\code{"Table of References"}}{reference list (\code{references}).}
#'   \item{\code{"Vector Data Table"}}{vector records (\code{vector_data}).}
#'   \item{\code{"Lineage History"}}{what changed about a lineage name and when:
#'     sequence updates, renames, merges and removals (\code{lineage_history}).
#'     Added to the release format on 2026-09-25; for a bundled release that
#'     predates it the table is returned empty, with its ten columns, so code
#'     written against it runs on every release. See
#'     \code{\link{lineage_history_map}} to follow old names forward.}
#'   \item{\code{"Genomic Resources"}}{genomes, transcriptomes and other large
#'     sequence records linked to a cyt b lineage: what kind, where deposited, how
#'     long, from which host and study (\code{genomic_resources}). Added to the
#'     release format on 2026-09-25 and first filled on 2026-10-14; for a bundled
#'     release that predates it the table is returned empty, with its eleven
#'     columns.}
#' }
#' Either the descriptive name above or its short \code{snake_case} key may be
#' supplied.
#'
#' Note: in the bundled 2026-03-23 release the source spreadsheets encoded some
#' blank cells in a way that made them read back as the literal string
#' \code{"LINEAGE_NAME"} (mostly in \code{grand_lineage_summary}, and a few in
#' \code{vector_data}). These have been reset to \code{NA} when building the
#' bundled data; all other values are verbatim from MalAvi.
#'
#' The returned tables also have their character columns whitespace-tidied: any
#' run of whitespace (including stray line breaks and tabs that the source
#' spreadsheets embed in free-text fields such as \code{SPECIES_NAME} or
#' \code{COMMENT}, and non-breaking spaces, which some \code{GENBANK_ACC} values
#' begin with) is collapsed to a single space and the ends are trimmed. This
#' only repairs formatting -- it never changes a value's identity -- but it keeps
#' e.g. host-name joins and distinct-species counts from being broken by an
#' invisible trailing newline. \code{NA} values are preserved.
#'
#' @param table Name of the table to return (see Details), or \code{"all"} to
#'   return a named list of all seven tables. Defaults to \code{"Hosts and Sites Table"}.
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
    "Lineage History"        = "lineage_history",
    "Genomic Resources"      = "genomic_resources",
    "hosts_and_sites"        = "hosts_and_sites",
    "grand_lineage_summary"  = "grand_lineage_summary",
    "morpho_species"         = "morpho_species",
    "references"             = "references",
    "vector_data"            = "vector_data",
    "lineage_history"        = "lineage_history",
    "genomic_resources"      = "genomic_resources"
  )

  if (!(table %in% c(names(lookup), "all"))) {
    stop('Please choose one of the following tables:\n',
         '"Hosts and Sites Table", "Grand Lineage Summary", "Morpho Species Summary",\n',
         '"Table of References", "Vector Data Table", "Lineage History",\n',
         '"Genomic Resources", or "all".',
         call. = FALSE)
  }

  db <- .malavi_load(version, "malavi_db_")
  mv <- .malavi_resolve_version(version)

  ## A release bundled before 2026-09-25 has no lineage_history element. It is
  ## returned as an empty table with the release format's columns rather than
  ## NULL or an error, so a caller can read it on any release; the attribute
  ## "malavi_absent" says the table was not in the bundle.
  if (is.null(db$lineage_history)) {
    db$lineage_history <- .malavi_empty_lineage_history()
  }
  ## The same for the Genomic Resources table (release format since 2026-09-25).
  if (is.null(db$genomic_resources)) {
    db$genomic_resources <- .malavi_empty_genomic_resources()
  }

  if (table == "all") {
    tabs <- db[c("hosts_and_sites", "grand_lineage_summary", "morpho_species",
                 "references", "vector_data", "lineage_history", "genomic_resources")]
    ## Tidy here too. Until version 1.1.2 the "all" branch returned the raw
    ## bundled tables while every single-table call tidied, so the same table
    ## differed in 3,143 cells depending on how it was asked for -- and the
    ## documentation described the tidy without qualification.
    tabs <- lapply(tabs, .clean_table_ws)
    ## stamp the bundled version onto each table (and the list) as a provenance
    ## attribute; the descriptive table label is deliberately not stamped so that
    ## e.g. extract_table("references") == extract_table("Table of References")
    return(.malavi_attach_meta(lapply(tabs, .malavi_attach_meta,
                                      malavi_version = mv),
                               malavi_version = mv))
  }
  .malavi_attach_meta(.clean_table_ws(db[[lookup[[table]]]]), malavi_version = mv)
}


## The Lineage History table's columns, as the release format declares them
## (malavi_rebuild: release_store.TABLES["lineage_history"]). Used to return an
## empty, correctly shaped table for a bundle that predates the table.
.malavi_lineage_history_columns <- c(
  "LINEAGE_NAME", "EVENT", "OLD_NAME", "NEW_NAME", "OLD_SEQUENCE",
  "NEW_SEQUENCE", "DATE", "REASON", "SOURCE", "RELEASE")

.malavi_empty_lineage_history <- function() {
  out <- as.data.frame(setNames(replicate(length(.malavi_lineage_history_columns),
                                          character(0), simplify = FALSE),
                                .malavi_lineage_history_columns),
                       stringsAsFactors = FALSE)
  attr(out, "malavi_absent") <- TRUE
  out
}


## The Genomic Resources table's columns, as the release format declares them
## (malavi_rebuild: release_store.TABLES["genomic_resources"]).
.malavi_genomic_resources_columns <- c(
  "LINEAGE_NAME", "RESOURCE_TYPE", "ACCESSION", "DATABASE", "LENGTH_BP",
  "DESCRIPTION", "HOST_SPECIES", "REFERENCE_NAME", "DOI", "SOURCE", "NOTES")

.malavi_empty_genomic_resources <- function() {
  out <- as.data.frame(setNames(replicate(length(.malavi_genomic_resources_columns),
                                          character(0), simplify = FALSE),
                                .malavi_genomic_resources_columns),
                       stringsAsFactors = FALSE)
  attr(out, "malavi_absent") <- TRUE
  out
}
