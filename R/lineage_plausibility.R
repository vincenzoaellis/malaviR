## Host and biogeographic plausibility of a lineage detection. Where lineage_qc()
## asks "is this sequence plausible?", this asks "is this lineage plausible in
## this host, in this place?" -- a different kind of evidence, drawn from the
## MalAvi host records rather than from the sequence.
##
## The whole function is a set of lookups into the Hosts and Sites table plus the
## per-region columns of the Grand Lineage Summary. Nothing is modeled; the
## output is a description of what has been recorded before.

## The per-region presence columns of the Grand Lineage Summary. Listed
## explicitly because the table also holds non-region columns (PASSERIFORMES,
## SUM_*), and intersected with the real column names at use so that a release
## that drops or renames one degrades rather than crashes.
.malavi_region_columns <- c(
  "EUROPE", "SUB_SAHARAN_AFRICA", "NORTH_AFRICA_AND_MIDDLE_EAST",
  "NORTH_AMERICA", "HAWAI", "CENTRAL_AMERICA", "SOUTH_AMERICA", "ASIA",
  "AUSTRALIA_AND_NEW_ZEALAND", "OCEANIA", "ANTARCTICA", "UNKNOWN_REGION"
)

## Normalize a user-supplied region to a Grand Lineage Summary column name:
## upper case, non-alphanumerics to underscore, plus the handful of spellings
## people reliably get "wrong" because MalAvi's own spelling is unusual.
.malavi_normalize_region <- function(region) {
  norm <- toupper(trimws(as.character(region)))
  norm <- gsub("[^A-Z0-9]+", "_", norm)
  norm <- gsub("^_|_$", "", norm)
  alias <- c(HAWAII = "HAWAI",
             MIDDLE_EAST = "NORTH_AFRICA_AND_MIDDLE_EAST",
             NORTH_AFRICA = "NORTH_AFRICA_AND_MIDDLE_EAST",
             AFRICA = "SUB_SAHARAN_AFRICA",
             AUSTRALIA = "AUSTRALIA_AND_NEW_ZEALAND",
             NEW_ZEALAND = "AUSTRALIA_AND_NEW_ZEALAND")
  hit <- match(norm, names(alias))
  norm[!is.na(hit)] <- unname(alias[hit[!is.na(hit)]])
  norm
}

#' Is a lineage plausible in this host, in this place?
#'
#' \strong{Experimental.} For each detection of a MalAvi lineage, report what has
#' been recorded before: how well known the lineage is, whether it has been found
#' in that host species (and in its genus, family and order), and whether it has
#' been found in that country or region. This is the host-and-biogeography
#' counterpart to \code{\link{lineage_qc}}, which screens the sequence itself.
#'
#' The intended use is deep-sequenced amplicon data, where low-abundance
#' detections are a mix of genuine low-level or migration-acquired infections and
#' laboratory artifacts (index hopping, cross-contamination, carry-over). Those
#' two look identical in read counts but not in their ecology: a genuine
#' low-level infection is usually a lineage that is already known from that
#' region, often from a related host, whereas a contaminant is typically a
#' lineage that is abundant elsewhere in the same run but has no business in that
#' host or that place. This function supplies that context; you supply the read
#' counts and the judgment.
#'
#' @section What "novel" does and does not mean:
#' Every column here describes \strong{sampling}, not biology. MalAvi records
#' where people have looked, and host-parasite records are heavily biased toward
#' well-studied host families and well-studied countries. A new host species for
#' a common, generalist lineage is unremarkable; a new host \emph{order} for a
#' lineage known from one study on another continent is worth a hard look. Read
#' novelty as "no prior record", never as "impossible", and never treat a flag on
#' its own as grounds to discard a detection.
#'
#' The evidence is drawn entirely from MalAvi: host taxonomy comes from the Hosts
#' and Sites table (so a host species that has never been sampled for
#' haemosporidians at all is reported as \code{host_not_in_malavi} rather than
#' being looked up elsewhere), and regions come from the per-region columns of
#' the Grand Lineage Summary. Note also that host records accumulate against a
#' lineage name, so the synonymy problem propagates here: a lineage whose
#' partial-sequence synonyms carry their own host records will look narrower than
#' it is (see \code{\link{synonymy_report}}).
#'
#' @param lineage Character vector of MalAvi lineage names. Both bare names
#'   (\code{"SGS1"}) and full alignment names (\code{"P_SGS1"},
#'   \code{"H_COLL2_Haemoproteus_pallidus"}) are accepted; they are passed
#'   through \code{\link{clean_names}}.
#' @param host Optional character vector of host species as binomials
#'   (e.g. \code{"Parus major"}), matching MalAvi's \code{SPECIES_NAME}. Length 1
#'   or the same length as \code{lineage}.
#' @param country Optional character vector of countries, matching MalAvi's
#'   \code{COUNTRY_NAME} (e.g. \code{"Sweden"}). Length 1 or the same length as
#'   \code{lineage}.
#' @param region Optional character vector of broad regions, matching the
#'   per-region columns of the Grand Lineage Summary. Case and punctuation are
#'   normalized, so \code{"North America"}, \code{"north_america"} and
#'   \code{"NORTH_AMERICA"} are the same; a few common alternatives
#'   (\code{"Hawaii"}, \code{"Africa"}, \code{"Australia"}) are mapped to
#'   MalAvi's spelling. Unrecognized regions raise an error listing the valid
#'   ones.
#' @param version MalAvi release to use; a date string or \code{"latest"}
#'   (default).
#' @return A \code{data.frame} of class \code{malavi_plausibility}, one row per
#'   detection, with
#'   \describe{
#'     \item{\code{lineage}}{the cleaned lineage name.}
#'     \item{\code{call}}{a one-word summary: \code{lineage_not_in_malavi},
#'       \code{previously_recorded}, \code{new_host}, \code{new_host_family},
#'       \code{new_location}, or \code{new_host_and_location}.}
#'     \item{\code{flags}}{semicolon-separated detail behind the call
#'       (see \emph{Flags}).}
#'     \item{\code{n_studies}, \code{n_host_records}, \code{n_countries}}{how well
#'       known the lineage is, from \code{\link{lineage_studies}}.}
#'     \item{\code{host}, \code{host_family}, \code{host_order}}{the host as
#'       supplied, and its family and order as recorded in MalAvi (\code{NA} if
#'       the host does not appear in MalAvi at all).}
#'     \item{\code{host_recorded}, \code{host_genus_recorded},
#'       \code{host_family_recorded}, \code{host_order_recorded}}{whether this
#'       lineage has been recorded in that host species, genus, family, order.}
#'     \item{\code{n_records_host}}{host records of this lineage in that species.}
#'     \item{\code{country}, \code{country_recorded}, \code{n_records_country}}{the
#'       country as supplied, whether this lineage has been recorded there, and
#'       how many records.}
#'     \item{\code{region}, \code{region_recorded}}{the normalized region, and
#'       whether the Grand Lineage Summary marks the lineage present in it.}
#'   }
#'   Columns for an argument you did not supply are \code{NA}.
#'
#' @section Flags:
#' \describe{
#'   \item{\code{lineage_not_in_malavi}}{no host record for this lineage at all.}
#'   \item{\code{single_study_lineage}}{reported by exactly one study, so its
#'     known host and geographic range is one study's sampling.}
#'   \item{\code{host_not_in_malavi}}{the host species appears nowhere in the
#'     Hosts and Sites table, so there is nothing to compare against -- it has
#'     probably just never been screened.}
#'   \item{\code{new_host_species}, \code{new_host_genus},
#'     \code{new_host_family}, \code{new_host_order}}{no prior record of this
#'     lineage at that taxonomic level.}
#'   \item{\code{new_country}, \code{new_region}}{no prior record of this lineage
#'     from that country or region.}
#' }
#'
#' @seealso \code{\link{lineage_qc}} (screens the sequence),
#'   \code{\link{lineage_studies}}, \code{\link{extract_table}},
#'   \code{\link{malavi_issues}}
#' @examples
#' ## a well-known generalist in a well-known host and place
#' lineage_plausibility("SGS1", host = "Parus major", country = "Sweden")
#'
#' ## the same lineage somewhere it has not been recorded
#' lineage_plausibility("SGS1", host = "Parus major", region = "Antarctica")
#'
#' ## several detections at once, as you would get from an amplicon run
#' lineage_plausibility(
#'   lineage = c("SGS1", "GRW04", "TUPHI01"),
#'   host    = "Parus major",
#'   country = "Sweden"
#' )
#' @export
lineage_plausibility <- function(lineage, host = NULL, country = NULL,
                                 region = NULL, version = "latest") {
  lineage <- clean_names(as.character(lineage))
  n <- length(lineage)
  if (n == 0) stop("`lineage` is empty.", call. = FALSE)

  ## Recycle the optional arguments to the number of detections. Only length 1
  ## and length n are allowed: silent partial recycling would quietly pair the
  ## wrong host with the wrong lineage, which is the one mistake this function
  ## must not make.
  recycle <- function(x, nm) {
    if (is.null(x)) return(rep(NA_character_, n))
    x <- as.character(x)
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) != n) {
      stop("`", nm, "` must be length 1 or the same length as `lineage` (",
           n, "), not ", length(x), ".", call. = FALSE)
    }
    x
  }
  host    <- recycle(host, "host")
  country <- recycle(country, "country")
  region  <- recycle(region, "region")

  hosts <- extract_table("Hosts and Sites Table", version = version)
  gls   <- extract_table("Grand Lineage Summary", version = version)
  st    <- lineage_studies(version = version)

  ## ---- what MalAvi has recorded for each lineage -------------------------
  ## Only the queried lineages matter, so cut the ~18k host records down to
  ## their rows before grouping. Indexing all ~5,000 lineages when the user
  ## asked about three is the difference between a second and ten.
  wanted   <- unique(stats::na.omit(lineage))
  rel_rows <- which(hosts$LINEAGE_NAME %in% wanted)
  rel      <- hosts[rel_rows, , drop = FALSE]

  ## Group the relevant records by lineage, then reduce each group to the sets
  ## of hosts, host taxa and countries it covers. Doing this per detection
  ## would rescan the records every time.
  by_lineage <- split(seq_len(nrow(rel)), rel$LINEAGE_NAME)

  set_of <- function(column) {
    lapply(by_lineage, function(i) unique(stats::na.omit(rel[[column]][i])))
  }
  host_species_of <- set_of("SPECIES_NAME")
  host_genus_of   <- set_of("GENUS_NAME")
  host_family_of  <- set_of("FAMILY_NAME")
  host_order_of   <- set_of("ORDER_NAME")
  country_of      <- set_of("COUNTRY_NAME")

  ## ---- host taxonomy, as MalAvi knows it ---------------------------------
  ## Family and order for a host species, taken from any record of that species
  ## anywhere in MalAvi (not just records of this lineage).
  host_key    <- hosts$SPECIES_NAME
  first_hit   <- match(host, host_key)
  host_family <- hosts$FAMILY_NAME[first_hit]
  host_order  <- hosts$ORDER_NAME[first_hit]
  ## Host genus is the first word of the binomial; that is how MalAvi's
  ## GENUS_NAME relates to its SPECIES_NAME, and it still works for a host that
  ## has never been sampled.
  host_genus  <- ifelse(is.na(host), NA_character_,
                        sub("^([^ ]+).*$", "\\1", trimws(host)))
  host_in_malavi <- !is.na(host) & !is.na(first_hit)

  ## ---- region presence, from the Grand Lineage Summary --------------------
  region_norm <- .malavi_normalize_region(region)
  region_cols <- intersect(.malavi_region_columns, names(gls))
  supplied    <- unique(stats::na.omit(region_norm))
  unknown     <- setdiff(supplied, region_cols)
  if (length(unknown) > 0) {
    stop("Unknown region: ", paste(unknown, collapse = ", "),
         ".\nValid regions: ", paste(region_cols, collapse = ", "),
         call. = FALSE)
  }
  ## the region columns are marked "1" where the lineage occurs and blank
  ## (NA) where it does not
  gls_lineage_row <- match(lineage, gls$LINEAGE_NAME)

  ## ---- assemble one row per detection -------------------------------------
  study_row <- match(lineage, st$lineage)
  in_malavi <- !is.na(study_row)

  ## logical helper: is `value` in the recorded set for this lineage?
  recorded_in <- function(sets, value) {
    vapply(seq_len(n), function(i) {
      if (is.na(value[i]) || is.na(lineage[i])) return(NA)
      s <- sets[[lineage[i]]]
      if (is.null(s)) return(FALSE)
      value[i] %in% s
    }, logical(1))
  }

  host_recorded        <- recorded_in(host_species_of, host)
  host_genus_recorded  <- recorded_in(host_genus_of,   host_genus)
  host_family_recorded <- recorded_in(host_family_of,  host_family)
  host_order_recorded  <- recorded_in(host_order_of,   host_order)
  country_recorded     <- recorded_in(country_of,      country)

  ## counts of records backing the host and country hits, over the relevant
  ## rows only (same reason as the grouping above)
  count_records <- function(column, value) {
    vapply(seq_len(n), function(i) {
      if (is.na(value[i])) return(NA_integer_)
      rows <- by_lineage[[lineage[i]]]
      if (is.null(rows)) return(0L)
      col <- rel[[column]][rows]
      sum(!is.na(col) & col == value[i])
    }, integer(1))
  }
  n_records_host    <- count_records("SPECIES_NAME", host)
  n_records_country <- count_records("COUNTRY_NAME", country)

  region_recorded <- vapply(seq_len(n), function(i) {
    if (is.na(region_norm[i]) || is.na(gls_lineage_row[i])) return(NA)
    val <- gls[[region_norm[i]]][gls_lineage_row[i]]
    !is.na(val) && nzchar(trimws(as.character(val)))
  }, logical(1))

  ## ---- flags and the one-word call ---------------------------------------
  flags <- vector("list", n)
  call  <- character(n)
  for (i in seq_len(n)) {
    f <- character(0)

    if (!is.na(host[i]) && !host_in_malavi[i]) f <- c(f, "host_not_in_malavi")

    if (!in_malavi[i]) {
      ## A lineage with no host record at all is trivially new in every host
      ## and every place, so the novelty flags would say nothing. Report only
      ## that it is unknown. (The logical columns still record FALSE, which is
      ## factually what the lookup found.)
      f <- c(f, "lineage_not_in_malavi")
    } else {
      if (isTRUE(st$n_studies[study_row[i]] == 1)) f <- c(f, "single_study_lineage")
      if (isFALSE(host_recorded[i]))        f <- c(f, "new_host_species")
      if (isFALSE(host_genus_recorded[i]))  f <- c(f, "new_host_genus")
      if (isFALSE(host_family_recorded[i])) f <- c(f, "new_host_family")
      if (isFALSE(host_order_recorded[i]))  f <- c(f, "new_host_order")
      if (isFALSE(country_recorded[i]))     f <- c(f, "new_country")
      if (isFALSE(region_recorded[i]))      f <- c(f, "new_region")
    }

    ## The call collapses the flags into the single distinction a user acts on:
    ## is the host new, is the place new, are both.
    new_host_any <- isFALSE(host_recorded[i])
    new_place    <- isFALSE(country_recorded[i]) || isFALSE(region_recorded[i])
    call[i] <-
      if (!in_malavi[i])                          "lineage_not_in_malavi"
      else if (isFALSE(host_family_recorded[i]))  "new_host_family"
      else if (new_host_any && new_place)         "new_host_and_location"
      else if (new_host_any)                      "new_host"
      else if (new_place)                         "new_location"
      else                                        "previously_recorded"

    flags[[i]] <- f
  }

  out <- data.frame(
    lineage              = lineage,
    call                 = call,
    flags                = vapply(flags, paste, character(1), collapse = "; "),
    n_studies            = st$n_studies[study_row],
    n_host_records       = st$n_host_records[study_row],
    n_countries          = st$n_countries[study_row],
    host                 = host,
    host_family          = host_family,
    host_order           = host_order,
    host_recorded        = host_recorded,
    host_genus_recorded  = host_genus_recorded,
    host_family_recorded = host_family_recorded,
    host_order_recorded  = host_order_recorded,
    n_records_host       = n_records_host,
    country              = country,
    country_recorded     = country_recorded,
    n_records_country    = n_records_country,
    region               = region_norm,
    region_recorded      = region_recorded,
    stringsAsFactors     = FALSE
  )
  rownames(out) <- NULL

  out <- .malavi_attach_meta(out, malavi_version = .malavi_resolve_version(version))
  class(out) <- c("malavi_plausibility", class(out))
  out
}

#' @export
print.malavi_plausibility <- function(x, ...) {
  cat("MalAvi host and biogeographic plausibility\n")
  meta <- .malavi_meta_line(x)
  if (!is.null(meta)) cat("  ", meta, "\n", sep = "")
  cat("  Novelty here means no prior record, not implausibility: MalAvi\n")
  cat("  records where people have looked. Weigh it with read abundance.\n\n")

  ## a compact view; the full table is still there as a data frame
  cols <- c("lineage", "call", "n_studies", "host", "host_recorded",
            "country", "country_recorded", "region", "region_recorded")
  print(as.data.frame(x)[, intersect(cols, names(x)), drop = FALSE])

  if (any(nzchar(x$flags))) {
    cat("\nFlags:\n")
    for (i in which(nzchar(x$flags))) {
      cat("  ", x$lineage[i], ": ", x$flags[i], "\n", sep = "")
    }
  }
  invisible(x)
}
