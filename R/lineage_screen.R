## Tools for identifying possible sequencing error or otherwise incorrect sequences. For example, a lineage reported by only one study that also
## carries singleton amino-acid changes (i.e., non-synonymous changes found in no other sequences) might be an error.
## lineage_studies() counts the studies per lineage (sort of a cheap helper); lineage_screen() joins those
## counts to a per-lineage mutation summary.

#' Count how many studies report each MalAvi lineage
#'
#' For every lineage, count the number of distinct studies (unique
#' \code{REFERENCE_NAME}) that have reported it in the Hosts and Sites table,
#' along with the number of host records and the number of distinct countries
#' (\code{COUNTRY_NAME}) it has been found in. Pair this with
#' \code{\link{lineage_screen}} to compare those counts with the number of
#' singleton non-synonymous mutations across the database.
#'
#' @param version Bundled MalAvi release to use; a date string or \code{"latest"}
#'   (default).
#' @param references If \code{TRUE}, add a \code{references} column that spells
#'   out \emph{which} studies reported each lineage (the distinct
#'   \code{REFERENCE_NAME} values, semicolon-separated), rather than just the
#'   count in \code{n_studies}. This is a convenience for tracing a lineage back
#'   to its source papers; it is off by default because the long strings make the
#'   table harder to read.
#' @return A \code{data.frame} ordered by lineage with columns
#'   \describe{
#'     \item{\code{lineage}}{the MalAvi lineage name (e.g. \code{"SGS1"}).}
#'     \item{\code{parasite_genus}}{the most frequent parasite genus recorded for
#'       the lineage.}
#'     \item{\code{n_studies}}{number of distinct studies reporting the lineage.}
#'     \item{\code{n_host_records}}{number of host records (rows) for the lineage.}
#'     \item{\code{n_countries}}{number of distinct countries the lineage has been
#'       recorded in.}
#'     \item{\code{references}}{(only if \code{references = TRUE}) the distinct
#'       study names, semicolon-separated.}
#'   }
#' @seealso \code{\link{lineage_screen}}, \code{\link{extract_table}}
#' @examples
#' studies <- lineage_studies()
#' head(studies[order(-studies$n_studies), ])     # most-reported lineages
#' studies[studies$lineage %in% c("SGS1", "GRW11", "GRW04"), ]
#' @export
lineage_studies <- function(version = "latest", references = FALSE) {
  hosts <- extract_table("Hosts and Sites Table", version = version)

  lineage <- trimws(hosts$LINEAGE_NAME)
  ref     <- trimws(hosts$REFERENCE_NAME)
  genus   <- trimws(hosts$PARASITE_GENUS)
  country <- trimws(hosts$COUNTRY_NAME)
  ref[is.na(ref) | ref == "" | ref == "N/A"] <- NA_character_
  country[is.na(country) | country == "" | country == "N/A"] <- NA_character_

  ## keep only rows with a real lineage name
  ok <- !is.na(lineage) & lineage != ""
  lineage <- lineage[ok]; ref <- ref[ok]; genus <- genus[ok]; country <- country[ok]

  ## group the row indices by lineage and summarize each group
  idx <- split(seq_along(lineage), lineage)
  n_studies   <- vapply(idx, function(i) length(unique(stats::na.omit(ref[i]))), integer(1))
  n_records   <- vapply(idx, length, integer(1))
  n_countries <- vapply(idx, function(i) length(unique(stats::na.omit(country[i]))), integer(1))
  pgenus      <- vapply(idx, function(i) .modal(genus[i]), character(1))

  out <- data.frame(lineage = names(idx), parasite_genus = pgenus,
                    n_studies = n_studies, n_host_records = n_records,
                    n_countries = n_countries,
                    row.names = NULL, stringsAsFactors = FALSE)
  if (references) {
    out$references <- vapply(idx, function(i)
      paste(sort(unique(stats::na.omit(ref[i]))), collapse = "; "), character(1))
  }
  out[order(out$lineage), , drop = FALSE]
}

#' Screen every MalAvi lineage for study support versus mutation burden
#'
#' \strong{Experimental.} Builds, in one pass over the bundled alignment, a
#' per-lineage table that pairs how many studies report a lineage
#' (\code{\link{lineage_studies}}) with how many \emph{singleton} substitutions the
#' lineage carries in its cytochrome \emph{b} sequence. A singleton substitution is
#' a base the lineage \strong{alone} carries at a well-covered alignment site: the
#' signature of a sequencing error, since a real variant is shared by more than
#' one lineage. Each singleton substitution is classified, against the consensus
#' codon (frame 1, genetic code 4) with the singleton base swapped in, as
#' synonymous, non-synonymous, or stop-creating.
#'
#' This makes the database-wide version of Staffan Bensch's observation
#' straightforward: lineages reported by a single study that also carry several
#' singleton non-synonymous changes are the strongest sequencing-error candidates.
#' For example, filter to \code{in_hosts_table & n_studies == 1} and compare
#' \code{n_singleton_nonsynonymous} against well-replicated lineages.
#'
#' @param version MalAvi release to use; a date string or \code{"latest"}
#'   (default).
#' @param reference Reference alignment to screen: a \code{DNAbin} alignment, a
#'   named character vector, or \code{NULL} (default) to use the bundled MalAvi
#'   alignment for \code{version}. Lineage names and genera are parsed from the
#'   sequence names, so this is meant for the MalAvi alignment naming scheme.
#' @param genus Restrict the screen to one parasite genus: \code{"all"} (default,
#'   pool every lineage), \code{"Plasmodium"}, \code{"Haemoproteus"}, or
#'   \code{"Leucocytozoon"}. Because the three genera are deeply divergent, a
#'   singleton (a base carried by one lineage alone) is judged relative to the
#'   chosen group's diversity: restricting to a genus sharpens the screen by not
#'   letting between-genus differences dominate the per-site consensus.
#' @param min_site_coverage Fraction of sequences that must carry an unambiguous
#'   base at a site for it to count toward singleton substitutions (default 0.5).
#'   This keeps sparsely sequenced alignment columns from masquerading as singleton
#'   substitutions.
#' @param studies If \code{TRUE} (default), join study counts from
#'   \code{\link{lineage_studies}} by lineage name.
#' @return A \code{data.frame}, one row per sequence in the alignment, with
#'   \describe{
#'     \item{\code{lineage}}{the MalAvi lineage name parsed from the sequence name.}
#'     \item{\code{parasite_genus}}{genus from the sequence-name prefix
#'       (Haemoproteus / Leucocytozoon / Plasmodium).}
#'     \item{\code{n_studies}, \code{n_host_records}, \code{n_countries}}{study,
#'       host-record, and distinct-country counts (only if \code{studies = TRUE});
#'       \code{NA} for lineages absent from the Hosts and Sites table.}
#'     \item{\code{in_hosts_table}}{whether the lineage has any host record (only
#'       if \code{studies = TRUE}).}
#'     \item{\code{seq_length}}{number of unambiguous bases in the sequence (a
#'       short value flags a partial sequence).}
#'     \item{\code{n_singleton_substitutions}}{bases the lineage alone carries at a
#'       well-covered site.}
#'     \item{\code{n_singleton_nonsynonymous}}{singleton substitutions that change the
#'       amino acid (stop-creating ones are included here).}
#'     \item{\code{n_singleton_synonymous}}{singleton substitutions that do not change
#'       the amino acid.}
#'     \item{\code{n_singleton_stop}}{singleton substitutions that create a stop codon
#'       (a subset of \code{n_singleton_nonsynonymous}).}
#'     \item{\code{seq_name}}{the original alignment sequence name.}
#'   }
#' @seealso \code{\link{lineage_studies}}, \code{\link{lineage_qc}}
#' @examples
#' \donttest{
#' library(dplyr)
#'
#' ## Staffan-style comparison: do single-study lineages carry more singleton
#' ## non-synonymous changes than well-replicated ones?
#' lineage_screen() %>%
#'   filter(in_hosts_table) %>%
#'   group_by(single_study = n_studies == 1) %>%
#'   summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#'
#' ## sharper within a genus (singletons judged against Plasmodium alone)
#' lineage_screen(genus = "Plasmodium") %>%
#'   filter(in_hosts_table) %>%
#'   group_by(single_study = n_studies == 1) %>%
#'   summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#'
#' ## restrict to a phylogenetic group: SGS1 and lineages within 3 bp of it
#' aln <- extract_alignment()
#' m    <- toupper(as.character(aln))                       # one row per lineage
#' name <- sub("^[A-Za-z]_([^_]+).*$", "\\1", rownames(m))  # bare lineage names
#' sgs1 <- m[match("SGS1", name), ]
#' is_base  <- function(x) x %in% c("A", "C", "G", "T")
#' n_diff   <- apply(m, 1, function(s) sum(is_base(s) & is_base(sgs1) & s != sgs1))
#' sgs1_grp <- name[n_diff <= 3]
#'
#' lineage_screen() %>%
#'   filter(in_hosts_table, lineage %in% sgs1_grp) %>%
#'   group_by(single_study = n_studies == 1) %>%
#'   summarize(n = n(), mean_nonsyn = mean(n_singleton_nonsynonymous))
#' }
#' @export
lineage_screen <- function(version = "latest", reference = NULL,
                           genus = c("all", "Plasmodium", "Haemoproteus",
                                     "Leucocytozoon"),
                           min_site_coverage = 0.5, studies = TRUE) {
  genus   <- match.arg(genus)
  charmat <- .qc_char_matrix(reference, version)
  parsed  <- .lineage_parse_names(rownames(charmat))

  ## Optionally restrict the screen to a single parasite genus. The three avian
  ## haemosporidian genera are deeply divergent, so pooling them dilutes
  ## within-genus conservation; judging a lineage's singletons against only its
  ## own genus is sharper. The site profile, consensus, and singleton detection
  ## are all then computed from the genus subset alone.
  if (genus != "all") {
    keep <- !is.na(parsed$genus) & parsed$genus == genus
    if (!any(keep)) {
      stop("No sequences in the reference belong to genus '", genus, "'.",
           call. = FALSE)
    }
    charmat <- charmat[keep, , drop = FALSE]
    parsed  <- parsed[keep, , drop = FALSE]
  }

  codemat <- .qc_code_matrix(charmat)
  site_profile <- .qc_site_profile(charmat)

  n_seq <- nrow(charmat)
  min_cov_count <- ceiling(min_site_coverage * n_seq)
  singletons <- .qc_singleton_substitutions(codemat, site_profile, min_cov_count)

  out <- data.frame(
    lineage        = parsed$lineage,
    parasite_genus = parsed$genus,
    seq_length     = rowSums(codemat > 0L),
    singletons,
    seq_name       = parsed$seq_name,
    stringsAsFactors = FALSE
  )

  if (studies) {
    st <- lineage_studies(version = version)
    m  <- match(out$lineage, st$lineage)
    out$n_studies      <- st$n_studies[m]
    out$n_host_records <- st$n_host_records[m]
    out$n_countries    <- st$n_countries[m]
    out$in_hosts_table <- !is.na(m)
    ## put the study columns up front, next to the lineage/genus identifiers
    out <- out[, c("lineage", "parasite_genus", "n_studies", "n_host_records",
                   "n_countries", "in_hosts_table", "seq_length",
                   "n_singleton_substitutions",
                   "n_singleton_nonsynonymous", "n_singleton_synonymous",
                   "n_singleton_stop", "seq_name")]
  }

  out <- out[order(out$lineage), , drop = FALSE]
  rownames(out) <- NULL
  out
}
