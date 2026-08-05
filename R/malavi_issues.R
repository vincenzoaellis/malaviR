## A registry of known problems in the MalAvi data as distributed. The bundled
## tables are shipped verbatim -- malaviR does not silently correct them -- so
## this file is how the package remembers what is wrong with them.
##
## Every entry is a CHECK plus a DESCRIBE:
##
##   * the CHECK re-derives the issue from the release the user has actually
##     loaded, returning the affected lineages. This is what keeps the list
##     honest over time: when a future release fixes something, the check stops
##     finding it and the issue drops out of the report, without anyone editing
##     this file.
##   * the DESCRIBE turns those lineages into the sentence the user reads. The
##     text is written from the check result rather than stored, so it can never
##     say "7 lineages" about a release that has nine.
##
## The report is deliberately plain: a heading, then one title and one sentence
## per issue. There is no status tag, no provenance line, no count column and no
## workaround code -- this is a list to read, not a table to program against.
##
## To add an issue: append a list to .malavi_issue_registry() with a title, a
## check that returns the affected lineages from `ctx`, a describe that turns
## them into a sentence, and a test.

## Small number words, so a description reads "two morphospecies" rather than
## "2 morphospecies". Falls back to the digits above ten.
.malavi_count_word <- function(n) {
  words <- c("one", "two", "three", "four", "five",
             "six", "seven", "eight", "nine", "ten")
  if (n >= 1 && n <= length(words)) words[n] else as.character(n)
}

## Render a series of binomials the way a biologist writes them: the genus in
## full the first time, abbreviated afterwards when it does not change.
## c("Haemoproteus minutus", "Haemoproteus asymmetricus") becomes
## "Haemoproteus minutus and H. asymmetricus".
.malavi_species_series <- function(species) {
  species <- species[!is.na(species) & nzchar(species)]
  if (length(species) == 0) return("")
  first_genus <- sub("^([A-Za-z]+).*$", "\\1", species[1])
  shown <- species
  if (length(species) > 1) {
    later      <- species[-1]
    same_genus <- sub("^([A-Za-z]+).*$", "\\1", later) == first_genus
    later[same_genus] <- sub("^([A-Za-z]+)", paste0(substr(first_genus, 1, 1), "."),
                             later[same_genus])
    shown <- c(species[1], later)
  }
  if (length(shown) == 1) return(shown)
  if (length(shown) == 2) return(paste(shown, collapse = " and "))
  paste0(paste(shown[-length(shown)], collapse = ", "),
         " and ", shown[length(shown)])
}

## A plain comma-separated list of lineage names, always in a stable order.
.malavi_lineage_list <- function(x) paste(sort(unique(x)), collapse = ", ")

## Build the data every check needs, once, rather than per issue.
.malavi_issue_context <- function(version) {
  gls   <- extract_table("Grand Lineage Summary", version = version)
  aln   <- extract_alignment(version = version)
  names <- rownames(aln)

  ## genus implied by the one-letter prefix of an alignment name, NA when the
  ## prefix is not one of the three recognized genus codes
  prefix_genus <- unname(c(P = "Plasmodium", H = "Haemoproteus",
                           L = "Leucocytozoon")[substr(names, 1, 1)])

  list(
    version      = version,
    gls          = gls,
    seq_names    = names,
    seq_lineage  = clean_names(names),
    prefix_genus = prefix_genus,
    ## GENUS_NAME as recorded for each alignment sequence; "N/A" left as-is so
    ## the missing-genus issue can see it
    table_genus  = gls$GENUS_NAME[match(clean_names(names), gls$LINEAGE_NAME)]
  )
}

## The registry itself. One list per issue.
.malavi_issue_registry <- function() {
  list(

    list(
      title = "Alignment genus prefix contradicts GENUS_NAME",
      check = function(ctx) {
        bad <- !is.na(ctx$prefix_genus) & !is.na(ctx$table_genus) &
          ctx$table_genus != "N/A" & ctx$prefix_genus != ctx$table_genus
        ctx$seq_names[bad]
      },
      describe = function(affected, ctx) {
        idx <- match(affected, ctx$seq_names)
        paste(paste0(
          affected, " has the ", ctx$prefix_genus[idx], " prefix (",
          substr(affected, 1, 1), "_) in the alignment, but in the Grand ",
          "Lineage Summary table it is listed as ", ctx$table_genus[idx]),
          collapse = ". ")
      }
    ),

    list(
      title = "No parasite genus listed",
      check = function(ctx) {
        na_genus <- !is.na(ctx$gls$GENUS_NAME) & ctx$gls$GENUS_NAME == "N/A"
        ctx$gls$LINEAGE_NAME[na_genus]
      },
      describe = function(affected, ctx) {
        paste0(
          length(affected), if (length(affected) == 1) " lineage" else " lineages",
          " in the Grand Lineage Summary table do not have a parasite genus ",
          "listed and are given the prefix N_ in the alignment file; ",
          if (length(affected) == 1) "lineage is: " else "lineages are: ",
          .malavi_lineage_list(affected))
      }
    ),

    list(
      title = "One lineage listed as more than one morphospecies",
      check = function(ctx) {
        ln <- ctx$gls$LINEAGE_NAME
        sort(unique(ln[duplicated(ln)]))
      },
      describe = function(affected, ctx) {
        parts <- vapply(affected, function(lin) {
          sp <- unique(ctx$gls$SPECIES_NAME[ctx$gls$LINEAGE_NAME == lin])
          sp <- sp[!is.na(sp) & nzchar(sp)]
          paste0("Lineage ", lin, " is listed as ", .malavi_count_word(length(sp)),
                 " morphospecies, ", .malavi_species_series(sp))
        }, character(1))
        paste(parts, collapse = ". ")
      }
    ),

    list(
      title = "Lineages in the alignment but not the table, or the reverse",
      check = function(ctx) {
        in_aln   <- unique(ctx$seq_lineage)
        in_table <- unique(ctx$gls$LINEAGE_NAME)
        c(setdiff(in_aln, in_table), setdiff(in_table, in_aln))
      },
      describe = function(affected, ctx) {
        in_aln     <- unique(ctx$seq_lineage)
        in_table   <- unique(ctx$gls$LINEAGE_NAME)
        only_aln   <- sort(setdiff(in_aln, in_table))
        only_table <- sort(setdiff(in_table, in_aln))
        parts <- character(0)
        if (length(only_aln) > 0)
          parts <- c(parts, paste0(
            "The following lineages are found in the alignment but not the ",
            "Grand Lineage Summary Table: ", .malavi_lineage_list(only_aln), "."))
        if (length(only_table) > 0)
          parts <- c(parts, paste0(
            "The following lineages are found in the Grand Lineage Summary ",
            "table but not the alignment: ", .malavi_lineage_list(only_table), "."))
        paste(parts, collapse = " ")
      }
    )
  )
}

#' Known issues in the MalAvi data
#'
#' Prints the list of known problems in the MalAvi database itself, as
#' distributed. \code{malaviR} ships the MalAvi tables \strong{verbatim} and does
#' not correct them, so this list is how the package reports what is known to be
#' wrong with them.
#'
#' Each issue is \strong{re-derived from the release you have loaded} rather than
#' stored as text, and the lineages named in each sentence are the ones found in
#' that release. An issue that a future MalAvi release fixes stops being found
#' and simply drops out of the report, without anything in the package changing.
#'
#' These are upstream data problems. Fixing them properly means fixing them in a
#' future MalAvi release.
#'
#' @param version MalAvi release to report on; a date string or \code{"latest"}
#'   (default).
#' @return The issue text, invisibly, as a character vector. Called for the
#'   report it prints.
#' @seealso \code{\link{extract_table}}, \code{\link{extract_alignment}}
#' @examples
#' malavi_issues()
#' @export
malavi_issues <- function(version = "latest") {
  resolved <- .malavi_resolve_version(version)

  cat("Known issues in the current MalAvi data release; MalAvi version ",
      resolved, "\n\n", sep = "")

  ctx   <- .malavi_issue_context(version)
  found <- character(0)

  for (issue in .malavi_issue_registry()) {
    affected <- unique(as.character(issue$check(ctx)))
    ## An issue the release no longer shows is not a known issue: leave it out
    ## rather than printing it as resolved.
    if (length(affected) == 0) next

    text  <- issue$describe(affected, ctx)
    found <- c(found, text)

    cat(issue$title, "\n", sep = "")
    writeLines(paste0("   ", strwrap(text, width = 74)))
    cat("\n")
  }

  if (length(found) == 0) cat("No known issues were found in this release.\n")

  invisible(found)
}
