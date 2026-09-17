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
    table_genus  = gls$GENUS_NAME[match(clean_names(names), gls$LINEAGE_NAME)],
    ## an environment, so a check that computes something expensive can leave
    ## it for its describe() to reuse (the list itself is copied, this is not)
    cache        = new.env(parent = emptyenv())
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
    ),

    list(
      title = "Parasite genus contradicts the nearest sequences",
      check = function(ctx) {
        out <- .malavi_genus_outliers_ctx(ctx)
        out$lineage
      },
      describe = function(affected, ctx) {
        ## The short form (Vincenzo, 2026-09-17): the count, the rule once, and
        ## the names. The per-lineage detail -- neighbour counts, the nearest
        ## sequence and its distance -- is in the table the check computes
        ## (.malavi_genus_outliers_ctx), not in the sentence, because the
        ## sentence is what the website's known-issues card shows and thirty
        ## clauses of it would be a wall.
        n <- length(affected)
        paste0(
          n, if (n == 1) " lineage is" else " lineages are",
          " listed under a genus that ", if (n == 1) "its" else "their",
          " sequences contradict: of the sequences within ",
          .malavi_count_word(.MALAVI_GENUS_MAX_MISMATCH),
          " mismatches of ", if (n == 1) "it" else "each",
          " (over at least ", .MALAVI_GENUS_MIN_COMPARABLE,
          " compared positions), three quarters or more belong to another genus. ",
          if (n == 1) "Lineage is: " else "Lineages are: ",
          .malavi_lineage_list(affected), ".")
      }
    )
  )
}

## ---- the genus-versus-neighbours check -------------------------------------
##
## A lineage's GENUS_NAME is typed in by hand; its sequence is not. When every
## sequence within a few substitutions of a lineage carries a different genus,
## the label is wrong far more often than the biology is (POEPAL01, entered as
## Haemoproteus, is SGS1 with one substitution and is "Plasmodium sp. PAPA01" at
## GenBank; 32 lineages in the 2026-09-15 release are like it). This check
## re-derives that list from the loaded release.
##
## The thresholds are deliberately conservative, so that only a lineage sitting
## INSIDE a cluster of another genus is reported:

## a lineage must have at least this many determined (A/C/G/T) positions, and a
## neighbour must share at least this many with it, before the comparison counts.
## Same value and reasoning as MIN_INFORMATIVE_FOR_IDENTITY elsewhere in the
## rebuild: MalAvi's own standard for calling two sequences the same lineage.
.MALAVI_GENUS_MIN_COMPARABLE <- 300L
## a neighbour is a sequence within this many mismatches over the comparable
## positions (pairwise deletion). Five is ~1 % of the barcode; genera differ by
## an order of magnitude more.
.MALAVI_GENUS_MAX_MISMATCH <- 5L
## at least this many neighbours, and at least this share of them of one genus
## other than the lineage's own. Three stops a single mislabelled neighbour
## from indicting a correctly labelled lineage; the share tolerates a couple of
## other mislabelled lineages in the same cluster.
.MALAVI_GENUS_MIN_NEIGHBORS <- 3L
.MALAVI_GENUS_MIN_SHARE <- 0.75

## The computation, on an upper-case character alignment (rows = sequences,
## columns = positions; A/C/G/T are determined, anything else is not) and a
## genus per row. Pure, so it can be tested on a hand-made alignment. Returns
## one row per outlier, in alignment order.
##
## Mismatch counts over comparable positions come from ape::dist.dna with
## pairwise deletion (its "known base" is exactly A/C/G/T, checked against the
## direct definition on 200 random pairs of the 2026-09-15 release, 0
## disagreements). Comparable counts are L - u_i - u_j + |u_i AND u_j| over the
## undetermined-position indicators, with the last term a matrix product over
## the rows that have any undetermined position. About 25 s for the release in
## plain R, which is why the release build stores the result in the bundle
## (see .malavi_genus_outliers_ctx).
.malavi_genus_outliers <- function(charmat, genus, row_names = rownames(charmat),
                                   min_comparable = .MALAVI_GENUS_MIN_COMPARABLE,
                                   max_mismatch   = .MALAVI_GENUS_MAX_MISMATCH,
                                   min_neighbors  = .MALAVI_GENUS_MIN_NEIGHBORS,
                                   min_share      = .MALAVI_GENUS_MIN_SHARE) {
  n <- nrow(charmat)
  empty <- data.frame(
    lineage = character(0), genus = character(0),
    n_neighbors = integer(0), n_neighbors_other = integer(0),
    neighbor_genus = character(0), nearest = character(0),
    nearest_genus = character(0), nearest_mismatches = integer(0),
    nearest_comparable = integer(0), stringsAsFactors = FALSE)
  if (n < 2L) return(empty)

  determined <- matrix(charmat %in% c("A", "C", "G", "T"), n, ncol(charmat))
  n_determined <- rowSums(determined)
  ## a genus that is unknown can neither be contradicted nor contradict
  known_genus <- !is.na(genus) &
    genus %in% c("Plasmodium", "Haemoproteus", "Leucocytozoon")

  ## mismatches: sites where both are determined and the bases differ
  dna <- ape::as.DNAbin(matrix(tolower(charmat), n, ncol(charmat)))
  mismatch <- ape::dist.dna(dna, model = "N", pairwise.deletion = TRUE,
                            as.matrix = TRUE)
  ## comparable: sites determined in both
  undetermined <- ncol(charmat) - n_determined
  comparable <- ncol(charmat) - outer(undetermined, undetermined, "+")
  partial <- which(undetermined > 0)
  if (length(partial) > 1L) {
    U <- (!determined[partial, , drop = FALSE]) + 0   # numeric for BLAS
    comparable[partial, partial] <- comparable[partial, partial] + tcrossprod(U)
  } else if (length(partial) == 1L) {
    comparable[partial, partial] <- comparable[partial, partial] + undetermined[partial]
  }

  rows <- list()
  for (i in seq_len(n)) {
    if (n_determined[i] < min_comparable || !known_genus[i]) next
    close <- which(comparable[i, ] >= min_comparable &
                   mismatch[i, ] <= max_mismatch & known_genus)
    close <- close[close != i]
    if (length(close) < min_neighbors) next
    tab <- sort(table(genus[close]), decreasing = TRUE)
    top <- names(tab)[1]
    if (top == genus[i] || tab[[1]] < min_share * length(close)) next
    ## the nearest neighbour: lowest mismatch rate, then most comparable
    ord <- order(mismatch[i, close] / comparable[i, close], -comparable[i, close])
    j <- close[ord[1]]
    rows[[length(rows) + 1L]] <- data.frame(
      lineage = row_names[i], genus = genus[i],
      n_neighbors = length(close), n_neighbors_other = as.integer(tab[[1]]),
      neighbor_genus = top, nearest = row_names[j], nearest_genus = genus[j],
      nearest_mismatches = as.integer(round(mismatch[i, j])),
      nearest_comparable = as.integer(round(comparable[i, j])),
      stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(empty)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

## The genus of each alignment row, as the check wants it: the Grand Lineage
## Summary's GENUS_NAME by lineage name, the alignment prefix only when the
## table has no row for the name, and NA for "N/A" (a missing genus is its own
## issue and can neither be contradicted nor contradict). Shared by the release
## build (data-raw/process_release.R) and the call-time fallback, so the table
## in the bundle and a fresh computation cannot differ in what "genus" means.
.malavi_genus_by_row <- function(alignment_names, gls) {
  lineage <- clean_names(alignment_names)
  genus <- gls$GENUS_NAME[match(lineage, gls$LINEAGE_NAME)]
  prefix <- unname(c(P = "Plasmodium", H = "Haemoproteus",
                     L = "Leucocytozoon")[substr(alignment_names, 1, 1)])
  genus[is.na(genus)] <- prefix[is.na(genus)]
  genus[!is.na(genus) & genus == "N/A"] <- NA_character_
  genus
}

## The table for one release bundle: rows of the bundle's alignment against the
## bundle's Grand Lineage Summary. This is what the release build stores as
## bundle$genus_outliers, and what the fallback computes when a bundle predates
## the check.
.malavi_genus_outliers_bundle <- function(bundle) {
  charmat <- toupper(as.character(bundle$alignment))
  gls <- bundle[["Grand Lineage Summary"]]
  if (is.null(gls)) gls <- bundle$grand_lineage_summary
  genus <- .malavi_genus_by_row(rownames(charmat), gls)
  .malavi_genus_outliers(charmat, genus, row_names = clean_names(rownames(charmat)))
}

## The same on the release in `ctx`, computed once per context: the check and
## the describe both need the table, and `ctx$cache` is an environment so the
## second call finds what the first stored. The bundle's stored table is used
## when the release build wrote one; otherwise it is computed here, slowly.
.malavi_genus_outliers_ctx <- function(ctx) {
  if (!is.null(ctx$cache$genus_outliers)) return(ctx$cache$genus_outliers)
  bundle <- .malavi_load(ctx$version)
  out <- bundle$genus_outliers
  if (is.null(out)) {
    message("This release bundle carries no genus-outlier table; computing it ",
            "from the alignment (about half a minute).")
    out <- .malavi_genus_outliers_bundle(bundle)
  }
  ctx$cache$genus_outliers <- out
  out
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
#' @return Invisibly, a two-column \code{data.frame} with the \code{title} and
#'   \code{text} of each issue found, in the order printed. Called mainly for the
#'   report it prints; the return value exists so the website can render the same
#'   list from the same source rather than repeating it in HTML.
#' @seealso \code{\link{extract_table}}, \code{\link{extract_alignment}}
#' @examples
#' malavi_issues()
#' @export
malavi_issues <- function(version = "latest") {
  resolved <- .malavi_resolve_version(version)

  cat("Known issues in the current MalAvi data release; MalAvi version ",
      resolved, "\n\n", sep = "")

  ctx    <- .malavi_issue_context(version)
  titles <- character(0)
  texts  <- character(0)

  for (issue in .malavi_issue_registry()) {
    affected <- unique(as.character(issue$check(ctx)))
    ## An issue the release no longer shows is not a known issue: leave it out
    ## rather than printing it as resolved.
    if (length(affected) == 0) next

    text   <- issue$describe(affected, ctx)
    titles <- c(titles, issue$title)
    texts  <- c(texts, text)

    cat(issue$title, "\n", sep = "")
    writeLines(paste0("   ", strwrap(text, width = 74)))
    cat("\n")
  }

  if (length(texts) == 0) cat("No known issues were found in this release.\n")

  invisible(data.frame(title = titles, text = texts, stringsAsFactors = FALSE))
}
