#' BLAST-like search of a sequence against MalAvi
#'
#' Finds the MalAvi lineages most similar to a query DNA sequence against the
#' database bundled in the package. This uses \pkg{DECIPHER}: the bundled,
#' pre-built inverted index is searched with \code{DECIPHER::SearchIndex} and the
#' hits are aligned to the query with \code{DECIPHER::AlignPairs}.
#'
#' Every hit the index returns is aligned, and \code{top_n} is applied to the
#' \emph{alignment} ranking. The index score is a k-mer heuristic that does not
#' reliably rank an exact match first -- before version 1.1.2 \code{top_n} was
#' applied to it, so \code{blast_malavi(seq, top_n = 1)} could miss the lineage
#' the query came from.
#'
#' Ranking cannot separate hits the alignment cannot separate. A short query is
#' often identical to many lineages over its own length -- a 200 bp fragment of
#' \code{P_SGS1} matches about 40 lineages at 100\% -- and \code{top_n} then shows
#' an arbitrary few of them. Raise \code{top_n} when the top hits all report the
#' same \code{PercentMatch}, and see \code{\link{synonymy_report}} for why
#' distinct MalAvi names can carry the same sequence.
#'
#' \pkg{DECIPHER} (>= 3.0) and \pkg{Biostrings} are required and must be
#' installed from Bioconductor:
#' \code{BiocManager::install(c("DECIPHER", "Biostrings"))}. DECIPHER >= 3.0
#' needs R >= 4.4.
#'
#' @param sequence A DNA sequence as a single character string. Whitespace and
#'   gap (\code{-}) characters are removed; the sequence may be upper or lower case.
#' @param top_n Number of top hits to return (default 5).
#' @param version MalAvi release to search, as a date string (e.g.
#'   \code{"2026-03-23"}) or \code{"latest"} (default).
#' @return A \code{data.frame} of hits, best first, with columns \code{Lineage},
#'   \code{ProportionMatch}, \code{PercentMatch}, \code{AlignmentLength},
#'   \code{Matches}, \code{Mismatches}, \code{Score}, \code{QueryGapLength},
#'   \code{ReferenceGapLength}, \code{ReferenceLineageLength}, and
#'   \code{ReferenceFullLength}. \code{AlignmentLength} counts alignment columns,
#'   gap columns included, so \code{Matches + Mismatches + QueryGapLength +
#'   ReferenceGapLength} equals it. \code{ReferenceLineageLength} is the position
#'   in the reference lineage where the alignment ends (as reported by the original
#'   MalAvi BLAST app), whereas \code{ReferenceFullLength} is the full length of the
#'   reference lineage sequence; the two differ when the query aligns to only part of
#'   a reference. If no hits are found, a one-row data frame of \code{NA}s is
#'   returned with a warning.
#' @seealso \code{\link{extract_alignment}}
#' @examples
#' \dontrun{
#' ## requires DECIPHER (>= 3.0) and Biostrings
#' seq <- paste(as.character(extract_alignment()[1, ]), collapse = "")
#' blast_malavi(seq, top_n = 5)
#' }
#' @export
blast_malavi <- function(sequence, top_n = 5, version = "latest") {

  if (!requireNamespace("DECIPHER", quietly = TRUE) ||
      !requireNamespace("Biostrings", quietly = TRUE)) {
    stop("blast_malavi() requires the Bioconductor packages 'DECIPHER' (>= 3.0) ",
         "and 'Biostrings'.\n  Install them with: ",
         "BiocManager::install(c(\"DECIPHER\", \"Biostrings\"))", call. = FALSE)
  }
  if (utils::packageVersion("DECIPHER") < "3.0.0") {
    stop("blast_malavi() requires DECIPHER >= 3.0 (you have ",
         utils::packageVersion("DECIPHER"), "). DECIPHER >= 3.0 needs R >= 4.4.",
         call. = FALSE)
  }

  ## clean and validate the query
  q <- toupper(gsub("\\s+", "", sequence))
  q <- gsub("-", "", q)
  if (nchar(q) == 0) stop("Please supply a non-empty query sequence.", call. = FALSE)
  if (!grepl("^[ACGTNURYSWKMBDHV]+$", q)) {
    stop("Query contains invalid DNA characters.", call. = FALSE)
  }
  query <- Biostrings::DNAStringSet(q)

  ## SearchIndex() and AlignPairs() exist only in DECIPHER >= 3.0; look them up
  ## dynamically so the package checks cleanly against older DECIPHER installs.
  SearchIndex <- getExportedValue("DECIPHER", "SearchIndex")
  AlignPairs  <- getExportedValue("DECIPHER", "AlignPairs")

  ## load the bundled, pre-built BLAST database + index for this release
  blast <- .malavi_load(version, "malavi_blast_")
  db <- blast$db

  ## search the inverted index
  hits <- SearchIndex(query, blast$index)
  if (nrow(hits) == 0) {
    warning("No hits found: check your input sequence")
    return(data.frame(Lineage = NA, ProportionMatch = NA, PercentMatch = NA,
                      AlignmentLength = NA, Matches = NA, Mismatches = NA,
                      Score = NA, QueryGapLength = NA, ReferenceGapLength = NA,
                      ReferenceLineageLength = NA, ReferenceFullLength = NA))
  }

  ## Align EVERY hit the index returned, and only then cut to top_n.
  ##
  ## SearchIndex() scores k-mer similarity. That is a fast filter, not a ranking:
  ## it does not put an exact full-length match first. For a full 479 bp SGS1
  ## query the index ranks SGS1 SECOND (522.7, behind 523.4 for P_PADOM07, a
  ## 447 bp reference), so cutting to top_n on the index score -- what this
  ## function did before version 1.1.2 -- meant blast_malavi(sgs1, top_n = 1)
  ## returned P_PADOM07 at 99.776% and never reported the 100% self-match. It
  ## gets worse as the query gets shorter: for a 200 bp partial of SGS1 the index
  ## ranks the source lineage 49th, behind 39 lineages that all align at 100%.
  ## "What lineage is this?" is the natural top_n = 1 call, and it was wrong.
  ##
  ## Only the alignment can rank hits, so every hit is aligned and the cut happens
  ## after. The cost is small compared with the search itself: AlignPairs over all
  ## 5,365 references takes about 0.7 s against about 4 s for SearchIndex.
  hits <- hits[order(hits$Score, decreasing = TRUE), ]

  ## align each hit to the query to get match statistics
  aln <- AlignPairs(pattern = query, subject = db, pairs = hits, type = "values")
  lineage <- names(db)[aln$Subject]
  ## PatternGapLength and SubjectGapLength are per-hit lists of gap lengths; sum
  ## each to a single number. Both are reported: with only the query side, the
  ## columns of a query carrying an insertion do not add up (Matches + Mismatches
  ## + QueryGapLength = 478 against an AlignmentLength of 481) and the user has no
  ## way to see that the missing 3 columns are gaps opened in the reference.
  query_gap     <- vapply(aln$PatternGapLength, function(x) sum(as.numeric(x)), numeric(1))
  reference_gap <- vapply(aln$SubjectGapLength, function(x) sum(as.numeric(x)), numeric(1))

  out <- data.frame(
    Lineage                = lineage,
    ProportionMatch        = paste(aln$Matches, aln$AlignmentLength, sep = "/"),
    PercentMatch           = round(aln$Matches / aln$AlignmentLength * 100, 3),
    AlignmentLength        = aln$AlignmentLength,
    Matches                = aln$Matches,
    Mismatches             = aln$Mismatches,
    Score                  = aln$Score,
    QueryGapLength         = query_gap,
    ReferenceGapLength     = reference_gap,
    ReferenceLineageLength = aln$SubjectEnd,
    ReferenceFullLength    = Biostrings::width(db)[aln$Subject],
    stringsAsFactors       = FALSE
  )
  ## rank by the alignment, then take the top_n the caller asked for
  out <- out[order(out$Score, decreasing = TRUE), ]
  out <- out[seq_len(min(top_n, nrow(out))), , drop = FALSE]
  rownames(out) <- NULL
  .malavi_attach_meta(out, malavi_version = .malavi_resolve_version(version),
                      method = "DECIPHER SearchIndex + AlignPairs")
}
