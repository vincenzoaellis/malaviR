#' BLAST-like search of a sequence against MalAvi
#'
#' Finds the MalAvi lineages most similar to a query DNA sequence against the
#' database bundled in the package. This uses \pkg{DECIPHER}: the bundled,
#' pre-built inverted index is searched with \code{DECIPHER::SearchIndex} and the
#' top hits are aligned to the query with \code{DECIPHER::AlignPairs}.
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
#'   \code{ReferenceLineageLength}, and \code{ReferenceFullLength}.
#'   \code{ReferenceLineageLength} is the position in the reference lineage where
#'   the alignment ends (as reported by the original MalAvi BLAST app), whereas
#'   \code{ReferenceFullLength} is the full length of the reference lineage
#'   sequence; the two differ when the query aligns to only part of a reference.
#'   If no hits are found, a one-row data frame of \code{NA}s is returned with a
#'   warning.
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
                      Score = NA, QueryGapLength = NA, ReferenceLineageLength = NA,
                      ReferenceFullLength = NA))
  }

  ## keep the top_n hits by score
  hits <- hits[order(hits$Score, decreasing = TRUE), ]
  n <- min(top_n, nrow(hits))
  hits <- hits[seq_len(n), ]

  ## align each hit to the query to get match statistics
  aln <- AlignPairs(pattern = query, subject = db, pairs = hits, type = "values")
  lineage <- names(db)[aln$Subject]
  ## PatternGapLength is a per-hit list of gap lengths; sum to a single number
  query_gap <- vapply(aln$PatternGapLength, function(x) sum(as.numeric(x)), numeric(1))

  out <- data.frame(
    Lineage                = lineage,
    ProportionMatch        = paste(aln$Matches, aln$AlignmentLength, sep = "/"),
    PercentMatch           = round(aln$Matches / aln$AlignmentLength * 100, 3),
    AlignmentLength        = aln$AlignmentLength,
    Matches                = aln$Matches,
    Mismatches             = aln$Mismatches,
    Score                  = aln$Score,
    QueryGapLength         = query_gap,
    ReferenceLineageLength = aln$SubjectEnd,
    ReferenceFullLength    = Biostrings::width(db)[aln$Subject],
    stringsAsFactors       = FALSE
  )
  out <- out[order(out$Score, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}
