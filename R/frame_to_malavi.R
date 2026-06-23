#' Place primer-trimmed ASVs into the MalAvi barcode frame
#'
#' Pads fully primer-trimmed amplicon sequence variants (ASVs) into the MalAvi
#' reference frame (479 bp for the bundled data) expected by
#' \code{\link{lineage_qc}}. The MalAvi cyt \emph{b}
#' barcode frame is the same length for all haemosporidian genera, but different
#' primer pairs bracket different sub-windows of it, so a clean primer-trimmed
#' ASV is \strong{not} the same length for every assay. This function is
#' parameterized by the sub-window the primers cover, so it works for any MalAvi
#' cyt \emph{b} primer set, not just the two built in.
#'
#' The covered window is given as 1-based inclusive positions in the reference
#' frame (\code{frame_start}, \code{frame_end}); everything outside it falls
#' under a primer footprint and is padded with \code{pad_char} (\code{"N"} by
#' default), which \code{lineage_qc()} treats as missing.
#' Bases under the primer footprint come from the primer, not the template, and
#' so are genuinely unobservable from these amplicons; representing them as
#' \code{N} is the honest placement (do not retain the primer bases).
#'
#' Two primer sets are built in as shortcuts (supply \code{primer = }):
#' \describe{
#'   \item{\code{"haem"} (HaemF / HaemR2)}{covers frame positions \strong{2-479}
#'     (clean ASV 478 bp); pad 1 \code{N} on the left.}
#'   \item{\code{"leuc"} (HaemFL / HaemR2L)}{covers frame positions \strong{2-477}
#'     (clean ASV 476 bp); the leuc reverse primer binds two bases further into
#'     the barcode, so pad 1 \code{N} left and 2 \code{N} right.}
#' }
#' For any other primer set, leave \code{primer = NULL} and pass the window
#' explicitly with \code{frame_start} / \code{frame_end}. The window must be
#' determined from the primers' actual binding positions (in-silico PCR against
#' reference mtDNA + mapping onto the cyt \emph{b} CDS); this function does not
#' infer it.
#'
#' Only ASVs of the expected clean length (\code{frame_end - frame_start + 1})
#' are framed. Anything else (including a sequence that is already the full
#' reference length) is treated as off-length: such sequences are typically
#' indels, chimeras, or off-target product and should be placed manually with
#' \code{\link{blast_malavi}} rather than forced into the frame -- forcing a
#' wrong-length sequence into the frame shifts its reading frame and produces
#' spurious stop codons in the codon-aware QC.
#'
#' @param seqs Character vector of fully primer-trimmed ASV sequences (any case;
#'   whitespace is stripped). Names, if present, are preserved on the output.
#' @param primer Optional built-in primer set: \code{"haem"} or \code{"leuc"}.
#'   Supply this \strong{or} \code{frame_start}/\code{frame_end}, not both. If
#'   \code{NULL} (default), the window must be given explicitly.
#' @param frame_start,frame_end 1-based inclusive positions, in the reference
#'   frame, of the first and last template base the primers leave intact (i.e.
#'   the sub-window the clean ASV occupies). Required when \code{primer} is
#'   \code{NULL}.
#' @param reference_length Length of the MalAvi reference frame (default 479 for
#'   the bundled data).
#' @param on_off_length How to handle ASVs that are not the expected clean
#'   length: \code{"set_na"} (default) returns \code{NA} for them (so the caller
#'   can drop them and screen only the framed ones, and report the rest as
#'   off-length); \code{"error"} stops; \code{"keep"} returns the original
#'   (un-framed) sequence.
#' @param pad_char Character used for padding (default \code{"N"}).
#' @return A character vector the same length as \code{seqs}: each clean ASV
#'   padded to \code{reference_length} bp, off-length ASVs handled per
#'   \code{on_off_length}. Off-length ASVs trigger a warning (or error) giving
#'   their count and the expected clean length.
#' @seealso \code{\link{lineage_qc}},
#'   \code{\link{blast_malavi}}
#' @examples
#' ## built-in leuc assay (clean ASV 476 bp) -> 479 bp frame
#' leuc_asv <- paste(rep("A", 476), collapse = "")
#' framed <- frame_to_malavi(leuc_asv, primer = "leuc")
#' nchar(framed)            # 479
#'
#' ## an arbitrary primer set covering frame positions 24-479 (clean ASV 456 bp)
#' asv456 <- paste(rep("C", 456), collapse = "")
#' frame_to_malavi(asv456, frame_start = 24, frame_end = 479)
#' @export
frame_to_malavi <- function(seqs,
                            primer = NULL,
                            frame_start = NULL,
                            frame_end = NULL,
                            reference_length = 479L,
                            on_off_length = c("set_na", "error", "keep"),
                            pad_char = "N") {
  on_off_length <- match.arg(on_off_length)
  if (!is.character(seqs))
    stop("`seqs` must be a character vector of ASV sequences.", call. = FALSE)
  reference_length <- as.integer(reference_length)

  ## --- resolve the covered window ------------------------------------------
  ## built-in primer frames (1-based inclusive positions in the 479 bp frame);
  ## see the cyt b primer-frame reference for how these were established.
  builtin <- list(
    haem = c(frame_start = 2L, frame_end = 479L),
    leuc = c(frame_start = 2L, frame_end = 477L)
  )
  have_primer <- !is.null(primer)
  have_window <- !is.null(frame_start) || !is.null(frame_end)
  if (have_primer && have_window) {
    stop("Supply either `primer` or `frame_start`/`frame_end`, not both.",
         call. = FALSE)
  }
  if (have_primer) {
    if (!is.character(primer) || length(primer) != 1L || !(primer %in% names(builtin))) {
      stop("`primer` must be one of: ", paste(shQuote(names(builtin)), collapse = ", "),
           ". For any other primer set leave `primer = NULL` and pass ",
           "`frame_start`/`frame_end`.", call. = FALSE)
    }
    frame_start <- builtin[[primer]][["frame_start"]]
    frame_end   <- builtin[[primer]][["frame_end"]]
  } else {
    if (is.null(frame_start) || is.null(frame_end)) {
      stop("Supply `primer` (\"haem\"/\"leuc\") or both `frame_start` and ",
           "`frame_end` for a custom primer set.", call. = FALSE)
    }
    frame_start <- as.integer(frame_start)
    frame_end   <- as.integer(frame_end)
  }

  ## --- validate the window -------------------------------------------------
  if (is.na(frame_start) || is.na(frame_end) ||
      frame_start < 1L || frame_end > reference_length || frame_start > frame_end) {
    stop(sprintf(paste0("Invalid window: frame_start=%s, frame_end=%s must satisfy ",
                        "1 <= frame_start <= frame_end <= reference_length (%d)."),
                 frame_start, frame_end, reference_length), call. = FALSE)
  }

  left_pad  <- frame_start - 1L
  right_pad <- reference_length - frame_end
  clean_len <- frame_end - frame_start + 1L

  ## --- frame the clean ASVs ------------------------------------------------
  s        <- toupper(trimws(seqs))
  is_clean <- nchar(s) == clean_len

  framed <- rep(NA_character_, length(s))
  framed[is_clean] <- paste0(strrep(pad_char, left_pad),
                             s[is_clean],
                             strrep(pad_char, right_pad))

  n_off <- sum(!is_clean)
  if (n_off > 0L) {
    msg <- sprintf(paste0("%d of %d sequence(s) are not the expected clean ",
                          "length (%d bp, = frame positions %d-%d) and were not ",
                          "framed. Off-length ASVs are likely indels, chimeras, ",
                          "or off-target product -- place them manually with ",
                          "blast_malavi()."),
                   n_off, length(s), clean_len, frame_start, frame_end)
    if (on_off_length == "error") stop(msg, call. = FALSE)
    warning(msg, call. = FALSE)
    if (on_off_length == "keep") framed[!is_clean] <- s[!is_clean]
    ## "set_na" (default): leave NA so the caller filters to screenable ASVs
  }

  names(framed) <- names(seqs)
  framed
}
