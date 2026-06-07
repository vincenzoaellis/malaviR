#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## process_release.R — turn a MalAvi release zip into the bundled package data.
##
## Recurring maintainer workflow when a new MalAvi release arrives:
##   1. Drop the new MalAvi_<date>.zip into data-raw/
##   2. Run this script on a compute node:
##        Rscript data-raw/process_release.R                # newest zip in data-raw/
##        Rscript data-raw/process_release.R data-raw/MalAvi_2026-03-23.zip
##   3. Commit the new inst/extdata/*.rds (+ the archived zip) and push.
##
## A release zip contains one folder with 5 .xlsx tables + 1 .fas alignment, all
## stamped with the release date. This script writes two files per release:
##   inst/extdata/malavi_db_<date>.rds    — list of 5 tables + alignment (DNAbin)
##   inst/extdata/malavi_blast_<date>.rds — list(db = ungapped DNAStringSet,
##                                               index = DECIPHER InvertedIndex)
##
## Decision: ship the latest release only -> older malavi_db_*/malavi_blast_*
## .rds files in inst/extdata are pruned (set KEEP below to keep more).
## The original .zip stays in data-raw/ for provenance.
## ---------------------------------------------------------------------------

## prepend the personal library for the *running* R version (harmless if absent)
rver <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
userlib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", rver)
if (dir.exists(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages({
  library(readxl)
  library(ape)
})

KEEP <- 1L  # number of most-recent releases to keep bundled in inst/extdata

## --- locate repo root + paths -------------------------------------------------
## Run from the repo root (or anywhere; we resolve relative to this file if possible).
repo   <- normalizePath(".")
raw    <- file.path(repo, "data-raw")
outdir <- file.path(repo, "inst", "extdata")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## --- pick the release zip -----------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  zip_path <- normalizePath(args[1])
} else {
  zips <- list.files(raw, pattern = "^MalAvi_.*\\.zip$", full.names = TRUE)
  if (length(zips) == 0) stop("No MalAvi_*.zip found in data-raw/")
  zip_path <- zips[order(basename(zips), decreasing = TRUE)][1]  # newest by date in name
}
date_tag <- sub("^MalAvi_(.*)\\.zip$", "\\1", basename(zip_path))
message("Release zip : ", basename(zip_path))
message("Date tag    : ", date_tag)

## --- extract to a temp dir ----------------------------------------------------
tmp <- file.path(tempdir(), paste0("malavi_", date_tag))
unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
utils::unzip(zip_path, exdir = tmp)
files <- list.files(tmp, recursive = TRUE, full.names = TRUE)

## helper: find the one file whose name starts with a given prefix
pick <- function(prefix, ext) {
  hit <- grep(sprintf("/%s.*\\.%s$", prefix, ext), files, value = TRUE, ignore.case = TRUE)
  if (length(hit) != 1) stop("Expected exactly one ", prefix, "*.", ext, " file; found ", length(hit))
  hit
}

## --- read the 5 tables + alignment -------------------------------------------
message("Reading tables...")
tables <- list(
  grand_lineage_summary = readxl::read_excel(pick("GrandLineageSummary", "xlsx")),
  hosts_and_sites       = readxl::read_excel(pick("Hosts_and_Sites",     "xlsx")),
  morpho_species        = readxl::read_excel(pick("MorphoSpecies",       "xlsx")),
  references            = readxl::read_excel(pick("References",          "xlsx")),
  vector_data           = readxl::read_excel(pick("VectorData",          "xlsx"))
)
tables <- lapply(tables, as.data.frame)  # plain data.frames, not tibbles
for (nm in names(tables)) message(sprintf("  %-22s %d rows x %d cols",
                                          nm, nrow(tables[[nm]]), ncol(tables[[nm]])))

message("Reading alignment...")
fas <- pick("MalAvi", "fas")
alignment <- ape::read.dna(fas, format = "fasta")
message(sprintf("  alignment              %d seqs x %d bp", nrow(alignment), ncol(alignment)))

db_bundle <- c(tables, list(alignment = alignment, version = date_tag))
db_out <- file.path(outdir, sprintf("malavi_db_%s.rds", date_tag))
saveRDS(db_bundle, db_out, compress = "xz")
message("Wrote ", db_out, " (", round(file.size(db_out) / 1e3), " KB)")

## --- build the DECIPHER BLAST index ------------------------------------------
## Mirrors MalAviBLAST/app.R: RemoveGaps() then IndexSeqs() on the ungapped db.
if (requireNamespace("DECIPHER", quietly = TRUE) &&
    requireNamespace("Biostrings", quietly = TRUE)) {
  message("Building DECIPHER BLAST index (this is the slow step)...")
  db_full <- Biostrings::readDNAStringSet(fas)
  db_ungapped <- DECIPHER::RemoveGaps(db_full)
  index <- DECIPHER::IndexSeqs(db_ungapped, sensitivity = 0.99,
                               percentIdentity = 99, patternLength = 479)
  blast_bundle <- list(db = db_ungapped, index = index, version = date_tag)
  blast_out <- file.path(outdir, sprintf("malavi_blast_%s.rds", date_tag))
  saveRDS(blast_bundle, blast_out, compress = "xz")
  message("Wrote ", blast_out, " (", round(file.size(blast_out) / 1e6, 1), " MB)")
} else {
  message("DECIPHER/Biostrings not available -- skipped BLAST index.")
  message("Install them and re-run, or build malavi_blast_", date_tag, ".rds separately.")
}

## --- prune to the KEEP most-recent releases ----------------------------------
prune <- function(pattern) {
  f <- list.files(outdir, pattern = pattern, full.names = TRUE)
  if (length(f) > KEEP) {
    drop <- f[order(basename(f), decreasing = TRUE)][-(seq_len(KEEP))]
    file.remove(drop)
    for (d in drop) message("Pruned old: ", basename(d))
  }
}
prune("^malavi_db_.*\\.rds$")
prune("^malavi_blast_.*\\.rds$")

message("DONE process_release.R for ", date_tag)
