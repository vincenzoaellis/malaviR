#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## build_taxonomy.R — build the bundled clootl taxonomy snapshot + host crosswalk.
##
## Run on a compute node whenever you want to refresh the taxonomy (e.g. after a
## new MalAvi release, or when clootl publishes a newer taxonomy year):
##   Rscript data-raw/build_taxonomy.R
##
## It writes two things:
##   R/sysdata.rda   — internal: clootl_ref (trimmed clootl taxonomy) + clootl_year
##   data/taxonomy.rda — exported dataset: MalAvi host -> eBird crosswalk
##
## We bundle a snapshot of the *current* clootl taxonomy (the latest year clootl
## serves) because older years/phylogenies become hard to access; matching to the
## bundled current taxonomy keeps results reproducible.
## ---------------------------------------------------------------------------

rver <- paste(R.version$major, sub("\\..*", "", R.version$minor), sep = ".")
userlib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", rver)
if (dir.exists(userlib)) .libPaths(c(userlib, .libPaths()))

suppressPackageStartupMessages(library(clootl))

## --- find the latest clootl taxonomy year that is actually available ----------
clootl_year <- NULL
for (yr in (as.integer(format(Sys.Date(), "%Y")) + 1):2021) {
  ok <- tryCatch({ taxonomyGet(taxonomy_year = yr); TRUE },
                 error = function(e) FALSE)
  if (ok) { clootl_year <- yr; break }
}
if (is.null(clootl_year)) stop("Could not retrieve any clootl taxonomy year.")
message("Using clootl taxonomy year: ", clootl_year)

tax <- taxonomyGet(taxonomy_year = clootl_year)

## --- trim to the columns match_taxonomy() needs; clean stray whitespace -------
## the order column has been named ORDER or ORDER1 in different clootl years
clean <- function(x) trimws(as.character(x))
order_col <- if ("ORDER" %in% names(tax)) tax$ORDER else tax$ORDER1
clootl_ref <- data.frame(
  SCI_NAME      = clean(tax$SCI_NAME),
  SPECIES_CODE  = clean(tax$SPECIES_CODE),
  ott_name      = clean(tax$ott_name),
  ORDER1        = clean(order_col),
  FAMILY        = clean(tax$FAMILY),
  IOC_name      = clean(tax$IOC_name),
  Birdlife_name = clean(tax$Birdlife_name),
  H_M_name      = clean(tax$H_M_name),
  stringsAsFactors = FALSE
)
clootl_year <- as.integer(clootl_year)

## --- legacy bridge: the original malaviR hand-curated key (MalAvi name ->
## corrected binomial). Kept only where the curator changed the name; used by
## match_taxonomy() as a last-resort lookup for names nothing else resolves. ----
legacy_env <- new.env()
load("data-raw/original_taxonomy_key.rda", envir = legacy_env)
oldkey <- legacy_env$taxonomy
oldkey$jetz <- gsub("_", " ", oldkey$Jetz.species)
keep <- !is.na(oldkey$jetz) & oldkey$jetz != oldkey$species & oldkey$species != ""
legacy_key <- stats::setNames(oldkey$jetz[keep], oldkey$species[keep])
legacy_key <- legacy_key[!duplicated(names(legacy_key))]

## --- manual overrides: the maintainer-curated key of MalAvi host names that
## have been hand-resolved to a current eBird species. Edit
## data-raw/manual_taxonomy.csv to add more; targets are checked against the
## bundled clootl taxonomy below. ----------------------------------------------
manual <- utils::read.csv("data-raw/manual_taxonomy.csv", stringsAsFactors = FALSE)
manual <- manual[!is.na(manual$ebird_species) & trimws(manual$ebird_species) != "", ]
manual_key <- stats::setNames(trimws(manual$ebird_species), trimws(manual$malavi_species))
manual_key <- manual_key[!duplicated(names(manual_key))]
bad <- manual_key[!manual_key %in% clootl_ref$SCI_NAME]
if (length(bad) > 0) {
  warning("manual_taxonomy.csv: ", length(bad),
          " target name(s) are not in the clootl taxonomy and will not match: ",
          paste(unique(bad), collapse = ", "))
}

## --- write internal data, then load the package so match_taxonomy() sees it ---
save(clootl_ref, clootl_year, legacy_key, manual_key,
     file = "R/sysdata.rda", compress = "xz")
message("Wrote R/sysdata.rda (clootl_ref: ", nrow(clootl_ref), " rows, year ",
        clootl_year, "; legacy_key: ", length(legacy_key),
        "; manual_key: ", length(manual_key), " entries)")

suppressMessages(devtools::load_all(quiet = TRUE))

## --- build the MalAvi host -> eBird crosswalk and save as exported data --------
res <- match_taxonomy()          # default: bundled MalAvi host species, latest release
taxonomy <- res$key
save(taxonomy, file = "data/taxonomy.rda", compress = "xz")

message("Wrote data/taxonomy.rda (", nrow(taxonomy), " host species)")
cat("\nmatch_type breakdown:\n")
print(table(taxonomy$match_type, useNA = "ifany"))

## --- audit: list the lowest-evidence matches for a maintainer to eyeball. The
## CSV is written next to this script as a local review aid (it is gitignored,
## not shipped); the same rows are printed below for convenience. ------------
audit <- malaviR:::.audit_taxonomy(taxonomy)
utils::write.csv(audit, "data-raw/taxonomy_audit.csv", row.names = FALSE)
cat("\nAudit: ", nrow(audit), " low-confidence rows to review ",
    "(weak genus reassignments + legacy bridges) written to ",
    "data-raw/taxonomy_audit.csv\n", sep = "")
print(audit, row.names = FALSE)
