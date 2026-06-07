#!/usr/bin/env Rscript
## Install the new malaviR dependencies into the personal R library.
## DECIPHER + Biostrings (Bioconductor) for offline BLAST; clootl (CRAN) for taxonomy.
## Run on a compute node (heavy compile). Logged to logs/install_deps.log.

userlib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library", "4.3")
dir.create(userlib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(userlib, .libPaths()))

options(Ncpus = max(1L, parallel::detectCores()),
        repos = c(CRAN = "https://cloud.r-project.org"))

message("libPaths: ", paste(.libPaths(), collapse = " | "))
message("Bioconductor: ", as.character(BiocManager::version()))

inst <- rownames(installed.packages())

## Bioconductor: DECIPHER (pulls Biostrings et al.)
if (!"DECIPHER" %in% inst) {
  message("Installing DECIPHER (Bioconductor)...")
  BiocManager::install("DECIPHER", lib = userlib, update = FALSE, ask = FALSE)
}

## CRAN: clootl (modern avian taxonomy/phylogeny + crosswalk)
if (!"clootl" %in% inst) {
  message("Installing clootl (CRAN)...")
  install.packages("clootl", lib = userlib)
}

## Report
inst2 <- rownames(installed.packages())
for (pk in c("DECIPHER", "Biostrings", "clootl")) {
  message(sprintf("  %-12s installed: %s", pk, pk %in% inst2))
}
message("DONE install_deps.R")
