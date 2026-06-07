#!/usr/bin/env bash
## Create a conda env with R >= 4.4 + DECIPHER 3.x for building the MalAvi BLAST
## index (IndexSeqs) and testing blast_malavi(). System R is 4.3.3 (DECIPHER 2.30,
## which lacks IndexSeqs/SearchIndex/AlignPairs). Run on a compute node.
set -eEuo pipefail

module load miniforge
source /opt/miniforge3/etc/profile.d/conda.sh

ENV_PREFIX="$HOME/conda_envs/malaviR"

conda create -y -p "$ENV_PREFIX" \
  -c conda-forge -c bioconda \
  r-base=4.5 \
  bioconductor-decipher \
  bioconductor-biostrings \
  r-ape \
  r-readxl \
  r-magrittr

echo "=== env created at $ENV_PREFIX ==="
conda run -p "$ENV_PREFIX" Rscript -e 'cat("R:", R.version.string, "\n"); suppressMessages(library(DECIPHER)); cat("DECIPHER:", as.character(packageVersion("DECIPHER")), "\n"); for(f in c("IndexSeqs","SearchIndex","AlignPairs")) cat(" ", f, exists(f), "\n")'
echo "DONE make_r45_env.sh"
