## Evaluate a two-parent breakpoint chimera detector against the current screen.
##
## Negatives: real MalAvi lineages, leave-one-out (the lineage removed from the reference).
##            None is a chimera, so any call is a false positive.
## Positives: synthetic chimeras -- two real lineages spliced at a breakpoint, with BOTH
##            parents left in the reference. That is the realistic case: a PCR chimera is
##            formed from templates that exist in the data.
##
## Run: Rscript data-raw/chimera_v2_eval.R
suppressMessages(devtools::load_all("."))
set.seed(20260820)

aln  <- extract_alignment()
nm   <- rownames(aln)
seqs <- toupper(vapply(seq_len(nrow(aln)),
                       function(i) paste(as.character(aln[i, ]), collapse = ""),
                       character(1)))
names(seqs) <- nm
clean <- nm[!grepl("[^ACGT]", seqs) & nchar(seqs) == 479L]
cat(sprintf("alignment %d x %d; %d clean full-length lineages\n",
            nrow(aln), ncol(aln), length(clean)))

code_vec <- function(s) malaviR:::.qc_code_vec(strsplit(s, "", fixed = TRUE)[[1]])
refcode_all <- do.call(rbind, lapply(seqs, code_vec))
rownames(refcode_all) <- nm

## ---- the candidate detector -------------------------------------------------------
## For every breakpoint b, the best two-parent explanation is
##   min_A mismatches(query[1..b], A) + min_B mismatches(query[b+1..L], B),
## computed for all references at once with cumulative sums. Delta is how much better
## that is than the best single parent. Guards that the old screen lacked:
##   * the two parents must be different lineages;
##   * each segment must be at least `min_segment` long, so a breakpoint a few bases from
##     an end cannot be manufactured out of trimming noise;
##   * the two parents must themselves differ over the query, by `min_parent_dist` --
##     "switching" between two lineages 1 bp apart is not evidence of anything;
##   * references must cover most of the query, so a thinly covered one cannot win by
##     having little to disagree over.
detect_v2 <- function(qcode, refcode, ref_names,
                      min_segment = 60L, min_parent_dist = 5L, min_cover = 0.9) {
  L <- length(qcode); n <- nrow(refcode)
  qmat  <- matrix(qcode, n, L, byrow = TRUE)
  known <- (refcode > 0L) & (qmat > 0L)
  mism  <- known & (refcode != qmat)

  keep <- rowSums(known) >= min_cover * sum(qcode > 0L)
  if (sum(keep) < 2L) return(NULL)
  mism <- mism[keep, , drop = FALSE]; known <- known[keep, , drop = FALSE]
  rn   <- ref_names[keep]; refc <- refcode[keep, , drop = FALSE]

  Pm <- t(apply(mism, 1L, cumsum))                 # mismatches in 1..b
  tot <- Pm[, L]
  d1  <- min(tot); a1 <- rn[which.min(tot)]

  bs <- seq.int(min_segment, L - min_segment)
  pre  <- Pm[, bs, drop = FALSE]                   # mismatches left of each breakpoint
  post <- tot - pre                                # and right of it
  best_pre_i  <- apply(pre,  2L, which.min)
  best_post_i <- apply(post, 2L, which.min)
  d2 <- pre[cbind(best_pre_i, seq_along(bs))] + post[cbind(best_post_i, seq_along(bs))]

  ok <- best_pre_i != best_post_i                  # two parents, not one
  if (!any(ok)) return(list(delta = 0, d1 = d1, best_single = a1))
  d2[!ok] <- Inf
  k <- which.min(d2)
  A <- best_pre_i[k]; B <- best_post_i[k]

  ## the parents must be distinguishable over the positions the query defines
  both <- refc[A, ] > 0L & refc[B, ] > 0L & qcode > 0L
  parent_dist <- sum(refc[A, both] != refc[B, both])

  list(delta = d1 - d2[k], d1 = d1, d2 = d2[k],
       breakpoint = bs[k], parent_a = rn[A], parent_b = rn[B],
       parent_dist = parent_dist,
       passes_parent_dist = parent_dist >= min_parent_dist,
       best_single = a1)
}

## ---- build the two labeled sets --------------------------------------------------
N <- 200L
neg_names <- sample(clean, N)
pos_pairs <- data.frame(a = sample(clean, N), b = sample(clean, N),
                        stringsAsFactors = FALSE)
pos_pairs <- pos_pairs[pos_pairs$a != pos_pairs$b, ]

splice <- function(a, b, bp) paste0(substr(seqs[[a]], 1, bp),
                                    substr(seqs[[b]], bp + 1L, 479L))

run <- function(qseq, drop = character(0)) {
  keep <- setdiff(nm, drop)
  detect_v2(code_vec(qseq), refcode_all[keep, , drop = FALSE], keep)
}

cat("\nscoring", length(neg_names), "negatives (leave-one-out)...\n")
neg <- lapply(neg_names, function(x) run(seqs[[x]], drop = x))

cat("scoring", nrow(pos_pairs), "synthetic chimeras (parents left in)...\n")
pos <- lapply(seq_len(nrow(pos_pairs)), function(i) {
  bp <- sample(120:360, 1L)
  run(splice(pos_pairs$a[i], pos_pairs$b[i], bp))
})

getd  <- function(x) if (is.null(x) || is.null(x$delta)) NA_real_ else x$delta
getpd <- function(x) if (is.null(x) || is.null(x$parent_dist)) NA_real_ else x$parent_dist
nd <- vapply(neg, getd, 0); pd <- vapply(pos, getd, 0)
npd <- vapply(neg, getpd, 0); ppd <- vapply(pos, getpd, 0)

cat("\ndelta on NEGATIVES: median", median(nd, na.rm = TRUE),
    " 95th", quantile(nd, .95, na.rm = TRUE), " max", max(nd, na.rm = TRUE), "\n")
cat("delta on POSITIVES: median", median(pd, na.rm = TRUE),
    " 5th",  quantile(pd, .05, na.rm = TRUE), " min", min(pd, na.rm = TRUE), "\n")

cat("\n threshold |  sensitivity  |  false positives\n")
cat("-----------+---------------+------------------\n")
for (th in c(3, 5, 8, 10, 12, 15, 20, 25, 30)) {
  sens <- mean(pd >= th & ppd >= 5, na.rm = TRUE)
  fpr  <- mean(nd >= th & npd >= 5, na.rm = TRUE)
  cat(sprintf("%10d | %6.1f%%       | %6.1f%%\n", th, 100 * sens, 100 * fpr))
}
