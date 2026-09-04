#' Match host species names to the clootl (eBird) avian taxonomy
#'
#' Aligns a set of bird species names to the avian taxonomy used by the
#' \pkg{clootl} package (the eBird/Clements taxonomy that underlies the
#' constantly updated avian phylogeny of McTavish et al. 2025). For each name it
#' returns the matching eBird species, the corresponding tip label in the clootl
#' phylogeny (\code{ott_name}), and the order and family, together with a
#' \code{match_type} describing how (or whether) it matched.
#'
#' Names are first looked up in a maintainer-curated override key
#' (\code{data-raw/manual_taxonomy.csv}) of MalAvi host names that have been
#' hand-resolved to a current eBird species; these are flagged
#' \code{match_type = "manual"}. Remaining names are matched against the eBird
#' scientific names, and, failing that, against the IOC, BirdLife, and Howard &
#' Moore synonyms carried by clootl (which are then resolved back to the eBird
#' name). Those synonym columns hold semicolon-joined lists wherever another
#' authority splits an eBird species, and every name in such a list is searched;
#' a host name found only inside one gets a \code{"-lump"} label (e.g.
#' \code{"synonym:BirdLife-lump"}), because eBird then treats that authority's
#' species as part of a broader one and the MalAvi host concept is narrower than
#' the name it maps to. Because the same synonym string is occasionally shared by
#' more than one eBird species, an ambiguous synonym is accepted only for the
#' candidate whose own epithet agrees with the host's, never by silently taking
#' the first.
#' Many MalAvi host names are
#' older binomials that no longer match any of those because the genus has since
#' been split or the specific epithet re-gendered (e.g. \emph{Anas clypeata} is
#' now \emph{Spatula clypeata}; \emph{Basileuterus basilicus} is now
#' \emph{Myiothlypis basilica}). To recover these, the specific epithet is
#' matched -- allowing for Latin gender agreement -- first within the same genus
#' (e.g. \emph{Saxicola maura} to \emph{Saxicola maurus}), and then, for genuine
#' genus transfers, within the host's
#' MalAvi family (or, if that family name is not used by clootl, within its
#' order), accepting the match only when it points to a single eBird species.
#' The same-genus step is tried first because a fixed genus is the strongest
#' identity signal and is not misled by a mislabeled MalAvi family; the
#' family/order constraint then
#' guards against epithet collisions between unrelated birds. That constraint is
#' only as good as MalAvi's family label, which is the least maintained field in
#' the release, so a genus-changing match is additionally checked against the
#' families clootl files the MalAvi genus in and declined when it falls outside
#' them. Names whose epithet
#' remains ambiguous are left unmatched rather than guessed. As a last step, host
#' names still unmatched are looked up in the hand-curated species key from the
#' original \code{malaviR} (which mapped many MalAvi names to corrected
#' binomials); the corrected name is then resolved to the current eBird name and
#' flagged \code{match_type = "legacy"}. These legacy matches come from a
#' hand-curated key made years ago against the Jetz \emph{et al.} (BirdTree)
#' taxonomy and may reflect taxonomic decisions that are now out of date, so they
#' are worth double-checking.
#'
#' Leading/trailing whitespace is removed before matching. Some MalAvi host names
#' are not identifiable binomials -- entries ending in \dQuote{sp.}, hybrids
#' written with \dQuote{ x }, or bare genus names -- and can never match; these
#' are flagged \code{match_type = "generic"} rather than forced to a species.
#'
#' The clootl taxonomy is bundled with \code{malaviR} as a dated snapshot, so no
#' internet connection or \pkg{clootl} installation is needed at run time. See
#' \code{\link{clootl_taxonomy_version}} for the bundled taxonomy year.
#'
#' @section Interpreting \code{match_type}:
#' The values differ in how much they should be trusted:
#' \describe{
#'   \item{\code{"exact"}}{the strongest match: the host name is itself a current
#'     eBird scientific name. It is a match of \emph{names}, not of species
#'     concepts: MalAvi records accrue under the name the original study used, so
#'     an old broad name that is still current for a narrower species matches
#'     exactly and silently. Of the 19 MalAvi records for \emph{Tyto alba}, 12 are
#'     from the New World or Australasia and belong to \emph{T. furcata} or
#'     \emph{T. javanica} under the bundled taxonomy. Exact matches need no
#'     taxonomic review, but records under a recently split name still need
#'     locality checking.}
#'   \item{\code{"manual"}, \code{"synonym:*"}, \code{"reassigned:*"},
#'     \code{"legacy"}}{resolved, but by a rule rather than an exact hit -- a
#'     maintainer override, a recognized synonym, an epithet/genus reassignment, or
#'     the old \code{malaviR} key. \strong{Inspect these} when the exact taxonomy
#'     matters for your analysis; \code{"legacy"} and weak reassignments are review
#'     targets, not necessarily errors.}
#'   \item{\code{"generic"}, \code{"none"}}{\strong{not} resolved to a species:
#'     \code{"generic"} names can never match (\dQuote{sp.}, hybrids, bare genera)
#'     and \code{"none"} simply did not match. Treat these as unresolved.}
#' }
#' The \code{differences} table returned below is exactly the set of non-\code{exact}
#' rows, i.e. the ones worth checking by hand.
#'
#' @param species Character vector of species names to match. If \code{NULL}
#'   (default), the unique host species in the bundled MalAvi
#'   \code{"Hosts and Sites Table"} are used, along with their MalAvi family and
#'   order.
#' @param version MalAvi release to take host names from when \code{species} is
#'   \code{NULL}; a date string or \code{"latest"} (default).
#' @param family,order Optional character vectors, the same length as
#'   \code{species}, giving each name's family and order. They are used only for
#'   the family/order-constrained epithet step (see Details) and are taken from
#'   MalAvi automatically when \code{species} is \code{NULL}. If you supply your
#'   own \code{species} without them, that recovery step is simply skipped.
#' @return A list with two data frames:
#'   \describe{
#'     \item{\code{key}}{one row per input species, with columns
#'       \code{malavi_species}, \code{ebird_species}, \code{ott_name},
#'       \code{order}, \code{family}, and \code{match_type} (one of
#'       \code{"manual"}, \code{"exact"}, \code{"synonym:IOC"},
#'       \code{"synonym:BirdLife"}, \code{"synonym:HowardMoore"},
#'       \code{"reassigned:genus"}, \code{"reassigned:family"},
#'       \code{"reassigned:order"}, \code{"legacy"},
#'       \code{"generic"}, or \code{"none"}).}
#'     \item{\code{differences}}{the subset of \code{key} that did not match an
#'       eBird name exactly (manual overrides, synonyms, reassignments, legacy
#'       matches, generics, and unmatched names) -- the rows worth checking by
#'       hand.}
#'   }
#' @references
#' McTavish EJ, Gerbracht JA, Holder MT, Iliff MJ, Lepage D, Rasmussen PC,
#' Redelings BD, Sanchez Reyes LL, Miller ET (2025). A complete and dynamic tree
#' of birds. Proceedings of the National Academy of Sciences 122(18):
#' e2409658122. \doi{10.1073/pnas.2409658122}
#' @seealso \code{\link{taxonomy}} for the pre-built crosswalk of MalAvi hosts,
#'   \code{\link{clootl_taxonomy_version}}
#' @examples
#' res <- match_taxonomy(c("Turdus merula", "Cyanistes caeruleus", "Anas sp."))
#' res$key
#' res$differences
#' @export
match_taxonomy <- function(species = NULL, version = "latest",
                           family = NULL, order = NULL) {
  ## record whether host names came from the bundled data *before* reassigning
  ## `species`, so the version stamp below is correct
  used_bundled <- is.null(species)
  if (is.null(species)) {
    hosts   <- extract_table("Hosts and Sites Table", version = version)
    info    <- .host_family_order(hosts)
    species <- info$species
    family  <- info$family
    order   <- info$order
  } else {
    species <- trimws(species)
    if (!is.null(family) && length(family) != length(species))
      stop("'family' must be the same length as 'species'.", call. = FALSE)
    if (!is.null(order) && length(order) != length(species))
      stop("'order' must be the same length as 'species'.", call. = FALSE)
    keep    <- !is.na(species) & species != ""
    species <- species[keep]
    family  <- if (is.null(family)) rep(NA_character_, length(species)) else trimws(family)[keep]
    order   <- if (is.null(order))  rep(NA_character_, length(species)) else trimws(order)[keep]
    dup     <- !duplicated(species)
    species <- species[dup]; family <- family[dup]; order <- order[dup]
  }

  ref <- clootl_ref  # bundled clootl taxonomy snapshot (internal data)
  ref$latin_family <- sub(" .*$", "", ref$FAMILY)  # clootl FAMILY is "Anatidae (Ducks...)"

  ## clootl's synonym columns, exploded to one name per row (they hold
  ## semicolon-joined lists wherever another authority splits an eBird species).
  ## Built once here and passed down, rather than rebuilt for every host name.
  syn_table <- .syn_table(ref)

  ## names that can never match a single species (".../ sp.", "... spp", hybrids,
  ## or bare genus names)
  generic <- grepl(" spp?\\.?$", species) | grepl(" x ", species) | !grepl(" ", species)

  ebird      <- rep(NA_character_, length(species))
  match_type <- rep(NA_character_, length(species))

  ## 1. manual overrides: a maintainer-curated key
  ##    (data-raw/manual_taxonomy.csv) of MalAvi host names hand-resolved to a
  ##    current eBird species. Applied first because it is authoritative.
  for (i in which(!generic)) {
    corrected <- manual_key[species[i]]
    if (is.na(corrected)) next
    res <- .resolve_name(corrected, family[i], order[i], ref, syn_table)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- "manual"
    }
  }

  ## 2. exact match to the eBird scientific name (names not already set above)
  hit  <- match(species, ref$SCI_NAME)
  fill <- which(is.na(ebird) & !is.na(hit))
  ebird[fill]      <- ref$SCI_NAME[hit[fill]]
  match_type[fill] <- "exact"

  ## 3. fall back to IOC / BirdLife / Howard & Moore synonyms, resolved to eBird
  ##    name. .syn_resolve() is ambiguity-aware: when the same synonym string is
  ##    carried by more than one eBird species it keeps only the one whose epithet
  ##    agrees, rather than silently taking the first row.
  todo <- which(is.na(ebird) & !generic)
  for (i in todo) {
    res <- .syn_resolve(species[i], ref, syn_table)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- res$type
    }
  }

  ## 4. same-genus epithet shift (genus unchanged, epithet re-gendered, e.g.
  ##    Saxicola maura -> Saxicola maurus). Tried before the family/order step
  ##    because a fixed genus is the strongest identity signal and is not fooled
  ##    by a wrong MalAvi family label.
  todo <- which(is.na(ebird) & !generic)
  for (i in todo) {
    res <- .same_genus_reassign(species[i], ref)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- res$type
    }
  }

  ## 5. family/order-constrained epithet match (recovers genus reassignments and
  ##    gender-agreement changes); only accepted when it resolves to one species
  todo <- which(is.na(ebird) & !generic)

  for (i in todo) {
    res <- .epithet_reassign(species[i], family[i], order[i], ref)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- res$type
    }
  }

  ## 6. legacy bridge: the original malaviR hand-curated key maps some MalAvi
  ##    host names to a corrected binomial; re-resolve that to the eBird name
  todo <- which(is.na(ebird) & !generic)
  for (i in todo) {
    corrected <- legacy_key[species[i]]
    if (is.na(corrected)) next
    res <- .resolve_name(corrected, family[i], order[i], ref, syn_table)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- "legacy"
    }
  }

  match_type[is.na(ebird) & generic]  <- "generic"
  match_type[is.na(ebird) & !generic] <- "none"

  ## The family/order-constrained epithet match (step 5) is the one step that needs
  ## information about the host beyond its name: it constrains its candidate pool by the
  ## host's family, falling back to its order, so without them it cannot run at all -- and
  ## it used to skip in silence. `match_taxonomy()` fetches both from MalAvi when `species`
  ## is NULL; a caller passing their own vector has to supply them too. Called the second
  ## way, 173 of MalAvi's own 2,339 host binomials come back "none" that this very function
  ## resolves when called the first way, and nothing explained the difference. That was
  ## read once as a missing rule in the function rather than a missing argument in the
  ## call, so it is worth saying out loud.
  ##
  ## Checked here rather than at step 5 deliberately: a name still unresolved at step 5 may
  ## yet be recovered by the legacy bridge at step 6, and warning earlier fired on names
  ## that resolved perfectly well.
  n_none <- sum(match_type == "none")
  if (n_none > 0 && all(is.na(family)) && all(is.na(order))) {
    warning(n_none, " name(s) did not resolve, and no `family`/`order` was supplied -- ",
            "the family/order-constrained epithet match cannot run without them and was ",
            "skipped, so genus reassignments (e.g. Grus leucogeranus -> Leucogeranus ",
            "leucogeranus) are reported as unmatched. Either pass family = and order = ",
            "alongside `species`, or call match_taxonomy() with no arguments to use ",
            "MalAvi's own host list, which carries them. For MalAvi host names the ",
            "resolved crosswalk is already shipped as `malaviR::taxonomy`.",
            call. = FALSE)
  }

  n_legacy <- sum(match_type == "legacy")
  if (n_legacy > 0) {
    message(n_legacy, " name(s) matched via the original malaviR hand-curated key ",
            "(match_type \"legacy\"). These are old, possibly out-of-date choices ",
            "-- please double-check them (see the 'differences' table).")
  }

  ## attach phylogeny tip label, order, family from the matched eBird species
  ridx <- match(ebird, ref$SCI_NAME)
  key <- data.frame(
    malavi_species = species,
    ebird_species  = ebird,
    ott_name       = ref$ott_name[ridx],
    order          = ref$ORDER1[ridx],
    family         = ref$FAMILY[ridx],
    match_type     = match_type,
    stringsAsFactors = FALSE
  )
  key <- key[order(key$malavi_species), ]
  rownames(key) <- NULL

  out <- list(key = key, differences = key[key$match_type != "exact", ])
  ## MalAvi version is only meaningful when host names came from the bundled data
  mv <- if (used_bundled) .malavi_resolve_version(version) else NA_character_
  out <- .malavi_attach_meta(out, malavi_version = mv,
                             clootl_version = clootl_year)
  class(out) <- c("malavi_taxonomy_match", class(out))
  out
}

#' @export
print.malavi_taxonomy_match <- function(x, ...) {
  key <- x$key
  cat("MalAvi -> clootl taxonomy match\n")
  meta <- .malavi_meta_line(x)
  if (!is.null(meta)) cat("  ", meta, "\n", sep = "")
  cat("  species:  ", nrow(key), "\n", sep = "")
  ## roll the detailed synonym:/reassigned: labels up to their family for a
  ## compact summary, but leave $key/$differences untouched
  fam <- sub(":.*$", "", key$match_type)
  tab <- sort(table(fam), decreasing = TRUE)
  cat("  by match_type:\n")
  for (nm in names(tab)) cat("    ", nm, ": ", tab[[nm]], "\n", sep = "")
  cat("\nNote: only match_type \"exact\" is a current-name hit; inspect the rest\n",
      "(see $differences) when taxonomy matters. \"generic\"/\"none\" are unresolved.\n",
      sep = "")
  invisible(x)
}
