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
#' name). Many MalAvi host names are
#' older binomials that no longer match any of those because the genus has since
#' been split or the specific epithet re-gendered (e.g. \emph{Anas clypeata} is
#' now \emph{Spatula clypeata}; \emph{Basileuterus basilicus} is now
#' \emph{Myiothlypis basilica}). To recover these, a final step matches the
#' specific epithet -- allowing for Latin gender agreement -- within the host's
#' MalAvi family (or, if that family name is not used by clootl, within its
#' order), accepting the match only when it points to a single eBird species.
#' This resolves most genus reassignments while the family/order constraint
#' guards against epithet collisions between unrelated birds; names whose epithet
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
#'       \code{"reassigned:family"}, \code{"reassigned:order"}, \code{"legacy"},
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
    res <- .resolve_name(corrected, family[i], order[i], ref)
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

  ## 3. fall back to IOC / BirdLife / Howard & Moore synonyms, resolved to eBird name
  syn_sources <- c(IOC = "IOC_name", BirdLife = "Birdlife_name", HowardMoore = "H_M_name")
  for (label in names(syn_sources)) {
    todo <- which(is.na(ebird) & !generic)
    if (length(todo) == 0) break
    idx <- match(species[todo], ref[[syn_sources[label]]])
    got <- which(!is.na(idx))
    if (length(got) > 0) {
      pos <- todo[got]
      ebird[pos]      <- ref$SCI_NAME[idx[got]]
      match_type[pos] <- paste0("synonym:", label)
    }
  }

  ## 4. family/order-constrained epithet match (recovers genus reassignments and
  ##    gender-agreement changes); only accepted when it resolves to one species
  todo <- which(is.na(ebird) & !generic)
  for (i in todo) {
    res <- .epithet_reassign(species[i], family[i], order[i], ref)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- res$type
    }
  }

  ## 5. legacy bridge: the original malaviR hand-curated key maps some MalAvi
  ##    host names to a corrected binomial; re-resolve that to the eBird name
  todo <- which(is.na(ebird) & !generic)
  for (i in todo) {
    corrected <- legacy_key[species[i]]
    if (is.na(corrected)) next
    res <- .resolve_name(corrected, family[i], order[i], ref)
    if (!is.na(res$ebird)) {
      ebird[i]      <- res$ebird
      match_type[i] <- "legacy"
    }
  }

  match_type[is.na(ebird) & generic]  <- "generic"
  match_type[is.na(ebird) & !generic] <- "none"

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

  list(key = key, differences = key[key$match_type != "exact", ])
}
