#' Match host species names to the clootl (eBird) avian taxonomy
#'
#' Aligns a set of host species names to the modern avian taxonomy used by the
#' \pkg{clootl} package (the eBird/Clements taxonomy that underlies the
#' constantly updated avian phylogeny of McTavish et al.). For each name it
#' returns the matching eBird species, the corresponding tip label in the clootl
#' phylogeny (\code{ott_name}), and the order and family, together with a
#' \code{match_type} describing how (or whether) it matched.
#'
#' Names are matched first against the eBird scientific names, and, failing that,
#' against the IOC, BirdLife, and Howard & Moore synonyms carried by clootl
#' (which are then resolved back to the eBird name). Leading/trailing whitespace
#' is removed before matching. Some MalAvi host names are not identifiable
#' binomials -- entries ending in \dQuote{sp.}, hybrids written with \dQuote{ x },
#' or bare genus names -- and can never match; these are flagged
#' \code{match_type = "generic"} rather than forced to a species.
#'
#' The clootl taxonomy is bundled with \code{malaviR} as a dated snapshot, so no
#' internet connection or \pkg{clootl} installation is needed at run time. See
#' \code{\link{clootl_taxonomy_version}} for the bundled taxonomy year.
#'
#' @param species Character vector of species names to match. If \code{NULL}
#'   (default), the unique host species in the bundled MalAvi
#'   \code{"Hosts and Sites Table"} are used.
#' @param version MalAvi release to take host names from when \code{species} is
#'   \code{NULL}; a date string or \code{"latest"} (default).
#' @return A list with two data frames:
#'   \describe{
#'     \item{\code{key}}{one row per input species, with columns
#'       \code{malavi_species}, \code{ebird_species}, \code{ott_name},
#'       \code{order}, \code{family}, and \code{match_type} (one of
#'       \code{"exact"}, \code{"synonym:IOC"}, \code{"synonym:BirdLife"},
#'       \code{"synonym:HowardMoore"}, \code{"generic"}, or \code{"none"}).}
#'     \item{\code{differences}}{the subset of \code{key} that did not match an
#'       eBird name exactly (synonyms, generics, and unmatched names) -- the rows
#'       worth checking by hand.}
#'   }
#' @seealso \code{\link{taxonomy}} for the pre-built crosswalk of MalAvi hosts,
#'   \code{\link{clootl_taxonomy_version}}
#' @examples
#' res <- match_taxonomy(c("Turdus merula", "Cyanistes caeruleus", "Anas sp."))
#' res$key
#' res$differences
#' @export
match_taxonomy <- function(species = NULL, version = "latest") {
  if (is.null(species)) {
    hosts <- extract_table("Hosts and Sites Table", version = version)
    species <- hosts$SPECIES_NAME
  }
  species <- trimws(species)
  species <- unique(species[!is.na(species) & species != ""])

  ref <- clootl_ref  # bundled clootl taxonomy snapshot (internal data)

  ## names that can never match a single species
  generic <- grepl(" sp\\.?$", species) | grepl(" x ", species) | !grepl(" ", species)

  ebird      <- rep(NA_character_, length(species))
  match_type <- rep(NA_character_, length(species))

  ## exact match to the eBird scientific name
  hit <- match(species, ref$SCI_NAME)
  ebird[!is.na(hit)]      <- ref$SCI_NAME[hit[!is.na(hit)]]
  match_type[!is.na(hit)] <- "exact"

  ## fall back to IOC / BirdLife / Howard & Moore synonyms, resolved to eBird name
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

  match_type[is.na(ebird) & generic]  <- "generic"
  match_type[is.na(ebird) & !generic] <- "none"

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
