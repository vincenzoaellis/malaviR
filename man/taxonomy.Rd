\name{taxonomy}
\alias{taxonomy}
\docType{data}
\title{
Host taxonomy in the MalAvi dataset
}
\description{
An object of class \code{data.frame} with the unique host species found in the MalAvi dataset (as of 14 September 2017) and their corresponding taxonomic classifications in the Jetz et al. phylogeny available on \url{birdtree.org}.
}
\usage{data("taxonomy")}
\format{
  A data frame with 1415 observations on the following 3 variables.
  \describe{
    \item{\code{species}}{a factor with levels equal to unique host species in the MalAvi dataset}
    \item{\code{Jetz.species}}{a factor with levels equal to the corresponding taxonomic classifications in the Jetz et al. dataset}
    \item{\code{match}}{a factor with levels \code{no} \code{yes}}
  }
}
\details{
The factor \code{Jetz.species} has underscores which correspond to how those names are stored as tip labels in the phylogenies available on \url{birdtree.org}.

Two host species in the MalAvi dataset do not have corresponding species in the Jetz et al. phylogeny:
\itemize{
  \item \emph{Anas sp.}
  \item \emph{Anthus cinnamomeus}
}

And 27 pairs and 2 triplets of species in the MalAvi dataset are recorded as single species in the Jetz et al. phylogeny:
\itemize{
  \item \emph{Acrocephalus_arundinaceus} \code{==} \emph{Acrocephalus arundinaceus / Acrocephalus orientalis}
  \item \emph{Agelastes_niger} \code{==} \emph{Agelastes niger / Phasidus niger}
  \item \emph{Andropadus_tephrolaemus} \code{==} \emph{Andropadus tephrolaemus / Arizelocichla tephrolaema}
  \item \emph{Aphelocoma_californica} \code{==} \emph{Aphelocoma californica / Aphelocoma woodhouseii}
  \item \emph{Arremon_torquatus} \code{==} \emph{Buarremon torquatus / Buarremon basilicus}
  \item \emph{Automolus_infuscatus} \code{==} \emph{Automolus infuscatus / Automolus paraensis}
  \item \emph{Basileuterus_coronatus} \code{==} \emph{Myiothlypis coronatus / Basileuterus coronatus / Myiothlypis coronata}
  \item \emph{Camaroptera_brachyura} \code{==} \emph{Camaroptera brachyura / Camaroptera brevicaudata}
  \item \emph{Cyanopica_cyanus} \code{==} \emph{Cyanopica cooki / Cyanopica cyanus}
  \item \emph{Ficedula_hypoleuca} \code{==} \emph{Ficedula hypoleuca / Ficedula speculigera}
  \item \emph{Hippolais_pallida} \code{==} \emph{Hippolais pallida / Iduna pallida}
  \item \emph{Icterus_cayanensis} \code{==} \emph{Icterus cayanensis / Icterus chrysocephalus}
  \item \emph{Larus_argentatus} \code{==} \emph{Larus argentatus / Larus mongolicus}
  \item \emph{Nectarinia_chloropygia} \code{==} \emph{Cinnyris chloropygius / Nectarinia chloropygia}
  \item \emph{Nectarinia_olivacea} \code{==} \emph{Cyanomitra olivacea / Nectarinia olivacea / Cyanomitra obscura}
  \item \emph{Nectarinia_oritis} \code{==} \emph{Cyanomitra oritis / Nectarinia oritis}
  \item \emph{Nectarinia_senegalensis} \code{==} \emph{Nectarinia senegalensis / Chalcomitra senegalensis}
  \item \emph{Nectarinia_venusta} \code{==} \emph{Nectarinia venusta / Cinnyris venustus}
  \item \emph{Petroica_australis} \code{==} \emph{Petroica australis / Petroica longipes}
  \item \emph{Phaeothlypis_fulvicauda} \code{==} \emph{Basileuterus fulvicauda / Phaeothlypis fulvicauda}
  \item \emph{Saxicola_torquatus} \code{==} \emph{Saxicola torquata / Saxicola maura}
  \item \emph{Serinus_mozambicus} \code{==} \emph{Serinus mozambicus / Crithagra mozambica}
  \item \emph{Sylvia_hortensis} \code{==} \emph{Sylvia hortensis / Sylvia crassirostris}
  \item \emph{Thraupis_sayaca} \code{==} \emph{Thraupis sayaca / Tangara sayaca}
  \item \emph{Thryothorus_genibarbis} \code{==} \emph{Pheugopedius genibarbis / Thryothorus genibarbis}
  \item \emph{Troglodytes_aedon} \code{==} \emph{Troglodytes musculus / Troglodytes aedon}
  \item \emph{Zosterops_borbonicus} \code{==} \emph{Zosterops borbonicus / Zosterops mauritianus}
  \item \emph{Zosterops_ficedulinus} \code{==} \emph{Zosterops ficedulinus / Zosterops feae}
  \item \emph{Zosterops_pallidus} \code{==} \emph{Zosterops capensis / Zosterops virens}
}

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
Bensch S, Hellgren O, and Perez-Tris J. 2009. MalAvi: a public database of malaria parasites and related
haemosporidians in avian hosts based on mitochondrial cytochrome \emph{b} lineages. Molecular Ecology Resources 9: 1353-1358.
\url{http://mbio-serv2.mbioekol.lu.se/Malavi/}

Jetz W, Thomas GH, Joy JB, Hartmann K, and Mooers AO. 2012. The global diversity of birds in space and time. Nature 491: 444-448.
\url{http://birdtree.org/}
}
\examples{
data(taxonomy)
}
\keyword{datasets}
