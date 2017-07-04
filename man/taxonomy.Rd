\name{taxonomy}
\alias{taxonomy}
\docType{data}
\title{
Host taxonomy in the MalAvi dataset
}
\description{
An object of class \code{data.frame} with the unique host species found in the MalAvi dataset (as of 08 June 2017) and their corresponding taxonomic classifications in the Jetz et al. phylogeny available on \url{birdtree.org}.
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
This is useful information if analyses with the MalAvi dataset require integrating host phylogenies.

The factor \code{Jetz.species} has underscores which correspond to how those names are stored as tip labels in the phylogenies available on \url{birdtree.org}.
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
