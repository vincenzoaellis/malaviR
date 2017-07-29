\name{sister_taxa}
\alias{sister_taxa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify sister taxa in a phylogeny
}
\description{
The function identifies sister taxa (tips on the phylogeny) descending from a particular node in a phylogeny.
}
\usage{
sister_taxa(tree, node)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
an object of class \code{phylo}.
}
  \item{node}{
node number (see details).
}
}
\details{
The node numbers in an object of class \code{phylo} can be visualized by calling \code{nodelabels()} after plotting the phylogeny. The \code{sister_taxa()} function uses these numbers, so the first node (i.e., the most basal node) number will always be equal to the number of taxa plus one.
}
\value{
Returns an object of class \code{data.frame}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Vincenzo A. Ellis \email{vincenzoaellis@gmail.com}
}
\note{
The sister taxa are grouped by the variable \code{sister.clade} as either a \code{1} or a \code{2} to indicate the different sides of the node in question. The actual numerical values of \code{sister.clade} are arbitrary, they could have been called "A" and "B".
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## simulate a phylogenetic tree using the rtree function in the ape package
library(ape)
tree <- rtree(n=10)

## check out node numbers
plot(tree)
nodelabels()

## extract all sister taxa from the root
sister_taxa(tree, 11)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
