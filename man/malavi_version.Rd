\name{malavi_version}
\alias{malavi_version}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Access current MalAvi version number
}
\description{
This function reads the current MalAvi version number and date of most recent update from the footer section of the MalAvi homepage.
}
\usage{
malavi_version()
}

\details{
The function does not take an argument. Simply call it and it will return the version info.
}
\value{
It returns a character string that can be found in the footer of the MalAvi homepage.
}
\references{
Bensch S, Hellgren O, and Perez-Tris J. 2009. MalAvi: a public database of malaria parasites and related
haemosporidians in avian hosts based on mitochondrial cytochrome \emph{b} lineages. Molecular Ecology Resources 9: 1353-1358.
\url{http://mbio-serv2.mbioekol.lu.se/Malavi/}
}
\author{
Vincenzo A. Ellis \email{vincenzoaellis@gmail.com}
}
\examples{
## just call the function without an argument
malavi_version()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
