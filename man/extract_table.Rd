\name{extract_table}
\alias{extract_table}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download MalAvi Tables
}
\description{
This function downloads the nine available MalAvi data tables.
}
\usage{
extract_table(table = "Hosts and Sites Table")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{table}{
One of nine table names in quotes or \code{"all"} (see details). Defaults to "Hosts and Sites Table".
}
}
\details{
The nine MalAvi tables are:
\itemize{
  \item "Hosts and Sites Table"
  \item "Table of References"
  \item "Grand Lineage Summary"
  \item "Parasite Summary Per Host"
  \item "Table of Lineage Names"
  \item "Morpho Species Summary"
  \item "Vector Data Table"
  \item "Other Genes Table"
  \item "Database Summary Report"
}

Alternatively, \code{"all"} will download all of the tables into a single list.
}
\value{
Each table will be in the form of an object of class \code{data.frame}; however \code{"all"} will result in a \code{list} with each table as an element.
}
\references{
Bensch S, Hellgren O, and Perez-Tris J. 2009. MalAvi: a public database of malaria parasites and related
haemosporidians in avian hosts based on mitochondrial cytochrome \emph{b} lineages. Molecular Ecology Resources 9: 1353-1358.
\url{http://mbio-serv2.mbioekol.lu.se/Malavi/}
}
\author{
Vincenzo A. Ellis \email{vincenzoaellis@gmail.com}
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## download the "Table of References"
refs <- extract_table("Table of References")
head(refs) # take a look
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
