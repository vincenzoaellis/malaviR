\name{extract_alignment}
\alias{extract_alignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Download MalAvi Sequence Alignments
}
\description{
This function downloads the three sequence alignment files available on the MalAvi website.
}
\usage{
extract_alignment(alignment = "all seqs")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{
One of the three sequence alignments (see details). Defaults to \code{"all seqs"} which is the complete alignment file.
}
}
\details{
This is a wrapper for the \code{read.dna} function from the \code{ape} package.

The three MalAvi alignments that can be called are:
\itemize{
  \item \code{"all seqs"} which includes all of the sequences in MalAvi
  \item \code{"long seqs"} which includes only sequences with all or nearly all 479bp
  \item \code{"morpho seqs"} which includes only sequences associated with a morphological species
}
}
\value{
Each alignment is returned in the form of an object of class \code{DNAbin}.
}
\references{
Bensch S, Hellgren O, and Perez-Tris J. 2009. MalAvi: a public database of malaria parasites and related
haemosporidians in avian hosts based on mitochondrial cytochrome \emph{b} lineages. Molecular Ecology Resources 9: 1353-1358.
\url{http://mbio-serv2.mbioekol.lu.se/Malavi/}

Paradis E., Claude J. & Strimmer K. 2004. APE: analyses of phylogenetics and evolution in R language. Bioinformatics 20: 289-290.
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
## download the "all seqs" alignment
all.seqs <- extract_alignment("all seqs")
all.seqs # take a look
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
