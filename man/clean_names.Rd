\name{clean_names}
\alias{clean_names}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Modify lineage names in sequence alignments to match the rest of the database
}
\description{
The sequence alignment files from MalAvi name the parasite lineages differently than the rest of the database. They add a letter indicating the parasite's genus to the beginning of the lineage name and they add the parasite's morphospecies (if one exists) to the end of the name. This function removes that extra information so that one can easily match lineage names from the sequence alignments to those in the rest of the database.
}
\usage{
clean_names(lin.names, keep.genus = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lin.names}{
The lineage names from the sequence alignment files to convert. Can be one or many.
}
  \item{keep.genus}{
If set to \code{TRUE} the genus of the lineage will be returned as well. Defaults to \code{FALSE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
If \code{keep.genus = FALSE} it will return an object (vector) of class \code{character}. If \code{keep.genus = TRUE} it will return an object of class \code{data.frame} with two columns (\code{parasiteGenus} and \code{Lineage_Name}).
}
\references{
%% ~put references to the literature/web site here ~
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
## download a sequence alignment file from MalAvi
seq.align <- extract_alignment("long seqs")

## lineage names are stored as rownames
lin.names <- rownames(seq.align)
head(lin.names)

## clean the names
lin.names.clean <- clean_names(lin.names)
head(lin.names.clean)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
