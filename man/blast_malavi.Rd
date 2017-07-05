\name{blast_malavi}
\alias{blast_malavi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BLAST MalAvi
}
\description{
A function for comparing a cytochrome \emph{b} DNA sequence with those in the MalAvi database.
}
\usage{
blast_malavi(sequence, evalue = 1e-80, hits = 5, print.alignments = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sequence}{
Sequence to BLAST.
}
  \item{evalue}{
Threshold for finding real hits. Should be a low number. Defaults to 1e-80.
}
  \item{hits}{
Number of hits to return. Defaults to 5.
}
  \item{print.alignments}{
Should the pairwise alignments be printed? Defaults to FALSE.
}
}
\details{
The function submits the sequence to the BLAST program on the MalAvi website and returns the result in a useful format.

\strong{hits} are "matching" sequences.
}
\value{
An object of class \code{data.frame} with resulting hits and information regarding how well they match the input sequence. If no hits were found a message will print asking you to check your input sequence.
}
\references{
Bensch S, Hellgren O, and Perez-Tris J. 2009. MalAvi: a public database of malaria parasites and related
haemosporidians in avian hosts based on mitochondrial cytochrome \emph{b} lineages. Molecular Ecology Resources 9: 1353-1358.
\cr

\url{http://mbio-serv2.mbioekol.lu.se/Malavi/blast.html}
}
\author{
Vincenzo A. Ellis \email{vincenzoaellis@gmail.com}
}
\note{
\code{Perfect.Match} in the output data frame refers to whether the pairwise alignments were identical irrespective of the coverage (e.g., an input sequence with low coverage could still generate, and often will generate, a perfect match with one or more sequences in MalAvi).
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## define a sequence. This is the Plasmodium parasite ACAGR1
ACAGR1 <- "GCAACTGGTGCTTCATTTGTATTTATTTTAACTTATTTACATATTTTAAGAGGATTAAATTATTCATATTCATATTTACCTTTATCATGGATATCTGGATTAATAATATTTTTAATATCTATAGTAACAGCTTTTATGGGTTACGTATTACCTTGGGGTCAAATGAGTTTCTGGGGTGCTACCGTAATAACTAATTTATTATATTTTATACCTGGACTAGTTTCATGGATATGTGGTGGATATCTTGTAAGTGACCCAACCTTAAAAAGATTCTTTGTACTACATTTTACATTTCCTTTTATAGCTTTATGTATTGTATTTATACATATATTCTTTCTACATTTACAAGGTAGCACAAATCCTTTAGGGTATGATACAGCTTTAAAAATACCCTTCTATCCAAATCTTTTAAGTCTTGATATTAAAGGATTTAATAATGTATTAGTATTATTTTTAGCACAAAGTTTATTTGGAATACT"
blast_malavi(ACAGR1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
