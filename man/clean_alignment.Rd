\name{clean_alignment}
\alias{clean_alignment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify and remove repeated haplotypes from a (MalAvi) sequence alignment
}
\description{
Several lineages in the MalAvi database differ by ambiguous base calls only (e.g., "N" or "Y") and thus represent repeated haplotypes. For phylogenetic analysis it might make sense to only include one representative of any repeated haplotype because there is no way to know if they represent one or two lineages. This function identifies such repeated haplotypes in an alignment and randomly selects one of their lineages to be representative of the haplotype. Using this selection, the function subsets the alignment so that all haplotypes are only represented once.
}
\usage{
clean_alignment(alignment, separate_by_genus = FALSE, haplotype_format_wide = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alignment}{
a DNA sequence alignment of class \code{DNAbin}.
}
  \item{separate_by_genus}{
if the alignment is a MalAvi alignment with uncleaned sequence names (see details) you can choose to output the cleaned alignments by parasite genus by setting to \code{TRUE}. Defaults to \code{FALSE}.
}
  \item{haplotype_format_wide}{
if the lineage names associated with each repeated haplotype should be in wide format (\code{TRUE}, easier to visualize) or long format (\code{FALSE}, easier to subset). Defaults to \code{TRUE}
}
}
\details{
In a MalAvi alignment the default sequence (i.e., lineage) names have extra information and typically begin with a letter that indicates the parasite genus. This information can be used to separate the alignments by parasite genus if \code{separate_by_genus} is set to \code{TRUE}.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns a list composed of the following elements:
\item{repeated_haplotypes}{A data frame (in wide or long format) of repeated haplotypes and associated sequence (lineage) names}
\item{selected_lineages}{A vector of randomly selected sequence (lineage) names chosen to represent each repeated haplotype}
\item{alignment_clean}{A sequence alignment of class \code{DNAbin} that has only a single representative for each haplotype. Alternatively \code{alignment_clean_Plasmodium}, \code{alignment_clean_Haemoproteus}, \code{alignment_clean_Leucocytozoon} if \code{separate_by_genus} is set to \code{TRUE}}
}
\author{
Vincenzo A. Ellis \email{vincenzoaellis@gmail.com}
}
\examples{
## load the long seqs alignment from MalAvi then clean it
long.seqs <- extract_alignment("long seqs")
long.seqs.clean <- clean_alignment(long.seqs)
long.seqs.clean
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
