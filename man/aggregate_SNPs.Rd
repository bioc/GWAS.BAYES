\name{aggregate_SNPs}
\alias{aggregate_SNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{aggregate_SNPs
%%  ~~function to do ... ~~
}
\description{
Aggregate SNPs and Y by Species
}
\usage{
aggregate_SNPs(SNPs,Y,na.rm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{
Standardized SNP data set where the values of each column are either 0 or 1
}
\item{Y}{
The phenotype response of interest
}
\item{na.rm}{
Logical value where TRUE removes NA's in response vector as well corresponding rows in SNP matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{SNPs}{reduced SNP dataset}
\item{Y}{reduced Y vector}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data("vignette_lm_dat")
Y <- vignette_lm_dat$Phenotype
SNPs <- vignette_lm_dat[,-1]
SNPs <- standardize(SNPs = SNPs,method = "major-minor",number_cores = 1)

aggregate_SNPs(SNPs = SNPs, Y = Y)

}
