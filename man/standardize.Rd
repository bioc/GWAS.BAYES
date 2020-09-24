\name{standardize}
\alias{standardize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{standardize
%%  ~~function to do ... ~~
}
\description{
Standardize the SNPs to the 0-1 scale
}
\usage{
standardize(SNPs,method=c("major-minor","alphabetical"),number_cores)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{
The SNP dataset with columns of the values A, C, T, and G
}
\item{method}{
The method in which to standarize. If "major-minor" is selected then the major allele gets the value 1 and the minor allele gets the value 0. If "alphabetical" is selected the first allele alphabetically gets a value of 0 and the second allele alphabetically gets the value of 1.
}
\item{number_cores}{
The number of cores on which to parallize over.
}
}
\value{
Returns a matrix of SNPs with the same dimension as th input, the columns are now defined as 0 and 1's.
}
\examples{
data("vignette_lm_dat")
Y <- vignette_lm_dat$Phenotype
SNPs <- vignette_lm_dat[,-1]

standardize(SNPs = SNPs,method = "major-minor",number_cores = 1)

}
