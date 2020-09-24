\name{preprocess_SNPs}
\alias{preprocess_SNPs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{preprocess_SNPs
%%  ~~function to do ... ~~
}
\description{
This functions takes raw SNP data and the associated phenotype response and returns a SNP dataset and phenotype response variable that can be used in the preselection function.
}
\usage{
preprocess_SNPs(SNPs,Y,MAF = 0.01,number_cores,na.rm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{
SNP data where each column is a SNP and the SNP column takes on the values A, C, T, or G.
}
\item{Y}{
The phenotype response of interest. Should be a numeric vector.
}
\item{MAF}{
The minor allele frequency at which to drop SNPs. Default is set to 0.01, meaning if the minor allele occurs less than 1 percent of the time in a given SNP, that given SNP will be dropped from the dataset.
}
\item{number_cores}{
The number of cores one would wish to parallelize over.
}
\item{na.rm}{
If there is NA's in the vector Y, set na.rm = TRUE and the Y values that are NA will be removed as well as the corresponding rows of the SNP matrix.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{SNPs}{A new SNP matrix. The matrix will be formatted so the minor allele's are coded as 0's and the major allele's are coded as 1's. This matrix will have columns dropped that have minor allele frequency less than the specified value. It will also aggregate over replications, so SNP's and the vector Y will be aggregated according to replications in the SNP matrix.}
\item{Y}{The new aggregated response vector Y. If you did not have any replications then this vector will be the exact same as the one entered.}
\item{SNPs_Dropped}{This will tell you which SNPs were dropped if the had minor allele frequency less than the specified value, it will be in the form of column index number. If no SNPs were dropped this will be the character string "None".}
}
\examples{
data("vignette_lm_dat")
Y <- vignette_lm_dat$Phenotype
SNPs <- vignette_lm_dat[,-1]

preprocess_SNPs(SNPs = SNPs,Y = Y,MAF = 0.01,number_cores = 1,na.rm = FALSE)
}
