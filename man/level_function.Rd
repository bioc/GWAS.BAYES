\name{level_function}
\alias{level_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{level_function
%%  ~~function to do ... ~~
}
\description{
Removes all SNPs that only have one level in it
}
\usage{
level_function(SNPs,MAF = 0.01)
}
\arguments{
  \item{SNPs}{
The standardized SNP data where the columns take on the values of 0 or 1
}
\item{MAF}{
The minor allele frequency at which to drop SNPs. Default is set to 0.01, meaning if the minor allele occurs less than 1 percent of the time in a given SNP, that given SNP will be dropped from the dataset.
}
}
\value{
\item{SNPs}{The SNP matrix where columns that were either all 1's or all 0's are removed}
\item{SNPs_Dropped}{A true/false vector with length ncol(SNPs), where the TRUE's indicate that the column was not dropped and the FALSE's indicate that the column was dropped}
}
\examples{
data("vignette_lm_dat")
Y <- vignette_lm_dat$Phenotype
SNPs <- vignette_lm_dat[,-1]
SNPs <- standardize(SNPs = SNPs,method = "major-minor",number_cores = 1)
list1 <- aggregate_SNPs(SNPs = SNPs, Y = Y)
SNPs <- list1[[1]]
Y <- list1[[2]]

level_function(SNPs, MAF = .01)
}
