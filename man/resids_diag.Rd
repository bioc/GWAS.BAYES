\name{resids_diag}
\alias{resids_diag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{resids_diag
%%  ~~function to do ... ~~
}
\description{
Tests residuals to see if they are normal. This looks at the model with all significant SNPs from the preselection phase.
}
\usage{
resids_diag(Y,SNPs,significant,kinship = FALSE,principal_components = FALSE,plot_it = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{Y}{
The phenotype response of interest
}
\item{SNPs}{
Standardized SNP data set where the values of each column are either 0 or 1
}
\item{significant}{
A vector of 0's and 1's where the 1's indicate a significant SNP. This is returned in the output of the preselection function.
}
\item{kinship}{
A kinship matrix, can be calculated from the rrBLUP package.
}
\item{principal_components}{
A matrix or vector of the principal components one would like to include in the analysis.
}
\item{plot_it}{
If TRUE a histogram of the residuals is returned.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{value 1}{The output of a Shapiro-Wilk test for the residuals. If the p-value is above .05, there is no evidence that the residuals are not normal. If the p-value is below .05 there is evidence that the residuals are not normal, and some transformation is suggested.}
\item{value 2}{A histogram of the residuals when plot_it = TRUE}
}
\examples{
data("vignette_lm_dat")
Y <- vignette_lm_dat$Phenotype
SNPs <- vignette_lm_dat[,-1]
fullPreprocess <- preprocess_SNPs(SNPs = SNPs,Y = Y,MAF = 0.01,number_cores = 1)
SNPs <- fullPreprocess$SNPs
Y <- fullPreprocess$Y
fullPreprocess$SNPs_Dropped
principal_comp <- pca_function(SNPs = SNPs,number_components = 1,plot_it = FALSE)
Significant_SNPs <- preselection(Y = Y, SNPs = SNPs,number_cores = 1, principal_components = principal_comp,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,info = FALSE)

resids_diag(Y = Y,SNPs = SNPs,significant = Significant_SNPs$Significant,kinship = FALSE,principal_components = principal_comp,plot_it = TRUE)
}
