\name{postGWAS}
\alias{postGWAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{postGWAS
%%  ~~function to do ... ~~
}
\description{
Performs GA model selection to identify the best model
}
\usage{
postGWAS(Y,SNPs,significant,number_cores,principal_components,maxiterations,runs_til_stop,kinship = FALSE,info = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The phenotype response on the reduced scale (aggregating phenotype by ecotype/taxa), this should be a matrix with 1 column.
}
 \item{SNPs}{
The SNP matrix on the reduced scale (aggregating phenotype by ecotype/taxa).
}
\item{significant}{
A vector of 0 and 1's where the 1's indicate what SNP's were found to be significant in the preselection function.
}
\item{number_cores}{
Number of cores to be passed on to the genetic algorithm to increase computational speed.
}
\item{principal_components}{
The principal component matrix on the reduced scale (aggregating phenotype by ecotype/taxa).
}
\item{maxiterations}{
This is the maximum number of iterations the Genetic Search algorithm will run.
}
\item{runs_til_stop}{
This is the number of consecutive iterations where the BIC is not improved before the genetic algorithm is stopped.
}
\item{kinship}{
Default is set at FALSE. If kinship model is desired, input a kinship matrix and this will search models with the kinship component.
}
\item{info}{
Default is set at FALSE. An information matrix where the first row is the chromosomes and the second row is the position information
}
}
\details{
This function will print out lines corresponding to the convergence of the genetic search algorithm.
}
\value{
A named matrix where the names correspond to the significant SNP's. This will usually out a matrix with a singular row, where the values of this row is 0 or 1 where 1 indicates significance in the final model and 0 indicates non significance in the final model. Sometimes this will output a matrix with multiple columns. This is because there is a SNP or multiple SNPs that when added to the model create rank deficiency issues in the model. Naturally rank deficient columns are forced out but the genetic algorithm is not smart enough to sort these.
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

postGWAS(Y = Y,SNPs = SNPs,number_cores = 1, significant = Significant_SNPs$Significant,principal_components = principal_comp,maxiterations = 100,runs_til_stop = 10,kinship = FALSE,info = FALSE)

}
