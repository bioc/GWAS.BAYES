\name{ga_modelselection_pcs_new}
\alias{ga_modelselection_pcs_new}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ga_modelselection_pcs_new
%%  ~~function to do ... ~~
}
\description{
Performs GA model selection to identify the best model when principal components are involved
}
\usage{
ga_modelselection_pcs_new(Y,X,regions,regionsnames,significant,number_cores,principal_components,maxiterations,runs_til_stop,kinship = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The phenotype response on the reduced scale (aggregating phenotype by ecotype/taxa), this should be a matrix with 1 column.
}
 \item{X}{
The SNP matrix on the reduced scale (aggregating phenotype by ecotype/taxa).
}
\item{significant}{
A vector of 0 and 1's where the 1's indicate what SNP's were found to be significant in the preselection function.
}
\item{regions}{
A matrix where each column represents a principal component for each region.
}
\item{regionsnames}{
A named list which highlights which SNPs fall into which region.
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
This is the numebr of consectutive iterations where the BIC is not improved before the genetic algorithm is stopped.
}
\item{kinship}{
Default is set at FALSE. If kinship model is desired, input a kinship matrix and this will search models with the kinship component.
}
}
\details{
This function will print out lines correpsonding to the convergence of the genetic search algorithm.
}
\value{
A named matrix where the names corespond to the significant SNP's. This will usually out a matrix with a singular row, where the values of this row is 0 or 1 where 1 indicates significance in the final model and 0 indicates non significance in the final model. Sometimes this will output a matrix with mulitple columns. This is because there is a SNP or multiple SNPs that when added to the model create rank deficiency issues in the model. Naturally rank deficient columns are forced out but the genetic algoritm is not smart enough to sort these.
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
