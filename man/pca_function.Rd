\name{pca_function}
\alias{pca_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{pca_function
%%  ~~function to do ... ~~
}
\description{
Create Principal Components from Standardized Set of SNPs.
}
\usage{
pca_function(SNPs,number_components,plot_it)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{
The SNP matrix where the columns consist of 0 and 1's.
}
\item{number_components}{
The number of principal components desired, if you don't know put a random value and use plot_it = TRUE.
}
\item{plot_it}{
A TRUE FALSE logical equality, if TRUE will plot the percent variation explained by the components, if FALSE will not create a plot. In both scenarios this function will return a matrix of principal components.
}
}
\details{
This will work with both the full SNP matrix and the reduced SNP matrix. If you use the full SNP matrix you will have to reduce it yourself and this is at a higher computational burden. If you use the reduced SNP matrix you will get the same values as if you aggregated the principal components from the full SNP matrix, but this will be much faster.
}
\value{
\item{Plot}{A plot of the percent variation explained by the components when plot_it = TRUE}
\item{Matrix}{A matrix that the number of columns is the number of principal components and the number of rows is the same number of rows as the inputted data matrix.}
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
fullPreprocess <- preprocess_SNPs(SNPs = SNPs,Y = Y,MAF = 0.01,number_cores = 1)
SNPs <- fullPreprocess$SNPs
Y <- fullPreprocess$Y
fullPreprocess$SNPs_Dropped

pca_function(SNPs = SNPs,number_components = 3,plot_it = TRUE)
}
