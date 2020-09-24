\name{cor_plot}
\alias{cor_plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
cor_plot(SNPs,significant,info = FALSE)
}
\description{
A function that creates correlation heatmaps for given significant SNPs from a SNP dataframe.
}
\usage{
cor_plot(SNPs,significant,info = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SNPs}{
A standardized SNP matrix where the columns take on the values of 0 or 1.
}
 \item{significant}{
A vector of 0's or 1's that contain which significant SNP's to look at.
}
\item{info}{
Default is FALSE. If information such as the chromosome and the position is known, then that can be entered here and will return a correlation heatmap with the labels chromosome - position. The format for this entry is the 2 row dataframe or matrix, where the first row is the chromosome and the second row is the position.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A correlation heatmap with the axis labels either SNP1 ... or chromosome - position
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
principal_comp <- pca_function(SNPs = SNPs,number_components = 1,plot_it = FALSE)
Significant_SNPs <- preselection(Y = Y, SNPs = SNPs,number_cores = 1, principal_components = principal_comp,frequentist = TRUE,controlrate = "bonferroni",threshold = .05,kinship = FALSE,info = FALSE)

cor_plot(SNPs = SNPs,significant = Significant_SNPs$Significant,info = FALSE)
}
