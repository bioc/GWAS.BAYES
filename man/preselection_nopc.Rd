\name{preselection_nopc}
\alias{preselection_nopc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{preselection_nopc
%%  ~~function to do ... ~~
}
\description{
Finds significant SNP's when no principal components are present.
}
\usage{
preselection_nopc(Y,X,number_cores,frequentist,controlrate,threshold,nullprob,alterprob,kinship = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
The reduced matrix of response values
}
\item{X}{
The reduced SNP matrix where th columns are either 1's or 0's.
}
\item{number_cores}{
The number of cores on which you would like to parallize this procedure
}
\item{frequentist}{
A logical value to see whether one would like to use a frequentist multiple comparison test or Bayesian False Discovery based on BIC's. The value of this affects whether values of the next parameters are needed.
}
\item{controlrate}{
Only used when frequentist = TRUE. This is for which multiple comparison method you would like to use. Examples are "bonferroni" and "BH". See p.adjust for a full list of methods.
}
\item{threshold}{
The value at which type 1 error rate is held at. .05 in most common literature. Used when frequentist is TRUE or FALSE
}
\item{nullprob}{
Used when frequentist = FALSE, the probability that is assigned to the null hypothesis.
}
\item{alterprob}{
Used when frequentist = FALSE, the probability that is assigned to the alternate hypothesis.
}
\item{kinship}{
The kinship matrix if a model with a kinship component is desired. If not set kinship = FALSE.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Frequentist Matrix}{The matrix of results when Frequentist = TRUE. The results are formated as a data.frame with the column Significant being 1 or 0 depending on if the SNP was significant (1 for significant). The P_values column will be the p-values that were calculated for each SNP.}
\item{Bayesian Matrix}{The matrix of results when Frequentist = FALSE. The results are formated as a data.frame with the column Significant being 1 or 0 depending on if the SNP was significant (1 for significant). The ApprPosteriorProbs column will be the Approximate Posterior Probabilities that were calculated for each SNP.}
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
