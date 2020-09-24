\name{Pval_function}
\alias{Pval_function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pval_function
%%  ~~function to do ... ~~
}
\description{
Performs multiple comparison corrections on p-values and returns significant SNP's. This is used internally for the preselection function.
}
\usage{
Pval_function(p_vals,n,thresh,control)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p_vals}{
A vector of p-values calculated by the preselection function.
}
\item{n}{
The original number of SNPs
}
\item{thresh}{
The type 1 error rate
}
\item{control}{
The multiple comparison correction one would like to apply.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The results are formated as a data.frame with the column Significant being 1 or 0 depending on if the SNP was significant (1 for significant). The P_values column will be the p-values that were calculated for each SNP.
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
