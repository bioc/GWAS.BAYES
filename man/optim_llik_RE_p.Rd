\name{optim_llik_RE_p}
\alias{optim_llik_RE_p}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{optim_llik_RE_p
%%  ~~function to do ... ~~
}
\description{
This will calculate the p-value and perform the optimization in the full kinship model case
}
\usage{
optim_llik_RE_p(x,y,d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Reduced design matrix with principal components, intercept and SNP of interest
}
\item{y}{
The reduced matrix of the response value of interest
}
\item{d}{
The spectral decomposition diagonal matrix of eigen values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a p-value for the specified data
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
