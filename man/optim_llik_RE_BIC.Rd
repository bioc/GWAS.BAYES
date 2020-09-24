\name{optim_llik_RE_BIC}
\alias{optim_llik_RE_BIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{optim_llik_RE_BIC
%%  ~~function to do ... ~~
}
\description{
Calculates the BIC in the full kinship model scenario
}
\usage{
optim_llik_RE_BIC(x,y,d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The reduced design matrix with principal components, intercept and SNP of interest
}
  \item{y}{
The reduced response matrix.
}
\item{d}{
The spectral decomposition diagonal matrix of eigen values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns the BIC for the model with this SNP
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
