\name{log_profile_likelihood_REML}
\alias{log_profile_likelihood_REML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{log_profile_likelihood_REML
%%  ~~function to do ... ~~
}
\description{
The log likelihood that needs to be optimized for the full kinship model.
}
\usage{
log_profile_likelihood_REML(x,t,y,d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The reduced design matrix with principal components, intercept and SNP of interest
}
\item{t}{
Tau value for the random effect term
}
\item{y}{
The reduced matrix for the response value of interest
}
\item{d}{
The spectral decomposition diagonal matrix of eigen values
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This returns the REML value
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
