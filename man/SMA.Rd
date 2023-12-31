% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BICOSS.R
\name{SMA}
\alias{SMA}
\title{Performs Single Marker Association tests for both Linear Mixed Models and Linear models.}
\usage{
SMA(Y, SNPs, kinship = FALSE, P3D = FALSE)
}
\arguments{
\item{Y}{The observed numeric phenotypes}

\item{SNPs}{The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.}

\item{kinship}{The observed kinship matrix, has to be a square positive semidefinite matrix. Defaulted as the identity matrix. The function used
to create the kinship matrix used in the BICOSS paper is A.mat() from package rrBLUP.}

\item{P3D}{Population previous determined, if TRUE BICOSS uses approximated variance parameters estimated from the baseline model when conducting
both the screening and the model selection steps. Setting P3D = TRUE is significantly faster. If FALSE, uses exact estimates of the variance
parameters all models in both the screening and model selection step.}
}
\value{
The p-values corresponding to every column provided in SNPs. These p-values can be used with any threshold of your choosing or with
p.adjust().
}
\description{
Performs Single Marker Association tests for both Linear Mixed Models and Linear models.
}
