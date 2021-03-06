% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moranI.R
\name{moranI}
\alias{moranI}
\alias{moranI.ate}
\title{Moran's I for temporal eigenfunctions}
\usage{
moranI(x, ...)

\method{moranI}{ate}(x, permute = FALSE, alternative = c("two.sided",
  "greater", "less"), nperm = 999, ...)
}
\arguments{
\item{x}{an R object of class \code{"ate"}.}

\item{...}{additional arguments passed to methods.}

\item{permute}{logical; test \eqn{I} via a permutation test}

\item{alternative}{character; the type of test to perform. The default is \code{"two.sided"}, which allows for both positive and negative spatial correlation (values of \eqn{I}). If you anticipate only one of positive or negative spatial association then you can specify \code{"greater"} or \code{"less"}, respectively, as the alternative hypothesis.}

\item{nperm}{numeric; number of permutations to perform in permutation test if requested.}
}
\value{
Numeric vector of Moran's I statistics
}
\description{
Calculate Moran's I for tempoeral eigenfunctions and assess significant of Moran's I values via a parametric test that assumes asymptotic normality or a permutation test.
}
\author{
Gavin L. Simpson
}
