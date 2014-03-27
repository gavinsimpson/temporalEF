\name{timelag}
\alias{timelag}
\alias{timelag.dist}
\title{Time lag analysis}
\usage{
timelag(x, ...)

\method{timelag}{dist}(x, ...)
}
\arguments{
  \item{x}{an R object. Only objects of class \code{"dist"}
  are currently supported.}

  \item{...}{additional arguments passed to other methods.}
}
\value{
Returns an object of class \code{"timelagAnalysis"}, a list
with components

\describe{ \item{data}{data frame; unpacked set of lags
(\code{Lag}) and compositional dissimilarities
(\code{Distance}).} }
}
\description{
Time lag analysis for a multivariate data set
}
\details{
Time lag analysis involves computing the compositional
dissimilarity between samples at lag 1, at lag 2, etc, then
the least squares slope between
\eqn{\sqrt{\mathrm{lag}}}{sqrt(lag)} and dissimilarity.
}
\examples{
## load analogue for Abernethy data set & distance()
require("analogue")

## Load Abernethy Forest data set
data("abernethy", package = "analogue")
## Load Abernethy Forest data set

## Remove the Depth and Age variables
abernethy2 <- abernethy[, -(37:38)]

## time lag analysis
tla <- timelag(as.dist(distance(abernethy2, method = "chord")))
head(tla[[1]])
}
\author{
Gavin L. Simpson
}

