\name{tef}
\alias{tef}
\alias{tef.default}
\title{Temporal eigenfunction ordination}
\usage{
tef(x, ...)

\method{tef}{default}(x, index, method = c("ate", "aem", "pctn", "pcnm"),
  ordination = c("rda"), ...)
}
\arguments{
  \item{x}{community data matrix or data frame}

  \item{index}{numeric; the time ordering of the samples
  from which temporal eigenfunctions will be computed}

  \item{method}{character; which method to use to create
  the temporal eigenfunctions}

  \item{ordination}{character; the name of the constrained
  ordination function, from the \pkg{vegan} package to use.
  Only \code{\link{rda}} is support currently.}

  \item{...}{additional arguments passed to other methods.
  Arguments are also passed to the function generating the
  temporal eigenfunctions and the constrained ordination
  function.}
}
\value{
An object of class \code{"tef"}, a list with the following
components:

\describe{ \item{ordination}{the fitted constrained
ordination.} \item{tefs}{the computed set of temporal
eigenfunctions.} }
}
\description{
Fits an ordination using temporal eigenfunctions as
constraints.
}
\details{
TODO
}
\author{
Gavin L. Simpson
}

