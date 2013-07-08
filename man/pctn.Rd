\name{pctn}
\alias{eigenvals.pctn}
\alias{pctn}
\alias{pctn.default}
\alias{print.pctn}
\alias{scores.pctn}
\title{Principal coordinates of temporal neighbours}
\usage{
  pctn(x, ...)

  \method{pctn}{default} (x, threshold, distfun = dist,
    ...)

  \method{print}{pctn} (x, digits = 3, ...)

  \method{scores}{pctn} (x, choices, ...)

  \method{eigenvals}{pctn} (x, ...)
}
\arguments{
  \item{x}{an R object. For \code{pctn} currently only a
  sorted vector of time points.}

  \item{threshold}{numeric; threshold beyond which the
  temporal separation of smaples is considered equal. The
  default if no value is supplied is to find the largest
  temporal separation between any two points. Separations
  greater than the threshold are given a notional
  separation of 4 times \code{threshold}.}

  \item{distfun}{function or character string naming a
  function that will be used to compute the temporal
  seperation between samples. Defaults to
  \code{\link{dist}} for the Euclidean distance. See
  Details for further information.}

  \item{...}{additional arguments passed to other methods
  or on to \code{distfun}.}

  \item{digits}{numeric; number of digits to display in
  output.}

  \item{choices}{numeric; vector indicating which
  eigenfunctions to return.}
}
\description{
  Computes the classic PCNM by the principal coordinate
  analysis of a truncated distance matrix, but for a
  one-dimensional process.
}
\details{
  The default distance coefficient used to compute temporal
  separation is the Euclidean distance. If you want to use
  a different coefficient, you can supply a suitable
  function to argument \code{distfun}. This should be a
  function that returns an object of class \code{"dist"} or
  a square symmetric matrix that can be coerced to one.
  Arguments can be passed to \code{distfun} via \ldots.
}
\examples{
tp <- seq_len(50)
mod <- pctn(tp)
mod
}

