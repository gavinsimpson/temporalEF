\name{makeWeights}
\alias{makeWeights}
\title{Build a temporal weighting matrix}
\usage{
  makeWeights(tp, FUN = NULL, link0 = TRUE,
    dimnames = TRUE, ...)
}
\arguments{
  \item{tp}{numeric vector of sorted time points}

  \item{FUN}{a function to apply}

  \item{link0}{logical; should the link from t[0] be
  included?}

  \item{dimnames}{logical; should dimnames by attachedd to
  the matrix?}

  \item{...}{optional arguments passed to \code{FUN}}
}
\value{
  a square, symmetric matrix of weights.
}
\description{
  Build a temporal weighting matrix
}
\examples{
tp <- seq_len(10)
makeWeights(tp)
}
\seealso{
  \code{\link{makeLinks}} for a matching link matrix
}
\keyword{utilities}

