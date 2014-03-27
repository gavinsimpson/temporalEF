\name{plot.moranI}
\alias{plot.moranI}
\title{Plot Moran's I statistics for temporal eigenfunctions}
\usage{
\method{plot}{moranI}(x, type = "b", xlab = "Eigenfunctions",
  ylab = "Moran's I", ...)
}
\arguments{
  \item{x}{an object of class \code{\link{moranI}}.}

  \item{type}{character; the type of plotting to use. See
  \code{\link{plot.default}}.}

  \item{xlab}{the label for the x-axis of the plot.}

  \item{ylab}{the label for the y-axis of the plot.}

  \item{...}{additional arguments passed to
  \code{\link{plot}}.}
}
\value{
A plot on the current device
}
\description{
A plot of Moran's I versus time.
}
\author{
Gavin L. Simpson
}

