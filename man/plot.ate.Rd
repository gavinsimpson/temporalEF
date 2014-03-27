\name{plot.ate}
\alias{plot.ate}
\title{Plot asymmetric temporal eigenfunctions}
\usage{
\method{plot}{ate}(x, pages = 1, ylim, ylab = nams, xlab = "",
  ask = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{"ate"}}

  \item{pages}{numeric; the number of pages over which to
  spread the plots of the individual eigenfunction}

  \item{ylim}{numeric vector of limits for the y-axis}

  \item{xlab,ylab}{x and y-axis labels}

  \item{ask}{logical; should plotting be paused between
  pages?}

  \item{...}{additional arguments passed to
  \code{\link{plot}}}
}
\value{
A plot on the currently active device
}
\description{
Plot asymmetric temporal eigenfunctions
}
\examples{
tp <- seq_len(50)
mod <- ate(tp)
plot(mod, pages = 2)
}
\author{
Gavin L. Simpson
}
\seealso{
\code{\link{ate}} for creating ATE objects
}
\keyword{hplot}

