\name{plot.pctn}
\alias{plot.pctn}
\title{Plot PCTN eigenfunctions}
\usage{
\method{plot}{pctn}(x, pages = 1, ylim, ylab = nams, xlab = "",
  ask = FALSE, ...)
}
\arguments{
  \item{x}{an object of class \code{"pctn"}}

  \item{pages}{numeric; the number of pages over which to
  spread the plots of the individual eigenfunctions}

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
Plot PCTN eigenfunctions
}
\examples{
tp <- seq_len(50)
mod <- pctn(tp)
plot(mod, pages = 2)
}
\author{
Gavin L. Simpson
}
\seealso{
\code{\link{pctn}} for creating PCTN objects
}
\keyword{hplot}

