##' @title Time lag analysis
##'
##' @description Time lag analysis for a multivariate data set
##'
##' @details Time lag analysis involves computing the compositional
##' dissimilarity between samples at lag 1, at lag 2, etc, then the
##' least squares slope between \eqn{\sqrt{\mathrm{lag}}}{sqrt(lag)}
##' and dissimilarity.
##'
##' @param x an R object. Only objects of class \code{"dist"} are
##' currently supported.
##' @param ... additional arguments passed to other methods.
##'
##' @return Returns an object of class \code{"timelagAnalysis"}, a list
##' with components
##'
##' \describe{
##'   \item{data}{data frame; unpacked set of lags (\code{Lag}) and
##'     compositional dissimilarities (\code{Distance}).}
##' }
##'
##' @rdname timelag
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @examples
##' ## load analogue for Abernethy data set & distance()
##' if (require("analogue")) {
##'
##' ## Load Abernethy Forest data set
##' data("abernethy", package = "analogue")
##' ## Load Abernethy Forest data set
##'
##' ## Remove the Depth and Age variables
##' abernethy2 <- abernethy[, -(37:38)]
##'
##' ## time lag analysis
##' dij <- as.dist(distance(abernethy2, method = "chord"))
##' tla <- timelag(dij)
##' head(tla[[1]])
##' }
`timelag` <- function(x, ...)
    UseMethod("timelag", x)

##' @rdname timelag
##'
##' @export
`timelag.dist` <- function(x, ...) {
    x <- as.matrix(x)
    lmax <- nrow(x) - 1
    lags <- seq_len(lmax)
    tl <- lapply(lags, function(i, m) m[row(m) == (col(m) + i)], m = x)
    out <- list(data = data.frame(Lag = rep(lags, sapply(tl, length)),
                Distance = unlist(tl, recursive = FALSE,
                use.names = FALSE)))
    class(out) <- "timelagAnalysis"
    out
}

##' @title A time lag analysis plot.
##'
##' @description A plot of \eqn{\sqrt{\mathrm{lag}}}{sqrt(lag)} against
##' compositional dissimilarity.
##'
##' @details TODO - smoothers etc.
##'
##' @param x an object of class \code{"timelagAnalysis"}, the result of
##' a call to \code{\link{timelag}}.
##' @param ... additional arguments passed to \code{\link{plot}}.
##'
##' @return Invisibly returns its input \code{x}.
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @examples
##' ## load analogue for Abernethy data set & distance()
##' if (require("analogue")) {
##'
##' ## Load Abernethy Forest data set
##' data("abernethy", package = "analogue")
##' ## Load Abernethy Forest data set
##'
##' ## Remove the Depth and Age variables
##' abernethy2 <- abernethy[, -(37:38)]
##'
##' ## time lag analysis
##' dij <- as.dist(distance(abernethy2, method = "chord"))
##' tla <- timelag(dij)
##' plot(tla, pch = 19)
##' }
`plot.timelagAnalysis` <- function(x, ...) {
    dat <- unclass(x[[1]])
    plot(Distance ~ sqrt(Lag), dat, ...)
    invisible(x)
}
