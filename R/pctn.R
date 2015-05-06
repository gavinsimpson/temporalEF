##' @title Principal coordinates of temporal neighbours
##'
##' @description Computes the classic PCNM by the principal coordinate
##' analysis of a truncated distance matrix, but for a one-dimensional
##' process.
##'
##' @param x an R object. For \code{pctn} currently only a sorted
##' vector of time points.
##' @param threshold numeric; threshold beyond which the temporal
##' separation of smaples is considered equal. The default if no value
##' is supplied is to find the largest temporal separation between any
##' two points. Separations greater than the threshold are given a
##' notional separation of 4 times \code{threshold}.
##' @param distfun function or character string naming a function that
##' will be used to compute the temporal seperation between samples.
##' Defaults to \code{\link{dist}} for the Euclidean distance. See
##' Details for further information.
##' @param ... additional arguments passed to other methods or on to
##' \code{distfun}.
##'
##' @details
##' The default distance coefficient used to compute
##' temporal separation is the Euclidean distance. If you want to use
##' a different coefficient, you can supply a suitable function to
##' argument \code{distfun}. This should be a function that returns an
##' object of class \code{"dist"} or a square symmetric matrix that can
##' be coerced to one. Arguments can be passed to \code{distfun} via
##' \ldots.
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(50)
##' mod <- pctn(tp)
##' mod
##'
##' @importFrom vegan spantree
##'
##' @rdname pctn
##' @export

`pctn` <- function(x, ...) {
    UseMethod("pctn")
}

##' @rdname pctn
##'
##' @export
`pctn.default` <- function(x, threshold, distfun = dist, ...) {
    tol <- sqrt(.Machine$double.eps)
    distfun <- match.fun(distfun, ...)
    dij <- distfun(x)
    if(!(inherits(dij, "matrix") || inherits(dij, "dist"))) {
        stop("'distfun' must return a matrix or an object of class 'dist'")
    }
    if(missing(threshold)) {
        threshold <- max(spantree(dij)$dist)
    }
    dij[dij > threshold] <- threshold * 4
    op <- options(warn = -1)
    on.exit(options(op))
    pcnm <- cmdscale(dij, k = length(x) - 1, eig = TRUE)
    want <- pcnm$eig > tol
    retval <- list(vectors = pcnm$points[, seq_len(sum(want))],
                   lambda = pcnm$eig[want], threshold = threshold,
                   tp = x, FUN = distfun)
    ## need to add dimnames to objects
    names(retval$lambda) <- colnames(retval$vectors) <-
        paste("EF", seq_along(retval$lambda), sep = "")
    rownames(retval$vectors) <- if(is.null(tnams <- names(x))) {
        paste("T", seq_len(nrow(retval$vectors)), sep = "")
    } else {
        tnams
    }
    class(retval) <- "pctn"
    retval
}

##' @rdname pctn
##'
##' @export
##' @param digits numeric; number of digits to display in output.
`print.pctn` <- function(x, digits = 3, ...) {
    cat("\n")
    writeLines(strwrap("Principal Coordinates of Temporal Neighbours"))
    cat("\n")
    writeLines(paste("No. of Eigenfunctions:", length(x$lambda)))
    writeLines("Eigenvalues:")
    print(x$lambda, digits = digits)
}

##' @rdname pctn
##'
##' @param choices numeric; vector indicating which eigenfunctions
##' to return.
##'
##' @export
`scores.pctn` <- function(x, choices, ...) {
    if(missing(choices)) {
        choices <- seq_len(ncol(x$vectors))
    }
    x$vectors[, choices]
}

##' @rdname pctn
##'
##' @export
`eigenvals.pctn` <- function(x, ...) {
    out <- x$lambda
    class(out) <- "eigenvals"
    out
}

##' @title Plot PCTN eigenfunctions
##'
##' @param x an object of class \code{"pctn"}
##' @param pages numeric; the number of pages over which to spread
##' the plots of the individual eigenfunctions
##' @param ylim numeric vector of limits for the y-axis
##' @param xlab,ylab x and y-axis labels
##' @param ask logical; should plotting be paused between pages?
##' @param ... additional arguments passed to \code{\link{plot}}
##'
##' @return A plot on the currently active device
##'
##' @seealso \code{\link{pctn}} for creating PCTN objects
##'
##' @keywords hplot
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(50)
##' mod <- pctn(tp)
##' plot(mod, pages = 2)
##'
##' @export
##'
`plot.pctn` <- function(x, pages = 1, ylim, ylab = nams, xlab = "",
                        ask = FALSE, ...) {
    len <- length(x$lambda)
    np <- ceiling(len / pages)
    prc <- n2mfrow(np)
    op <- par(mar = c(3,4,1,1), ask = ask)
    on.exit(par(op))
    on.exit(layout(1), add = TRUE)
    layout(matrix(seq_len(prod(prc)), ncol = prc[2], nrow = prc[1]))
    nams <- colnames(x$vectors)
    if(missing(ylim))
        ylim <- range(x$vectors)
    for(i in seq_len(len)) {
        plot(x$tp, x$vectors[,i], type = "l", ylab = ylab[i], xlab = xlab,
             ylim = ylim, ...)
    }
    invisible(x)
}
