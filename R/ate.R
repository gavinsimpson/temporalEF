##' Asymmetric temporal eigenfunctions
##'
##' Generate a set of asymmetric temporal eigenfunctions
##'
##' @details The asymmetric eigenvector map (AEM) is a recently proposed
##' method for describing orthongonal spatial functions for use in
##' multivariate ordination. AEMs are asymmetric because they apply a
##' directionality to the spatial dependencies modelled by the
##' eigenfunctions. Asymmetric temporal eigenfunctions (ATEs) implement
##' the AEM idea to model patterns of temporal dependence; i.e. a single
##' spatial dimension.
##'
##' @rdname ate
##'
##' @param x an R object. For \code{ate} currently only a sorted vector
##' of time points. For the \code{print} method an object of class
##' \code{"ate"}.
##' @param N numeric; the number of eigenvectors to return. If not
##' supplied, N is taken from the appropriate dimension of \code{x};
##' for the default method, this is the length of \code{x}.
##' @param weight logical; should a weighting matrix be generated and
##' applied to the link matrix?
##' @param FUN a function to be applied to the weighting matrix.
##' Ignored if \code{weight} is \code{FALSE}.
##' @param link0 logical; should the link from t[0] be included?
##' @param ... additional arguments passed to \code{FUN}.
##'
##' @importFrom vegan scores eigenvals
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @examples
##' tp <- seq_len(10)
##' tefs <- ate(tp)
##' tefs.I <- moranI(tefs)
##' plot(tefs.I)

`ate` <- function(x, ...) {
    UseMethod("ate")
}

##' @rdname ate
##'
##' @export
`ate.default` <- function(x, N, weight = FALSE, FUN = NULL,
                          link0 = TRUE, ...) {
    ## check we have at least 2 time points
    lenx <- length(x)
    if (lenx < 2) {
        stop("At least two time points are required for ATE.")
    }
    if(missing(N))
        N <- lenx
    if(!link0)
        N <- N-1
    ## create the links matrix, store twice, once for links & once
    ## for computations as latter may be weighted
    links <- lmat <- makeLinks(x, link0 = link0)
    if(weight) {
        w <- makeWeightVec(x, FUN = FUN, link0 = link0, ...)
        lmat <- lmat * w
    } else {
        w <- rep(1, N)
        attr(w, "link0") <- link0
    }
    ## centre lmat
    lmat <- sweep(lmat, 2, colMeans(lmat), "-")
    ## decomposition
    SVD <- svd(lmat, nu = N, nv = 0)
    ## keep non-zero SVs
    eps <- sqrt(.Machine$double.eps)
    k <- SVD$d > eps
    lambda <- SVD$d[k]^2 / N
    vectors <- SVD$u[, k, drop = FALSE]
    ## normalize the AETs
    vectors <- sweep(vectors, 2, sqrt(colSums(vectors^2)), "/")
    ## add some names
    names(lambda) <- colnames(vectors) <- paste("EF", seq_along(lambda),
                                                sep = "")
    rownames(vectors) <-
        if(is.null(tnams <- names(x))) {
            paste("T", seq_len(nrow(vectors)), sep = "")
        } else {
            tnams
        }
    retval <- list(vectors = vectors, tp = x, weights = w, links = links,
                   lambda = lambda, FUN = FUN)
    class(retval) <- "ate"
    retval
}

##' @rdname ate
##'
##' @export
##' @param digits numeric; number of digits to display in output.
`print.ate` <- function(x, digits = 3, ...) {
    ev <- eigenvals(x)
    cat("\n")
    writeLines(strwrap("Asymmetric Temporal Eigenfunctions"))
    cat("\n")
    writeLines(paste("No. of Eigenfunctions:", length(ev)))
    writeLines("Eigenvalues:")
    print(ev, digits = digits)
}

##' @rdname ate
##'
##' @param choices numeric; vector indicating which eigenfunctions
##' to return.
##'
##' @export
`scores.ate` <- function(x, choices, ...) {
    efs <- eigenfuns(x)
    if(missing(choices)) {
        choices <- seq_len(ncol(efs))
    }
    efs[, choices, drop = FALSE]
}

##' @rdname ate
##'
##' @export
`eigenvals.ate` <- function(x, ...) {
    out <- x$lambda
    class(out) <- "eigenvals"
    out
}

##' Plot asymmetric temporal eigenfunctions
##'
##' A multi-panel layout showing the calculated tempoeral eigenfunctions.
##'
##' @param x an object of class \code{"ate"}
##' @param pages numeric; the number of pages over which to spread
##' the plots of the individual eigenfunction
##' @param ylim numeric vector of limits for the y-axis
##' @param xlab,ylab x and y-axis labels
##' @param ask logical; should plotting be paused between pages?
##' @param ... additional arguments passed to \code{\link{plot}}
##'
##' @return A plot on the currently active device
##'
##' @seealso \code{\link{ate}} for creating ATE objects
##'
##' @keywords hplot
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(50)
##' mod <- ate(tp)
##' plot(mod, pages = 2)
##'
##' @export
##'
`plot.ate` <- function(x, pages = 1, ylim, ylab = nams, xlab = "",
                       ask = FALSE, ...) {
    ev <- eigenvals(x)
    len <- length(ev)
    np <- ceiling(len / pages)
    prc <- n2mfrow(np)
    op <- par(mar = c(3,4,1,1), ask = ask)
    on.exit(par(op))
    on.exit(layout(1), add = TRUE)
    layout(matrix(seq_len(prod(prc)), ncol = prc[2], nrow = prc[1]))
    efs <- eigenfuns(x)
    nams <- colnames(efs)
    if(missing(ylim))
        ylim <- range(efs)
    for(i in seq_len(len)) {
        plot(x$tp, efs[,i], type = "l", ylab = ylab[i], xlab = xlab,
             ylim = ylim, ...)
    }
    invisible(x)
}
