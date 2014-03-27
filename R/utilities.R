##' @title Build an asymmetric link matrix
##'
##' @param tp numeric vector of sorted time points
##' @param link0 logical; should the link from t[0] be included?
##' @param dimnames logical; should dimnames be attached to the link
##' matrix?
##'
##' @return binary matrix of links between time points. The object is
##' a square asymmetric matrix with \code{length(tp)} rows and columns.
##' Rows represent the \code{tp} time points, and columns are the
##' network links between timepoints.
##' @seealso \code{\link{makeWeights}} for a matching weight matrix.
##' @export
##' @keywords utilities
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(10)
##' makeLinks(tp)
`makeLinks` <- function(tp, link0 = TRUE, dimnames = TRUE) {
    nt <- length(tp)
    st <- seq_len(nt)
    lmat <- matrix(0, ncol = nt, nrow = nt)
    lmat[lower.tri(lmat, diag = TRUE)] <- 1
    if(dimnames)
        dimnames(lmat) <- list(paste("t", st, sep = ""),
                               paste("l", st - 1, sep = ""))
    ## drop the link0?
    if(!link0) {
        lmat <- lmat[,-1]
    }
    class(lmat) <- c("linkMatrix", "matrix")
    attr(lmat, "link0") <- link0
    lmat
}

##' Build a temporal weighting matrix
##'
##' @title Temporal weight matrix
##' @param tp numeric vector of sorted time points
##' @param FUN a function to apply
##' @param link0 logical; should the link from t[0] be included?
##' @param dimnames logical; should dimnames by attached to the
##' matrix?
##' @param ... optional arguments passed to \code{FUN}
##' @return a square, symmetric matrix of weights.
##' @seealso \code{\link{makeLinks}} for a matching link matrix
##' @export
##' @keywords utilities
##'
##' @rdname makeWeights
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(10)
##' makeWeights(tp)
`makeWeights` <- function(tp, FUN = NULL, link0 = TRUE, dimnames = TRUE,
                          ...) {
    if(is.unsorted(tp, strictly = TRUE))
        stop("'tp' must be strictly sorted in increasing temporal order.")
    w <- as.matrix(dist(tp)) + 1
    w[upper.tri(w)] <- 0
    if(!is.null(FUN)) {
        FUN <- match.fun(FUN)
        w <- FUN(w, ...)
    }
    if(dimnames) {
        st <- seq_along(tp)
        dimnames(w) <- list(paste("t", st, sep = ""),
                            paste("l", st - 1, sep = ""))
    } else {
        dimnames(w) <- NULL
    }
    ## drop the link0?
    if(!link0) {
        w <- w[,-1]
    }
    class(w) <- c("weightMatrix", "matrix")
    attr(w, "link0") <- link0
    w
}

##' @title Weight vector for time points
##'
##' @description A vector of weights based on the inverse of time
##' duration between time points.
##'
##' @param tp numeric vector of sorted time points
##' @param FUN a function to apply
##' @param link0 logical; should the link from t[0] be included?
##' @param names logical; should names by attached to the vector?
##' @param ... optional arguments passed to \code{FUN}
##'
##' @return a matrix of weights of length \code{length(tp) - 1}
##'
##' @export
##'
##' @keywords utilities
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(10)
##' makeWeightVec(tp)
`makeWeightVec` <- function(tp, FUN = NULL, link0 = TRUE, names = TRUE,
                            ...) {
    if(is.unsorted(tp, strictly = TRUE))
        stop("'tp' must be strictly sorted in increasing temporal order.")
    w <- 1 / diff(tp)
    if (!is.null(FUN)) {
        FUN <- match.fun(FUN)
        w <- FUN(w, ...)
    }
    if (link0) {
        w <- c(1,w)
    }
    if(names) {
        st <- seq_along(w) - 1
        names(w) <- paste("l", st, sep = "")
    } else {
        names(w) <- NULL
    }
    class(w) <- c("weightVec", "numeric")
    attr(w, "link0") <- link0
    w
}

##' @title Diagonal weight matrix for time points
##'
##' @description A diagonal matrix of weights based on the inverse
##' of time duration between time points.
##'
##' @param tp numeric vector of sorted time points
##' @param FUN a function to apply
##' @param link0 logical; should the link from t[0] be included?
##' @param dimnames logical; should dimnames by attached to the
##' matrix?
##' @param ... optional arguments passed to \code{FUN}
##'
##' @return a diagonal matrix of weights with dimensions \code{length(tp) - 1}
##'
##' @export
##'
##' @keywords utilities
##'
##' @author Gavin L. Simpson
##'
##' @examples
##' tp <- seq_len(10)
##' makeWeightMat(tp)
`makeWeightMat` <- function(tp, FUN = NULL, link0 = TRUE, dimnames = TRUE,
                            ...) {
    w <- makeWeightVec(tp, FUN = FUN, link0 = link0, names = dimnames, ...)
    W <- diag(w)
    if (dimnames) {
        nams <- names(w)
        dimnames(W) <- list(nams, nams)
    }
    class(W) <- c("weightMat", "matrix")
    attr(w, "link0") <- link0
    W
}
