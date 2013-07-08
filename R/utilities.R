##' Build an asymmetric link matrix
##'
##' @param tp numeric vector of sorted time points
##' @param link0 logical; should the link from t[0] be included?
##' @param dimnames logical; should dimnames be attached to the link
##' matrix?
##' @return binary matrix of links between time points. The object is
##' a square asymmetric matrix with \code{length(tp)} rows and columns.
##' Rows represent the \code{tp} time points, and columns are the
##' network links between timepoints.
##' @seealso \code{\link{makeWeights}} for a matching weight matrix.
##' @export
##' @keywords utilities
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
    lmat
}

##' Build a temporal weighting matrix
##'
##' @param tp numeric vector of sorted time points
##' @param FUN a function to apply
##' @param link0 logical; should the link from t[0] be included?
##' @param dimnames logical; should dimnames by attachedd to the
##' matrix?
##' @param ... optional arguments passed to \code{FUN}
##' @return a square, symmetric matrix of weights.
##' @seealso \code{\link{makeLinks}} for a matching link matrix
##' @export
##' @keywords utilities
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
    w
}