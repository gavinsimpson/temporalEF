##' @title Temporal eigenfunction ordination
##'
##' @description Fits an ordination using temporal eigenfunctions as
##' constraints.
##'
##' @details TODO
##'
##' @rdname tef
##'
##' @param x community data matrix or data frame
##' @param index numeric; the time ordering of the samples from which
##' temporal eigenfunctions will be computed
##' @param method character; which method to use to create the temporal
##' eigenfunctions
##' @param ordination character; the name of the constrained ordination
##' function, from the \pkg{vegan} package to use. Only \code{\link{rda}}
##' is support currently.
##' @param ... additional arguments passed to other methods. Arguments
##' are also passed to the function generating the temporal
##' eigenfunctions and the constrained ordination function.
##'
##' @return An object of class \code{"tef"}, a list with the following
##' components:
##'
##' \describe{
##'   \item{ordination}{the fitted constrained ordination.}
##'   \item{tefs}{the computed set of temporal eigenfunctions.}
##' }
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @importFrom vegan rda

`tef` <- function(x, ...) {
    UseMethod("tef", x)
}

##' @rdname tef
##'
##' @method tef default
##' @S3method tef default
`tef.default` <- function(x, index, method = c("ate", "aem", "pctn", "pcnm"),
                          ordination = c("rda"), ...) {
    method <- match.arg(method)
    ## allow the standard name AEM for ATE
    if (identical(method, "aem"))
        method <- "ate"
    ## allow the standard name AEM for ATE
    if (identical(method, "pcnm"))
        method <- "pctn"

    fun <- switch(method,
                  ate = ate,
                  pctn = pctn)
    tefs <- fun(index, ...)

    ## fit the constrained ordination
    ordination <- match.arg(ordination)
    oFun <- match.fun(ordination)
    ord <- oFun(X = x, Y = tefs$vectors, ...)

    ## return object
    out <- list(ordination = ord, tefs = tefs)
    class(out) <- "tef"
    out
}
