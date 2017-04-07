##' Extract temporal eigenfunctions
##'
##' Simple extractor function to access temporal eigenfunctions from various objects. Currently methods for class(es) `"ate"`, `"pctn"`, and `"tef"` are provided.
##'
##' @param x the object from which to extract temporal eigenfunctions
##' @param ... additional arguments passed to methods
##'
##' @return A data frame with one column per eigenfunction
##'
##' @author Gavin L. Simpson
##'
##' @rdname eigenfuns
##' @export
`eigenfuns` <- function(x, ...) {
    UseMethod("eigenfuns")
}

##' @rdname eigenfuns
##' @param take numeric; which eigenfunctions should be returned
##' @export
`eigenfuns.ate` <- function(x, take = NULL, ...) {
    efs <- as.data.frame(x$vectors)

    if (!is.null(take)) {
        nc <- NCOL(efs)
        if (any(miss <- !take %in% seq_len(nc))) {
            stop(paste("Requested TEFs (",
                       paste(take[miss], collapse = ", "),
                       ") are not available.", sep = ""))
        }
        efs <- efs[, take, drop = FALSE]
    }

    efs
}

##' @export
##' @rdname eigenfuns
`eigenfuns.pctn` <- function(x, take = NULL, ...) {
    efs <- as.data.frame(x$vectors)

    if (!is.null(take)) {
        nc <- NCOL(efs)
        if (any(miss <- !take %in% seq_len(nc))) {
            stop(paste("Requested TEFs (",
                       paste(take[miss], collapse = ", "),
                       ") are not available.", sep = ""))
        }
        efs <- efs[, take, drop = FALSE]
    }

    efs
}

##' @export
##' @rdname eigenfuns
`eigenfuns.tef` <- function(x, take = NULL, ...) {
    efs <- as.data.frame(x$tefs)

    if (!is.null(take)) {
        nc <- NCOL(efs)
        if (any(miss <- !take %in% seq_len(nc))) {
            stop(paste("Requested TEFs (",
                       paste(take[miss], collapse = ", "),
                       ") are not available.", sep = ""))
        }
        efs <- efs[, take, drop = FALSE]
    }

    efs
}
