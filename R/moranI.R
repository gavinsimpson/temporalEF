##' @title Moran's I for temporal eigenfunctions
##'
##' @param x an R object of class \code{"ate"}.
##' @param ... additional arguments passed to methods.
##'
##' @return Numeric vector of Moran's I statistics
##'
##' @rdname moranI
##'
##' @author Gavin L. Simpson
##'
##' @export
##'
##' @importFrom permute shuffleSet
##'
`moranI` <- function(x, ...) {
    UseMethod("moranI")
}

##' @rdname moranI
##'
##' @method moranI ate
##' @S3method moranI ate
`moranI.ate` <- function(x, ...) {
    v <- x$vectors
    n <- NROW(v)
    w <- x$weights
    if (attr(w, "link0")) {
        w <- w[-1]
    }
    stat <- apply(v, 2, calcMoranI, w = w, ...)
    stat <- do.call(rbind, stat)

    ## expected value of Moran's I
    expI <- -1 / (n - 1)

    out <- list(statistic = stat, expected = expI)
    class(out) <- "moranI"
    out
}

`calcMoranI` <- function(y, w, scale = FALSE, norm = FALSE,
                         test = c("none", "parametric", "permutation"),
                         alternative = c("two.sided", "greater", "less"),
                         nperm = 999) {
    funI <- function(ybar, W, n) {
        num <- sum(W * (ybar %o% ybar))
        den <- sum(ybar^2)
        obsI <- (n / s0) * (num / den) ## observed I, G & K Eqn (1)
        obsI
    }

    n <- length(y) ## number of sites

    ## make a matrix from weights; first upper off-diagonal
    W <- matrix(0, ncol = n, nrow = n)
    W[row(W) == col(W) - 1] <- w

    rsum <- rowSums(W)
    if (norm) {
        ind <- rsum > 0
        W[ind, ] <- W[ind, ] / rsum[ind]
    }

    s0 <- sum(w) ## G & K Eqn (2)

    ybar <- y - mean(y)  ## centre y

    obsI <- funI(ybar, W, n)

    if (scale) {
        limI <- (n / s0) * (sd(rsum * ybar) / sd(ybar)) ## G & K Eqn (4)
        obsI <- obsI / limI
    }

    ## expected value of Moran's I
    expI <- -1 / (n - 1)

    ## test Moran's I
    test <- match.arg(test)
    alternative <- match.arg(alternative)

    pval <- NULL
    ## parametric, assuming I is distributed Gaussian
    if (isTRUE(all.equal(test, "parametric"))) {
        alternative <- if (alternative %in% c("less","greater")) {
            c("less","greater")[(obsI > expI) + 1]
        }
        s1 <- 0.5 * sum((W + t(W))^2)                ## G & K Eqn (7)
        s2 <- sum((rowSums(W) + colSums(W))^2)       ## G & K Eqn (8)
        k <- (sum(ybar^4) / n) / (sum(ybar^2) / n)^2 ## kurtosis G & K Eqn (9)
        sdev <- (n * ((n^2 - 3 * n + 3) * s1 - n * s2 + 3 * s0^2) -
                 k * (n * (n - 1) * s1 - 2 * n * s2 + 6 * s0^2)) /
                     ((n - 1) * (n - 2) * (n - 3) * s0^2)
        sdev <- sqrt(sdev - (1 / ((n - 1)^2)))       ## G & K Eqn (6)
        ## pval
        if (isTRUE(all.equal(alternative, "two.sided"))) {
            pval <- if (obsI <= expI) {
                2 * pnorm(obsI, mean = expI, sd = sdev, lower.tail = FALSE)
            } else {
                2 * pnorm(obsI, mean = expI, sd = sdev)
            }
        } else if (isTRUE(all.equal(alternative, "greater"))) {
            pval <- pnorm(obsI, mean = expI, sd = sdev, lower.tail = FALSE)
        } else {
            pval <- pnorm(obsI, mean = expI, sd = sdev)
        }
    }

    ## permutation-based
    if (isTRUE(all.equal(test, "permutation"))) {
        alternative <- if (isTRUE(all.equal(alternative, "two.sided"))) {
            c("less","greater")[(obsI > expI) + 1]
        }
        perm <- t(shuffleSet(n, nset = nperm))
        perm[] <- y[perm]
        perm <- t(perm)
        permI <- apply(perm, 1, funI, W = W, n = n)
        if (scale) { ## scale but obsI was already scaled
            permI <- permI / limI
        }
        permI <- c(obsI, permI)
        pval <- if (isTRUE(all.equal(alternative, "less"))) {
            sum(permI <= obsI) / (nperm + 1)
        } else if (isTRUE(all.equal(alternative, "greater"))) {
            sum(permI >= obsI) / (nperm + 1)
        }
    }

    list(statistic = obsI, p.value = pval)
}
##' @title Plot Moran's I statistics for temporal eigenfunctions
##'
##' @description A plot of Moran's I versus time.
##' @param x an object of class \code{\link{moranI}}.
##' @param type character; the type of plotting to use. See
##' \code{\link{plot.default}}.
##' @param xlab the label for the x-axis of the plot.
##' @param ylab the label for the y-axis of the plot.
##' @param ... additional arguments passed to \code{\link{plot}}.
##' @return A plot on the current device
##' @author Gavin
`plot.moranI` <- function(x, type = "b", xlab = "Eigenfunctions",
                          ylab = "Moran's I", ...) {
    statistic <- x$statistic[,1]
    xvals <- seq_along(statistic)
    plot(xvals, statistic, type = "n", ..., axes = FALSE,
         ylab = ylab, xlab = xlab)
    abline(h = x$expected, col = "red")
    points(xvals, statistic, type = type)
    axis(side = 2)
    axis(side = 1, labels = paste0("EF", xvals), at = xvals)
    box()
    invisible(x)
}
