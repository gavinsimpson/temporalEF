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
##' @param permute logical; test \eqn{I} via a permutation test
##' @param alternative character; the type of test to perform. The default is \code{"two.sided"}, which allows for both positive and negative spatial correlation (values of \eqn{I}). If you anticipate only one of positive or negative spatial association then you can specify \code{"greater"} or \code{"less"}, respectively, as the alternative hypothesis.
##' @param nperm numeric; number of permutations to perform in permutation test if requested.
##'
##' @export
`moranI.ate` <- function(x, permute = FALSE,
                         alternative = c("two.sided", "greater", "less"),
                         nperm = 999, ...) {
    v <- x$vectors
    n <- NROW(v)
    w <- x$weights
    if (attr(w, "link0")) {
        w <- w[-1]
    }
    alternative <- match.arg(alternative)
    stat <- apply(v, 2, calcMoranI, w = w, permute = permute, ...)
    stat <- do.call("rbind", stat)

    ## expected value of Moran's I
    expI <- -1 / (n - 1)

    out <- list(statistic = stat, expected = expI)
    class(out) <- "moranI"
    out
}

##' @export
##' @importFrom stats symnum
`print.moranI` <- function(x,
                           digits = max(3, getOption("digits") - 2),
                           eps.Pvalue = .Machine$double.eps,
                           na.print = "NA",
                           signif.stars = getOption("show.signif.stars"),
                           signif.legend = signif.stars, ...) {
    cat("\n")
    writeLines(strwrap("Moran's I of Temporal Eigenfunctions"))
    cat("\n")
    stat <- unlist(x$statistic[,1])
    pval <- unlist(x$statistic[,2])
    statistic <- cbind("Moran's I" = stat, "p value" = pval)
    statistic[,1] <- format(round(unlist(statistic[,1]), digits = digits),
                            digits = digits)
    statistic[,2] <- format.pval(pval, digits = digits, eps = eps.Pvalue)
    Signif <- symnum(pval, corr = FALSE, na = FALSE,
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " "))
    statistic <- cbind(statistic, format(Signif))
    print.default(statistic, quote = FALSE, right = TRUE, na.print = na.print,
                  ...)
    if (signif.stars && signif.legend) {
        if ((w <- getOption("width")) < nchar(sleg <- attr(Signif, "legend")))
            sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
        cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
            4 + max(nchar(sleg, "bytes") - nchar(sleg)))
    }
    invisible(x)
}

##' @importFrom stats pnorm sd
`calcMoranI` <- function(y, w, scale = FALSE, norm = FALSE,
                         permute = FALSE,
                         alternative = "two.sided",
                         nperm = 999) {
    ## Implementation based on equations in Gittleman & Kot (1990) Systematic
    ## Zoology 39(3):227--241
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
    pval <- NULL
    ## parametric, assuming I is distributed Gaussian
    if (!isTRUE(permute)) {
        if (alternative %in% c("less","greater")) {
            alternative <- c("less","greater")[(obsI > expI) + 1]
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
            pval <- 2 * pnorm(-abs(obsI), mean = expI, sd = sdev)
        } else if (isTRUE(all.equal(alternative, "greater"))) {
            pval <- pnorm(obsI, mean = expI, sd = sdev, lower.tail = FALSE)
        } else {
            pval <- pnorm(obsI, mean = expI, sd = sdev)
        }
    } else {
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
        } else {
            sum(abs(perm) >= abs(obsI))
        }
    }

    list(statistic = obsI, p.value = pval)
}
##' @title Plot Moran's I statistics for temporal eigenfunctions
##'
##' @description A plot of Moran's I versus time.
##' @param x an object of class \code{\link{moranI}}.
##' @param alpha numeric; level of significance
##' @param type character; the type of plotting to use. See
##' \code{\link{plot.default}}.
##' @param xlab the label for the x-axis of the plot.
##' @param ylab the label for the y-axis of the plot.
##' @param pch the plotting character to use
##' @param bg fill colour of points
##' @param col border colour of points
##' @param ... additional arguments passed to \code{\link{plot}}.
##'
##' @return A plot on the current device
##'
##' @author Gavin L. Simpson
##'
##' @export
`plot.moranI` <- function(x, alpha = 0.05,
                          type = "b",
                          xlab = "Eigenfunctions", ylab = "Moran's I",
                          pch = 21, bg = "red", col = "black", ...) {
    statistic <- x$statistic[,1]
    signif <- x$statistic[,2] > alpha
    n <- NROW(statistic)
    col <- rep(col, length.out = n)
    bg <- rep(bg, length.out = n)
    bg[signif] <- "transparent"
    xvals <- seq_along(statistic)
    plot(xvals, statistic, type = "n", ..., axes = FALSE,
         ylab = ylab, xlab = xlab,
         panel.first = abline(h = x$expected, col = "red"))
    points(xvals, statistic, type = type, pch = pch, col = col, bg = bg, ...)
    axis(side = 2)
    axis(side = 1, labels = paste0("EF", xvals), at = xvals)
    box()
    invisible(x)
}
