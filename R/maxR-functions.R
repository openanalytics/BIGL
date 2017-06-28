#' Summary of maxR object
#'
#' @param object Object of \code{"maxR"} class
#' @param ... Further arguments
#' @export
summary.maxR <- function(object, ...) {

  ans <- list()
  ans$call <- object$Call
  ans$points <- outsidePoints(object$Ymean)
  ans$totals <- data.frame("Call" = object$Call,
                           "Syn" = sum(object$Ymean$call == "Syn"),
                           "Ant" = sum(object$Ymean$call == "Ant"),
                           "Total" = nrow(object$Ymean))

  class(ans) <- append("summary.maxR", class(ans))
  ans
}

#' Print summary of maxR object
#'
#' @param x Summary of \code{"maxR"} object
#' @inheritParams summary.maxR
#' @export
print.summary.maxR <- function(x, ...) {

  if (x$call != "None") {

    x$points$`p-value` <- ifelse(x$points$`p-value` < 2e-16, "<2e-16",
                                 round(x$points$`p-value`, 5))

    cat("\nEvidence for effects in data: ", x$call, "\n", sep="")
    cat("Points with significant deviations from the null: \n")
    print(x$points)
    rownames(x$totals) <- ""
    cat("\nOverall maxR summary:\n")
    print(x$totals)
  } else {
    cat("MaxR test did not detect any deviations from the null.\n")
  }
}


#' Plot of maxR object
#'
#' @param x Output of \code{\link{maxR}}. This can also be \code{"maxR"}
#'   element in the output of \code{\link{fitSurface}}.
#' @param plevels Probability levels used to generate a color scale
#' @param cutoff Probability cutoff to use for range of colors
#' @param maxshow Forced value for range of colors
#' @param ... Further arguments that are passed to \code{\link{format}} function
#'   for formatting of axis labels
#' @inheritParams plotResponseSurface
#' @inheritParams graphics::title
#' @importFrom graphics axis filled.contour points title
#' @importFrom grDevices extendrange rgb
#' @export
plot.maxR <- function(x,
                      main = "Contour plot for maxR",
                      xlab = "Dose (Compound 1)", ylab = "Dose (Compound 2)",
                      colorPalette = c("blue", "white", "red"),
                      logScale = TRUE, zTransform = function(z) { z },
                      plevels = c(0.7, 0.8, 0.9, 0.95, 0.99, 0.999),
                      cutoff = max(plevels), maxshow = NULL, ...) {

  uniqueDoses <- with(x$Ymean,
                      list("d1" = sort(unique(d1)),
                           "d2" = sort(unique(d2))))
  doseGrid <- expand.grid(uniqueDoses)

  log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
  transformF <- if (logScale) log10T else function(z) z

  maxRvalues <- x$Ymean$R
  maxRvalues <- maxRvalues[match(with(doseGrid, paste(d1, d2)),
                                 with(x$Ymean, paste(d1, d2)))]
  maxRvalues[is.na(maxRvalues)] <- 0

  if (is.null(maxshow)) {
    ## We show all values above cutoff in the same color
    maxshow <- quantile(attr(x$Ymean, "distr"), cutoff)
    if (is.null(maxshow)) {
      maxshow <- 3.5
      warning("No `maxshow` parameter specified, so 3.5 is used.")
    }
  }

  origMaxRvalues <- maxRvalues # Keep unchanged for point sizes
  maxRvalues[maxRvalues > maxshow] <- maxshow
  maxRvalues[maxRvalues < -maxshow] <- -maxshow

  ## Levels to show on color-scale
  zlevels <- sort(quantile(attr(x$Ymean, "distr"),
                           c(plevels[plevels < cutoff], cutoff)) %o% c(-1,1))

  filled.contour(
    x = transformF(uniqueDoses$d1),
    y = transformF(uniqueDoses$d2),
    z = matrix(maxRvalues, sapply(uniqueDoses, length)),
    levels = zlevels,
    color.palette = colorRampPalette(colorPalette, space = "rgb"),
    plot.axes = {
      axis(1, at = transformF(uniqueDoses$d1),
           labels = format(uniqueDoses$d1, ...), cex.axis = 0.8)
      axis(2, at = transformF(uniqueDoses$d2),
           labels = format(uniqueDoses$d2, ...), cex.axis = 0.8)
      points(expand.grid(x = transformF(uniqueDoses$d1),
                         y = transformF(uniqueDoses$d2)),
             pch = 20, col = rgb(0, 0, 0, 0.3),
             ## Size proportional to maxR statistic
             cex = matrix(abs(origMaxRvalues)/2,
                          sapply(uniqueDoses, length)))
    },
    key.axes = {
      axis(4, at = c(0, zlevels), tck = -0.1,
           mgp = c(3, 0.3, 0), cex.axis = 0.7,
           labels = c(
             paste0("\u2191",
                    if (colorPalette[1] == "blue") "Antagonism" else "Synergy",
                    "\n\n \n\n",
                    "\u2193",
                    if (colorPalette[1] == "blue") "Synergy" else "Antagonism"),
             paste0("\u2265", cutoff),
             rev(plevels[plevels < cutoff]),
             plevels[plevels < cutoff],
             paste0("\u2265", cutoff))
           )
    },
    key.title = title(main = "p-values", line = 1, cex.main = 1),
    xlim = extendrange(transformF(uniqueDoses$d1)),
    ylim = extendrange(transformF(uniqueDoses$d2)),
    zlim = maxshow*c(-1, 1),
    main = main,
    xlab = xlab, ylab = ylab,
    bty = "n"
  )

}

#' List non-additive points
#'
#' List all points with corresponding p-values declared non-additive by the
#' maxR statistical test.
#'
#' @param maxR maxR statistics table returned by \code{Ymean} component from the
#'   output of \code{\link{maxR}} function. This can also be \code{"maxR"}
#'   element in the output of \code{\link{fitSurface}} function.
#' @param B Iterations to use for the distribution of the maxR statistic. This
#'   is only used if \code{Ymean} dataframe does not have a \code{"distr"} attribute
#'   attached as is normally done when using \code{\link{fitSurface}} or \code{\link{maxR}}
#'   function.
#' @return Returns a dataframe listing only dose combinations that exhibit
#'   significant deviations from the expected response surface.
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 2)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   CP <- CPBootstrap(data, fitResult, null_model = "loewe", B.CP = 5)
#'   statM <- maxR(data, fitResult, null_model = "loewe", CP = CP)
#'   outsidePoints(statM$Ymean)
outsidePoints <- function(maxR, B = 10000) {
  if (is.null(attr(maxR, "distr"))) {
    n1 <- nrow(maxR)
    df0 <- attr(maxR, "df0")
    stopifnot(length(n1) == 1)
    stopifnot(length(df0) == 1)

    sim1 <- abs(matrix(rt(B*n1, df = df0), ncol = n1, byrow = TRUE))
    distr <- apply(sim1, 1, max)  #q1 <- quantile(apply(sim1, 1, max), cutoff)
    f <- ecdf(distr)
  } else f <- attr(maxR, "distr")

  resDF <- maxR[maxR$sign, ]
  resDF$`p-value` <- 1-f(resDF$absR)
  if (nrow(resDF) > 0) resDF[, c(1:2, 4, 8, 7)] else NULL
}
