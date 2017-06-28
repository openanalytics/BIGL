#' Method for plotting response surface objects
#'
#' @param x Output of \code{\link{fitSurface}}
#' @param color Character indicating on what values surface coloring will be
#'   based.
#'
#'   If \code{color = "z-score"}, surface coloring will be based on median of
#'   standartized off-axis Z-scores. Median function can be replaced by other
#'   function using an optional \code{colorfun} argument which will be passed to
#'   \code{plotResponseSurface}. Color breaks are determined here by standard
#'   deviation of off-axis Z-scores. For \code{color = "maxR"}, coloring will be
#'   based on values of maxR statistic and the quantile of its distribution
#'   (bootstrapped or not). If \code{color = "occupancy"}, coloring will be
#'   based on calculated occupancy rate for the respective dose combination.
#'
#' @param ... Further parameters passed to \code{\link{plotResponseSurface}}.
#'   \code{colorBy} argument in this method is computed automatically and thus
#'   cannot be passed to \code{\link{plotResponseSurface}}.
#' @export
plot.ResponseSurface <- function(x, color = c("z-score", "maxR", "occupancy"), ...) {

  color <- match.arg(color)
  inputs <- as.list(substitute(list(...)))[-1L]

  ## Green is synergy, red is antagonism
  if(!exists("colorPalette", inputs)) {
    inputs$colorPalette <- c("red", rep("grey70", 2), "blue")
    if (x$fitResult$coef["b"] > x$fitResult$coef["m1"]) {
      inputs$colorPalette <- rev(inputs$colorPalette)
    }
  }

  if (color == "z-score") {
    boundary <- sd(x$offAxisTable[["z.score"]])
    inputs$colorBy <- x$offAxisTable[, c("d1", "d2", "z.score")]
    if (!exists("breaks", inputs)) inputs$breaks <- c(-Inf, -boundary, 0, boundary, Inf)
    if (!exists("main", inputs)) inputs$main <- "Z-scores"
  } else if (color == "maxR") {
    inputs$colorBy <- x$maxR$Ymean[, c("d1", "d2", "R")]
    q <- attr(x$maxR$Ymean, "q")
    if (!exists("breaks", inputs)) inputs$breaks <- c(-Inf, -q, 0, q, Inf)
    if (!exists("main", inputs)) inputs$main <- "maxR"
  } else if (color == "occupancy") {
    inputs$colorBy <- x$occupancy
    inputs$colorPalette <- c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5")
    if (!exists("breaks", inputs)) inputs$breaks <- c(0, 0.25, 0.5, 0.75, 1)
    if (!exists("main", inputs)) inputs$main <- "Occupancy rate"
  }

  inputs$data <- x$data
  inputs$fitResult <- x$fitResult
  inputs$transforms <- x$transforms
  inputs$null_model <- x$null_model

  do.call(plotResponseSurface, inputs)

}

#' Method for plotting of contours based on maxR statistics
#'
#' @param x Output of \code{\link{fitSurface}}
#' @param ... Further parameters passed to \code{\link{plot.maxR}}
#' @export
contour.ResponseSurface <- function(x, ...) {

  if (!exists("maxR", x))
    stop("maxR statistics were not found.")

  plot(x$maxR, ...)

}

#' Summary of \code{ResponseSurface} object
#'
#' @param object Output of \code{\link{fitSurface}}
#' @param ... Further parameters
#' @export
summary.ResponseSurface <- function(object, ...) {

  ans <- list()
  ans$marginalFit <- summary(object$fitResult)
  ans$null_model <- object$null_model
  ans$shared_asymptote <- object$fitResult$shared_asymptote

  if (!is.null(object$meanR)) ans$meanR <- summary(object$meanR)
  if (!is.null(object$maxR)) ans$maxR <- summary(object$maxR)

  ans$occup <- mean(object$occupancy$occupancy)

  class(ans) <- "summary.ResponseSurface"
  ans
}

#' Print method for the summary function of \code{ResponseSurface} object
#'
#' @param x Summary of \code{ResponseSurface} object
#' @param ... Further parameters
#' @export
print.summary.ResponseSurface <- function(x, ...) {

  cat("Null model: ")
  if (x$null_model == "loewe" & x$shared_asymptote == TRUE)
    cat("Standard Loewe Additivity")
  else if (x$null_model == "loewe" & x$shared_asymptote == FALSE)
    cat("Generalized Loewe Additivity")
  else if (x$null_model == "hsa" & x$shared_asymptote == TRUE)
    cat("Highest Single Agent with shared maximal response")
  else if (x$null_model == "hsa" & x$shared_asymptote == FALSE)
    cat("Highest Single Agent with differing maximal response")
  else
    cat(x$null_model)

  cat("\n")
  cat("Mean occupancy rate:", x$occup)
  cat("\n\n")
  print(x$marginalFit)
  cat("\n")

  if (!is.null(x$meanR)) print(x$meanR)
  if (!is.null(x$maxR)) print(x$maxR)

  if (is.null(x$meanR) & is.null(x$maxR)) {
    cat("\n\n")
    cat("No test statistics were computed.")
  }

  cat("\n")
}


#' Predicted values of the response surface according to the given null model
#'
#' @param object Output of \code{\link{fitSurface}}
#' @param ... Further parameters
#' @export
fitted.ResponseSurface <- function(object, ...) {

  doseInput <- object$data[, c("d1", "d2")]
  parmInput <- coef(object$fitResult)

  switch(object$null_model,
         "loewe" = generalizedLoewe(doseInput, parmInput)$response,
         "hsa" = hsa(doseInput, parmInput))

}
