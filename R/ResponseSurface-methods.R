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
#'   If \code{color = "effect-size"}, coloring will be
#'   based on effect size for the respective dose combination.
#'
#' @param ... Further parameters passed to \code{\link{plotResponseSurface}}.
#'   \code{colorBy} argument in this method is computed automatically and thus
#'   cannot be passed to \code{\link{plotResponseSurface}}.
#' @export
plot.ResponseSurface <- function(x, color = c("z-score", "maxR", "occupancy", "effect-size"), ...) {

  color <- match.arg(color)
  inputs <- as.list(substitute(list(...)))[-1L]

  ## Blue is synergy, red is antagonism
  if(!exists("colorPalette", inputs)) {
    inputs$colorPalette <- c("red", rep("grey70", 2), "blue")
    if (x$fitResult$coef["b"] >= x$fitResult$coef["m1"] && 
        x$fitResult$coef["b"] >= x$fitResult$coef["m2"]) {
      inputs$colorPalette <- rev(inputs$colorPalette)
    }
    # TODO: what to do in the 'undefined' case - agonist+antagonist or both flat?
  }

  #TODO include the name of the `color` to be used in the legend of `plotResponseSurface()`
  
  if (color == "z-score") {
    boundary <- sd(x$offAxisTable[["z.score"]])
    inputs$colorBy <- x$offAxisTable[, c("d1", "d2", "z.score")]
    if (!exists("breaks", inputs)) inputs$breaks <- c(-Inf, -boundary, 0, boundary, Inf)
    if (!exists("main", inputs)) inputs$main <- "Z-scores"
  } else if (color == "maxR") {
    inputs$colorBy <- x$maxR$Ymean[, c("d1", "d2", "R")]
    q <- attr(x$maxR$Ymean, "q")
    #TODO add `q` to the inputs so we can use it in plotResponseSurface
    if (!exists("breaks", inputs)) inputs$breaks <- c(-Inf, -q, 0, q, Inf)
    if (!exists("main", inputs)) inputs$main <- "maxR"
  } else if (color == "occupancy") {
    inputs$colorBy <- x$occupancy
    inputs$colorPalette <- c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5")
    if (!exists("breaks", inputs)) inputs$breaks <- c(0, 0.25, 0.5, 0.75, 1)
    if (!exists("main", inputs)) inputs$main <- "Occupancy rate"
  } else if (color == "effect-size") {
    if(is.null(x$confInt))
      stop("No confidence intervals were calculated")
    if (!exists("main", inputs)) inputs$main <- "Effect size"
    synOut <- x$maxR$Ymean
    names(synOut)[names(synOut) == "call"] <- "synCall"
    effectOut <- x$confInt$offAxis
    names(effectOut)[names(effectOut) == "call"] <- "effectCall"
    effectOut$d1 <- as.numeric(gsub("(.+)_.+", "\\1", rownames(effectOut)))
    effectOut$d2 <- as.numeric(gsub(".+_(.+)", "\\1", rownames(effectOut)))
    x_new <- merge(synOut, effectOut, by = c("d1","d2"))
    inputs$colorBy <- x_new[, c("d1", "d2", "effectCall")]
    if (!exists("breaks", inputs)) inputs$breaks <- seq_len(4)
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
#' @param colorBy String indicating the characteristic to use for coloring ("maxR" or "effect-size"). By default, "maxR".
#' @param ... Further parameters passed to \code{\link{plot.maxR}}
#' @export
contour.ResponseSurface <- function(x, colorBy = "maxR", ...) {
  
  if (!exists("maxR", x))
    stop("maxR statistics were not found.")

  cpdNames <- if (!is.null(x$names)) x$names else c("Compound 1", "Compound 2")
  args <- list(...)
  if (!exists("xlab", args))
    args$xlab <- paste0("Dose (", cpdNames[[1]], ")")
  if (!exists("ylab", args))
    args$ylab <- paste0("Dose (", cpdNames[[2]], ")")
  
  ## Blue is synergy, red is antagonism
  if (!exists("colorPalette", args)) {
    args$colorPalette <- c("red", "white", "blue")
    names(args$colorPalette) <- c("Ant", "None", "Syn")
    if (x$fitResult$coef["b"] >= x$fitResult$coef["m1"] && 
        x$fitResult$coef["b"] >= x$fitResult$coef["m2"]) {
      args$colorPalette <- rev(args$colorPalette)
    }
    # TODO: what to do in the 'undefined' case - agonist+antagonist or both flat?
  }
  
  if (colorBy == "maxR") {
    args$x <- x$maxR
  } else if (colorBy == "effect-size") {
    args$x <- x
  }
  class(args$x) <- c(colorBy, setdiff(class(args$x), c("maxR", "effect-size")))
  do.call(plot, args)
  
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

  ans$occup <- if (!is.null(object$occupancy)) mean(object$occupancy$occupancy) else NULL
  ans$method <- object$method
  object$confInt$cutoff = object$cutoff
  ans$CI = summary(object$confInt)

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
  else if (x$null_model == "bliss" & x$shared_asymptote == TRUE)
    cat("Bliss independence with shared maximal response")
  else if (x$null_model == "bliss" & x$shared_asymptote == FALSE)
    cat("Bliss independence with differing maximal response")
  else if (x$null_model == "loewe2" & x$shared_asymptote == TRUE)
    cat("Standard Loewe Additivity") # FIXME: check
  else if (x$null_model == "loewe2" & x$shared_asymptote == FALSE)
    cat("Alternative generalization of Loewe Additivity")
  else
    cat(x$null_model)

  cat("\n")
  cat("Variance assumption used:", dQuote(x$method))
  if (!is.null(x$occup)) {
    cat("\n")
    cat("Mean occupancy rate:", x$occup)
  }
  cat("\n\n")
  print(x$marginalFit)
  cat("\n")

  if (!is.null(x$meanR)) print(x$meanR)
  if (!is.null(x$maxR)) print(x$maxR)

  if (is.null(x$meanR) & is.null(x$maxR)) {
    cat("\n\n")
    cat("No test statistics were computed.")
  }

  if(!is.null(x$CI)) {
    cat("\nCONFIDENCE INTERVALS\n")
    print(x$CI)
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
         "hsa" = hsa(doseInput, parmInput),
         "bliss" = Blissindependence(doseInput, parmInput),
         "loewe2" = harbronLoewe(doseInput, parmInput))

}
