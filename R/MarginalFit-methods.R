#' Summary of \code{MarginalFit} object
#'
#' @param object Output of \code{\link{fitMarginals}} function
#' @param ... Further arguments
#' @export
summary.MarginalFit <- function(object, ...) {

  ans <- list()
  ans$coef <- matrix(object$coef[c("h1", "h2", "m1", "m2", "e1", "e2")],
                     ncol = 2, byrow = TRUE)
  colnames(ans$coef) <- if (!is.null(object$names)) object$names else 
        c("Compound 1", "Compound 2")
  rownames(ans$coef) <- c("Slope", "Maximal response", "log10(EC50)")
  # show log10(EC50) in the summary output
  ans$coef["log10(EC50)", ] <- log10(exp(ans$coef["log10(EC50)", ]))
  # show slope always positive
  ans$coef["Slope", ] <- abs(ans$coef["Slope", ])
  ans$baseline <- object$coef["b"]
  ans$vcov <- object$vcov
  ans$n <- nrow(object$data)
  ans$formula <- object$model$formula
  ans$transforms <- object$transforms

  class(ans) <- "summary.MarginalFit"
  ans
}

#' Print method for summary of \code{MarginalFit} object
#'
#' @param x Summary of \code{MarginalFit} object
#' @param ... Further arguments
#' @export
print.summary.MarginalFit <- function(x, ...) {

  cat("Formula:", x$formula)
  cat("\n")
  cat("Transformations:", ifelse(is.null(x$transforms), "No", "Yes"))
  cat("\n\n")
  print(round(x$coef, 3))
  cat("\n")
  cat("Common baseline at:", round(x$baseline, 3))
  cat("\n")

  if (is.character(x$vcov))
    warning(x$vcov)
  else if (any(eigen(x$vcov)$values < 0))
    warning("Hessian is not positive definite. Estimates might be unstable.")
  else if (any(sqrt(na.omit(diag(x$vcov)[c("e1", "e2")])) > 0.3))
    warning("Variance of EC50 estimates may be large.")

  if (x$n < 12)
    warning("There are less than 12 points in total. Estimates might be imprecise.")

}

#' Coefficients from marginal model estimation
#'
#' @inheritParams summary.MarginalFit
#' @export
coef.MarginalFit <- function(object, ...) {
  object$coef
}

#' Predict values on the dose-response curve
#'
#' @inheritParams summary.MarginalFit
#' @param newdata An optional data frame in which to look for \code{d1} and
#'   \code{d2} variables with which to predict. If omitted, the fitted values
#'   are used. Doses that are passed to this function must correspond to
#'   marginal data, i.e. at least one of the doses must be zero.
#' @export
predict.MarginalFit <- function(object, newdata, ...) {

  if (missing(newdata)) {

    fitted(object)

  } else {

    dose1 <- newdata[,"d1"]
    dose2 <- newdata[,"d2"]
    response <- rep(NA, nrow(newdata))

    ## If combination data is passed, throw an error
    if (any(dose1 > 1e-12 & dose2 > 1e-12))
      stop("Predictions are only available for marginal data.")

    resp <- with(as.list(object$coef), {
      response[dose2 == 0] <- L4(dose1[dose2 == 0], h1, b, m1, e1)
      response[dose1 == 0] <- L4(dose2[dose1 == 0], h2, b, m2, e2)
      response
    })

    if (!is.null(object$transforms$BiolT))
      resp <- with(object$transforms, BiolT(resp, compositeArgs))
    if (!is.null(object$transforms$PowerT))
      resp <- with(object$transforms, PowerT(resp, compositeArgs))

    resp
  }

}

#' Compute fitted values from monotherapy estimation
#'
#' @inheritParams summary.MarginalFit
#' @export
fitted.MarginalFit <- function(object, ...) {
  predict(object, newdata = object$data)
}

#' Estimate of coefficient variance-covariance matrix
#'
#' @inheritParams summary.MarginalFit
#' @export
vcov.MarginalFit <- function(object, ...) {
  object$vcov
}

#' Residuals from marginal model estimation
#'
#' @inheritParams summary.MarginalFit
#' @export
residuals.MarginalFit <- function(object, ...) {

  if (is.null(object$transforms))
    PowerT <- function(x) x
  else
    PowerT <- function(x) object$transforms$PowerT(x, object$transforms$compositeArgs)

  PowerT(object$data$effect) - fitted.MarginalFit(object)
}

#' Residual degrees of freedom in marginal model estimation
#'
#' @inheritParams summary.MarginalFit
#' @export
df.residual.MarginalFit <- function(object, ...) {
  object$df
}

#' Plot monotherapy curve estimates
#'
#' @param x Output of \code{\link{fitMarginals}} function or a
#'   \code{"MarginalFit"} object
#' @inheritParams summary.MarginalFit
#' @param ncol Number of plots per row
#' @param logScale Whether x-axis should be plotted on a logarithmic scale
#' @param smooth Whether to draw a smooth fitted curve (deafult), or 
#'   line segments connecting predicted points only
#' @param dataScale Whether to draw plot on original data scale in case when 
#'   transformations were used for fitting. Default (FALSE) is to plot on the 
#'   \code{coef(x)} scale
#' @return Returns a \code{ggplot} object. It can be consequently modified by
#'   using standard operations on \code{ggplot} objects (if \code{ggplot2}
#'   package is loaded).
#' @importFrom ggplot2 aes_string facet_wrap geom_line geom_point ggplot
#'   scale_x_log10 theme_bw xlab ylab
#' @importFrom scales trans_new
#' @importFrom stats fitted
#' @export
plot.MarginalFit <- function(x, ncol = 2, logScale = TRUE, smooth = TRUE, dataScale = FALSE, ...) {

  data <- as.data.frame(x$data)

  transformF <- function(z, comp) {
    eps <- tapply(z, comp, function(x) min(x[x != 0]))
    z + 0.5 * eps[comp]
  }

  labnames <- c("Response", 
      if (!is.null(x$names)) x$names else c("Compound 1", "Compound 2"))
  if (!is.null(attr(x$data, "orig.colnames"))) {
    labnames <- unlist(attr(x$data, "orig.colnames"))
  }

  ## Reorder the data so that non-zero drug 1 observations are stacked
  ## above non-zero drug 2 observations
  dat <- rbind(data[!data$d2, ], data[!data$d1, ])
  if (!is.null(x$transforms$InvBiolT) & !dataScale) {
    dat$effect <- with(x$transforms,
                       InvBiolT(dat$effect, compositeArgs))
  }
  ## Assign the appropriate Compound 1/2 label to the row scale the doses
  dat$comp <- rep(labnames[2:3], c(sum(!data$d2), sum(!data$d1)))
  # make sure given order is unchanged
  dat$comp <- factor(dat$comp, levels = unique(dat$comp))
  dat$dose <- with(dat, ifelse(!d2, d1, d2))
  if (logScale) dat$dose <- with(dat, transformF(dose, comp))

  # predicted smooth curve
  curveDat <- unique(dat[, c("d1", "d2", "comp", "dose")])
  if (smooth) {
    # make finer grid for smooth prediction lines
    
    gridDat <- rbind(
        expand.grid(d1 = makeGrid(curveDat$d1[curveDat$comp == labnames[2]], log = logScale), d2 = 0, comp = labnames[2]),
        expand.grid(d1 = 0, d2 = makeGrid(curveDat$d2[curveDat$comp == labnames[3]], log = logScale), comp = labnames[3]))

    gridDat$dose <- with(gridDat, ifelse(!d2, d1, d2))
    if (logScale) gridDat$dose <- with(gridDat, transformF(dose, comp))
    curveDat <- gridDat    
  }
  
  curveDat$predicted <- predict(x, curveDat)
  if (!is.null(x$transforms$InvPowerT)) {
    curveDat$predicted <- with(x$transforms,
        InvPowerT(curveDat$predicted, compositeArgs))
  }
  if (!is.null(x$transforms$InvBiolT) & !dataScale) {
    curveDat$predicted <- with(x$transforms,
        InvBiolT(curveDat$predicted, compositeArgs))
  }
  
  # draw a dotted line from 0 to the first non-0 dose
  curveDat$type <- FALSE
  if (any(curveDat$d1 + curveDat$d2 == 0, na.rm = TRUE)) {
    minD1 <- min(dat$d1[dat$d1 != 0], na.rm = TRUE)
    minD2 <- min(dat$d2[dat$d2 != 0], na.rm = TRUE)
    curveDat$type <- (curveDat$d1<=minD1+.Machine$double.eps & curveDat$d2 == 0) | (curveDat$d2<=minD2+.Machine$double.eps & curveDat$d1 == 0)
    # for non-smooth curves, we need to duplicate first non-0 dose to avoid line
    # breakage, as we are actually drawing 2 different lines 
    if (!smooth) {
      auxDat <- rbind(curveDat[abs(curveDat$d1-minD1) < .Machine$double.eps, ], curveDat[abs(curveDat$d2-minD2) < .Machine$double.eps, ])
      auxDat$type <- !auxDat$type
      curveDat <- rbind(curveDat, auxDat)
    }
  }

  p <- ggplot() +
    geom_line(data = curveDat, aes_string(x = "dose", y = "predicted", linetype = "type")) +
    geom_point(data = dat, aes_string(x = "dose", y = "effect")) +
    facet_wrap(~ comp, ncol = ncol, scales = "free_x") +
    scale_linetype_manual(values = c("TRUE" = "dotted", "FALSE" = "solid"), guide = "none") +
    xlab("Dose") + ylab("Effect") + theme_bw()
  if (logScale) p <- p + scale_x_log10()

  p
}

makeGrid <- function(x, n = 100, log = TRUE) {
  xRange <- range(x, na.rm = TRUE)
  minNonZero <- min(x[x != 0], na.rm = TRUE)
  if (log) {
    c(if (xRange[1] == 0) 0 else NULL, 10^seq(from = log10(minNonZero), to = log10(xRange[2]), length.out = n))
  } else 
    c(seq(from = xRange[1], to = xRange[2], length.out = n))
}
