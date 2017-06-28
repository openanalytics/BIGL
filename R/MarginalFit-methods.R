#' Summary of \code{MarginalFit} object
#'
#' @param object Output of \code{\link{fitMarginals}} function
#' @param ... Further arguments
#' @export
summary.MarginalFit <- function(object, ...) {

  ans <- list()
  ans$coef <- matrix(object$coef[c("h1", "h2", "m1", "m2", "e1", "e2")],
                     ncol = 2, byrow = TRUE)
  colnames(ans$coef) <- c("Compound 1", "Compound 2")
  rownames(ans$coef) <- c("Slope", "Maximal response", "ln(EC50)")
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

  if (any(eigen(x$vcov)$values < 0))
    warning("Hessian is not positive definite. Estimates might be unstable.")
  else if (any(sqrt(diag(x$vcov)[5:6]) > 0.3))
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

#' Compute fitted values from monotherapy estimation
#'
#' @inheritParams summary.MarginalFit
#' @export
fitted.MarginalFit <- function(object, ...) {

  dose1 <- object$data[,"d1"]
  dose2 <- object$data[,"d2"]
  response <- rep(NA, nrow(object$data))

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
#' @return Returns a \code{ggplot} object. It can be consequently modified by
#'   using standard operations on \code{ggplot} objects (if \code{ggplot2}
#'   package is loaded).
#' @importFrom ggplot2 aes_string facet_wrap geom_line geom_point ggplot
#'   scale_x_log10 theme_bw xlab ylab
#' @importFrom scales trans_new
#' @importFrom stats fitted
#' @export
plot.MarginalFit <- function(x, ncol = 2, logScale = TRUE, ...) {

  data <- x$data

  ## Transformation object for ggplot2
  log10eps_trans <- function() {
    trans_new("log10eps",
              function(z) log10(z + 0.5 * min(z[z != 0])),
              function(z) 10^z - min(10^z)*1/3, domain = c(0, Inf))
  }

  transformF <- function(z, comp) {
    eps <- tapply(z, comp, function(x) min(x[x != 0]))
    z + 0.5 * eps[comp]
  }

  labnames <- c("Response", "Compound 1", "Compound 2")
  if (!is.null(attr(x$data, "orig.colnames"))) {
    labnames <- unlist(attr(x$data, "orig.colnames"))
  }

  ## Predict response for monotherapy observations only
  data <-  data[!data$d1 | !data$d2, ]
  data$predicted <- fitted(x)

  ## Reorder the data so that non-zero drug 1 observations are stacked
  ## above non-zero drug 2 observations
  dat <- rbind(data[!data$d2, ], data[!data$d1, ])
  if (!is.null(x$transforms$PowerT)) {
    dat$effect <- with(x$transforms,
                       PowerT(dat$effect, compositeArgs))
  }
  ## Assign the appropriate Compound 1/2 label to the row scale the doses
  dat$comp <- rep(labnames[2:3], c(sum(!data$d2), sum(!data$d1)))
  dat$dose <- with(dat, ifelse(!d2, d1, d2))
  if (logScale) dat$dose <- with(dat, transformF(dose, comp))

  p <- ggplot(dat, aes_string(x = "dose", y = "effect")) +
    geom_line(aes_string(x = "dose", y = "predicted")) +
    geom_point() +
    facet_wrap(~ comp, ncol = ncol, scales = "free_x") +
    xlab("Dose") + ylab("Effect") + theme_bw()
  if (logScale) p <- p + scale_x_log10()

  p
}
