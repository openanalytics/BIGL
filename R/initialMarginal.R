#' Estimate initial values for fitting marginal dose-response curves
#'
#' This is a wrapper function which, when a dose-response dataframe is provided,
#' returns start value estimates for both compounds that could be supplied to
#' \code{\link{fitMarginals}} function. This function is also used by
#' \code{\link{fitMarginals}} if no initials values were supplied.
#'
#' Note that this function returns \code{e1} and \code{2} which are
#' log-transformed inflection points for respective compounds.
#'
#' @param ... Further parameters that are currently not used
#' @return Named vector with parameter estimates. Parameter names are consistent
#'   with parameter names in \code{\link{fitMarginals}}. \code{h1} and \code{h2}
#'   are Hill's slope coefficients for each of the compounds, \code{m1} and
#'   \code{m2} are their maximal response levels whereas \code{b} is the shared
#'   baseline. Lastly, \code{e1} and \code{e2} are log-transformed EC50 values.
#' @note Returns starting value for e = log(EC50).
#' @inheritParams fitMarginals
#' @inheritParams fitSurface
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   transforms <- getTransformations(data)
#'   initialMarginal(data, transforms)
initialMarginal <- function(data, transforms = NULL, ...) {

  ## Extract marginal data for each of the compounds
  df1 <- with(data[data$d2 == 0,],
              data.frame("dose" = d1, "effect" = effect))
  df2 <- with(data[data$d1 == 0,],
              data.frame("dose" = d2, "effect" = effect))

  sg1 <- GetStartGuess(df1, transforms)
  sg2 <- GetStartGuess(df2, transforms)

  mean.sg.cc <- mean(c(sg1["cc"], sg2["cc"]))
  mean.sg.dd <- mean(c(sg1["dd"], sg2["dd"]))

  pg <- c(h1 = sg1[["bb"]], h2 = sg2[["bb"]],
          b = mean.sg.cc,
          m1 = sg1[["dd"]], m2 = sg2[["dd"]],
          e1 = sg1[["ee"]], e2 = sg2[["ee"]])
  pg[is.na(pg)] <- 1e-05

  pg
}

#' Estimate initial values for dose-response curve fit
#'
#' @param df Dose-response dataframe containing \code{"dose"} and
#'   \code{"effect"} columns
#' @inheritParams fitSurface
#' @importFrom MASS rlm
#' @importFrom stats lm
GetStartGuess <- function(df, transforms = NULL) {

  x <- df$dose
  resp <- df$effect

  ind <- !is.na(resp)
  resp <- resp[ind]
  x <- x[ind]
  ox <- order(x)
  resp <- resp[ox]
  x <- x[ox]

  ## Work out whether response is increasing or decreasing with increasing
  ## dose of the drug
  if (!is.null(transforms$PowerT))
    resp <- with(transforms, PowerT(resp, compositeArgs))
  mod <- lm((resp - resp[1]) ~ I(seq_along(x) - 1) - 1)

  ## Based on positive/negative relationship, assign asymptotes accordingly
  if (coef(mod) > 0){
    ## Maximum response at x = Inf
    c0 <- min(resp)
    d0 <- max(resp)
  } else {
    ## Maximum response at x = 0
    c0 <- max(resp)
    d0 <- min(resp)
  }

  ## Given previous min/max values, response variable is squeezed into 0/1 so
  ## as to treat 0/1 as asymptotes. It is also slightly shrunk to allow
  ## inclusion of 1 and 0 and then logit-transformed.
  ## r <- (resp - min(x) * 0.99)/(max(resp * 1.01) - min(resp) * 0.99)
  r <- (resp - c0)/(d0 - c0)
  r <- r*0.999 + 0.0005
  r <- log(abs(r/(1 - r)))

  ## bb = Hill coefficient / slope
  ## cc = baseline response
  ## dd = asymptote
  ## ee = inflection point / EC50

  d <- log(x)
  ## ee represents the midpoint of the logistic curve that is fit to the
  ## monotherapy curve. bb is the slope coefficient of the curve
  theFit <- suppressWarnings(
    rlm(r ~ d, subset = !is.na(r) &
                        is.finite(r) &
                        !is.na(d) &
                        is.finite(d))$coef)
  bb <- theFit[[2]]
  ee <- -theFit[[1]]/bb
  ## If the midpoint is beyond range, set the slope coefficient to 1.
  if (exp(ee) > max(x)) {
    theFit <- suppressWarnings(
      rlm(r ~ offset(d), subset = !is.na(r) &
                                  is.finite(r) &
                                  !is.na(d) &
                                  is.finite(d))$coef)
    bb <- 1
    ee <- -theFit[[1]]
  }

  ## Use previous min/max for baseline and asymptote values
  cc <- c0  # median(resp[x == min(x)])
  dd <- d0  # median(resp[x == max(x)])

  ## If transformation functions were provided, baseline and asymptote
  ## parameters are accordingly transformed.
  if (!is.null(transforms$PowerT)) {
    cc <- with(transforms, InvPowerT(cc, compositeArgs))
    dd <- with(transforms, InvPowerT(dd, compositeArgs))
  }
  if (!is.null(transforms)) {
    eps <- min(abs(resp)[abs(resp) > 0])/2
    ## Transform baseline response
    cc0 <- with(transforms, InvBiolT(cc, compositeArgs))
    if (!is.finite(cc0)) {
      cc <- with(transforms,
                 InvBiolT(cc + eps, compositeArgs))
    } else cc <- cc0
    ## Transform asymptotic response
    dd0 <- with(transforms, InvBiolT(dd, compositeArgs))
    if (!is.finite(dd0)) {
      dd <- with(transforms,
                 InvBiolT(dd + eps, compositeArgs))
    } else dd <- dd0
  }

  return(c("bb" = abs(bb),
           "cc" = cc,
           "dd" = dd,
           "ee" = ee))
}
