#' Generate data from parameters of marginal monotherapy model
#'
#' This function is used to generate data for bootstrapping of the null
#' distribution for various estimates. Optional arguments such as specific
#' choice of sampling vector or corrections for heteroskedasticity can be
#' specified in the function arguments.
#'
#' @param pars Coefficients of the marginal model along with their appropriate
#'   naming scheme. These will typically be estimated using
#'   \code{\link{fitMarginals}}. Futhermore, \code{pars} can simply be a
#'   \code{MarginalFit} object and \code{transforms} object will be
#'   automatically extracted.
#' @param sigma Standard deviation to use for randomly generated error terms. This
#'   argument is unused if \code{error = 4} so that sampling error vector is
#'   provided.
#' @param data Data frame with dose columns \code{("d1", "d2")} to generate the
#'   effect for. Only \code{"d1"} and \code{"d2"} columns of the dose-response
#'   dataframe should be passed to this argument. \code{"effect"} column should
#'   not be passed and if it is, the column will be replaced by simulated data.
#' @param null_model Specified null model for the expected response surface.
#'   Currently, allowed options are \code{"loewe"} for generalized Loewe model,
#'   \code{"hsa"} for Highest Single Agent model, \code{"bliss"} for Bliss additivity,
#'   and \code{"loewe2"} for the alternative Loewe generalization.
#' @param error Type of error for resampling. \code{error = 1} (Default) adds
#'   normal errors to the simulated effects, \code{error = 2} adds errors sampled
#'   from a mixture of two normal distributions, \code{error = 3} generates errors
#'   from a rescaled chi-square distribution. \code{error = 4} will use bootstrap.
#'   Choosing this option, the error terms will be resampled from the vector
#'   specified in \code{sampling_errors}.
#' @param sampling_errors Sampling vector to resample errors from. Used only if
#'   \code{error = 4}.
#' @param wild_bootstrap Whether special bootstrap to correct for
#'   heteroskedasticity should be used. If \code{wild_bootstrap = TRUE}, errors
#'   are generated from \code{sampling_errors} multiplied by a random variable
#'   following Rademacher distribution. Argument is used only if \code{error = 4}.
#' @param ... Further arguments
#' @inheritParams fitSurface
#' @return Dose-response dataframe with generated data including \code{"effect"}
#'   as well as \code{"d1"} and \code{"d2"} columns.
#' @export
#' @examples
#'   coefs <- c("h1" = 1, "h2" = 1.5, "b" = 0,
#'              "m1" = 1, "m2" = 2, "e1" = 0.5, "e2" = 0.1)
#'
#'   ## Dose levels are set to be integers from 0 to 10
#'   generateData(coefs, sigma = 1)
#'
#'   ## Dose levels are taken from existing dataset with d1 and d2 columns
#'   data <- subset(directAntivirals, experiment == 1)
#'   generateData(data = data[, c("d1", "d2")], pars = coefs, sigma = 1)
generateData <- function(pars, sigma, data = NULL,
                         transforms = NULL,
                         null_model = c("loewe", "hsa", "bliss", "loewe2"),
                         error = 1, sampling_errors = NULL,
                         wild_bootstrap = FALSE, ...) {
  
  ## Argument matching
  null_model <- match.arg(null_model)
  
  if (inherits(pars, "MarginalFit")) {
    transforms <- pars$transforms
    pars <- pars$coef
  }
  
  if (is.null(data)) data <- expand.grid("d1" = rep(0:10, each = 2),
                                         "d2" = rep(0:10, each = 2))
  
  if ("effect" %in% colnames(data)) {
    warning("effect column is unneeded for generateData() function and will be dropped.")
    data <- data[, c("d1", "d2")]
  }
  
  ## Use identity transformation if no transform functions are supplied
  if (is.null(transforms)) {
    idF <- function(z, ...) z
    transforms <- list("PowerT" = idF,
                       "InvPowerT" = idF,
                       "BiolT" = idF,
                       "compositeArgs" = NULL)
  }
  
  ySim <- switch(null_model,
                 "loewe" = generalizedLoewe(data, pars, asymptotes = 2)$response,
                 "hsa" = hsa(data[, c("d1", "d2")], pars),
                 "bliss" = Blissindependence(data[, c("d1", "d2")], pars),
                 "loewe2" = harbronLoewe(data[, c("d1", "d2")], pars))
  ySim <- with(transforms,
               PowerT(BiolT(ySim, compositeArgs), compositeArgs))
  
  with(as.list(pars), {
    errors <- switch(as.character(error),
                     ## Normal
                     "1" = { rnorm(length(ySim), 0, sigma) },
                     ## Two normals
                     "2" = { ru <- sample(1:2, replace = TRUE, size = length(ySim))
                     mus <- c(-sigma, sigma)
                     sigmas <- c(sigma/2, sigma/3)
                     rnorm(length(ySim), mus[ru], sigmas[ru]) },
                     ## Distribution with right-tail outliers
                     "3" = { sigma*(rchisq(length(ySim), df=4)-4)/8 },
                     ## Resampling from defined vector
                     "4" = { if (!wild_bootstrap) {
                       errors_test <- sample(sampling_errors, nrow(data), replace = TRUE)
                       errors_test
                     } else {
                       ## Use Rademacher distribution to account for heteroskedasticity
                       errors_test <- sampling_errors*
                         (2*rbinom(length(sampling_errors), size = 1, prob = 0.5)-1)
                       errors_test
                     }},
                     stop("Unavailable error type.")
    )
    
    ySim <- with(transforms, InvPowerT(ySim + errors, compositeArgs))
    data.frame("effect" = ySim, data)
  })
}

#' Estimate CP matrix with bootstrap
#'
#' This function is generally called from within \code{\link{fitSurface}}.
#'
#' @param ... Further parameters that will be passed to \code{\link{generateData}}
#' @inheritParams fitSurface
#' @inheritParams generateData
#' @importFrom stats aggregate var
#' @importFrom parallel clusterApply
#' @return Estimated CP matrix
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 5)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   CPBootstrap(data, fitResult, null_model = "loewe", B.CP = 5)
CPBootstrap <- function(data, fitResult,
                        transforms = fitResult$transforms,
                        null_model = c("loewe", "hsa", "bliss", "loewe2"), B.CP, ...) {

  ## Argument matching
  null_model <- match.arg(null_model)

  if (B.CP < 2) stop("Number of iterations for bootstrapping CP needs to exceed 2.")
  sigma0 <- fitResult$sigma

  pred <- sapply(seq_len(B.CP), function(b) {
    simModel <- simulateNull(data = data, fitResult = fitResult,
                             transforms = transforms,
                             null_model = null_model, ...)
    dataCP <- simModel$data
    fitResultCP <- simModel$fitResult

    predOffAxis <- predictOffAxis(data = dataCP, fitResult = fitResultCP,
                                  null_model = null_model,
                                  transforms = transforms)

    ## If multiple predictions with same x-y coordinates are available,
    ## average them out.
    pred <- aggregate(predicted ~ d1 + d2,
                      predOffAxis$offaxisZTable, mean)$predicted/sigma0

    return(pred)
  })

  CP <- var(t(pred))
  return(CP)
}

#' Data generating function used for constructing null distribution of meanR and
#' maxR statistics
#'
#' This function uses \code{\link{simulateNull}} and simulates all necessary
#' steps to calculate null distribution which will furtherly be used in either
#' \code{\link{meanR}} or \code{\link{maxR}} functions.
#'
#' @param ... Further arguments that will be passed to
#'   \code{\link{generateData}} function
#' @inheritParams fitSurface
#' @inheritParams generateData
bootstrapData <- function(data, fitResult,
                          transforms = fitResult$transforms,
                          null_model = c("loewe", "hsa", "bliss", "loewe2"), ...) {

  ## Argument matching
  null_model <- match.arg(null_model)

  simModel <- simulateNull(data = data, fitResult = fitResult,
                           transforms = transforms,
                           null_model = null_model, ...)
  dataB <- simModel$data
  fitResultB <- simModel$fitResult

  respS <- predictOffAxis(dataB, fitResultB,
                          null_model = null_model,
                          transforms = transforms)

  Rb <- aggregate(effect - predicted ~ d1 + d2,
                  data = respS$offaxisZTable, mean)[[3]]

  repsb <- aggregate(effect ~ d1 + d2,
                     data = respS$offaxisZTable, length)$effect

  ## Covariance matrix of the predictions
  n1b <- length(Rb)
  stopifnot(n1b == length(repsb))

  ## Estimators of the residual variance
  df0b <- fitResultB$df
  MSE0b <- fitResultB$sigma^2
  
  dat_offB  <- dataB[dataB$d1 != 0 & dataB$d2 != 0, ]
  off_varB  <- aggregate(effect ~ d1 + d2, data = dat_offB, var)[["effect"]]
  mse_offb  <- mean(off_varB)
  off_meanB <- aggregate(effect ~ d1 + d2, data = dat_offB, mean)[["effect"]]
  
  linmodB <- lm(off_varB ~ off_meanB)
  PredvarB <- predict(linmodB)

  out <- list("Rb" = Rb, "MSE0b" = MSE0b, "fitResult" = fitResultB,
              "n1b" = n1b, "repsb" = repsb, "Predvarb" = PredvarB,
              "mse_offb" = mse_offb)

  return(out)
}

#' Simulate data from a given null model and monotherapy coefficients
#'
#' @param ... Further parameters that will be passed to
#'   \code{\link{generateData}}
#' @inheritParams fitSurface
#' @inheritParams generateData
#' @return List with \code{data} element containing simulated data and
#'   \code{fitResult} element containing marginal fit on the simulated data.
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   simulateNull(data, fitResult, null_model = "hsa")
simulateNull <- function(data, fitResult,
                         transforms = fitResult$transforms,
                         null_model = c("loewe", "hsa", "bliss", "loewe2"), ...) {

  ## Argument matching
  null_model <- match.arg(null_model)

  method <- fitResult$method
  coefFit0 <- fitResult$coef
  sigma0 <- fitResult$sigma
  model <- fitResult$model

  control <- {
    if (method %in% c("nls", "nlslm"))
      list("maxiter" = 200)
  }

  ## Parameter estimates may at times return an error due to non-convergence. If
  ## necessary, repeat the step until it functions properly and 1000 times at
  ## most.
  counter <- 0
  initPars <- coefFit0
  repeat {

    simData <- generateData(pars = coefFit0, sigma = sigma0,
                            data = data[, c("d1", "d2")],
                            transforms = transforms,
                            null_model = null_model, ...)

    ## In cases where added errors put the response into negative domain, revert
    ## it back to the positive one. Usually, values of such observations tend to
    ## be quite small.
    simData$effect <- abs(simData$effect)

    ## construct a list of arguments, including ... passed to original 
    ## `fitMarginals` call (saved as `extraArgs`)
    paramsMarginal <- list("data" = simData, "method" = method, 
        "start" = initPars, "model" = model, "transforms" = transforms,
        "control" = control)
    if (!is.null(fitResult$extraArgs) && is.list(fitResult$extraArgs))
      # use `modifyList` here, since `control` could be user-defined
      paramsMarginal <- modifyList(paramsMarginal, fitResult$extraArgs)
    
    simFit <- try({
      do.call(fitMarginals, paramsMarginal)
    }, silent = TRUE)

    counter <- counter + 1
    initPars <- NULL
    if (counter > 1000)
      stop(paste("Data simulation process failed. ",
                 "Check that transformation functions correspond ",
                 "to the marginal model."))
    if (!inherits(simFit, "try-error")) break

  }

  return(list("data" = simData,
              "fitResult" = simFit))

}
