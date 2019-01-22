#' Fit response surface model and compute meanR and maxR statistics
#'
#' This function computes predictions for off-axis dose combinations according
#' to the BIGL or HSA null model and, if required, computes appropriate meanR
#' and maxR statistics. Function requires as input dose-response dataframe and
#' output of \code{\link{fitMarginals}} containing estimates for the monotherapy
#' model. If transformation functions were used in monotherapy estimation, these
#' should also be provided.
#'
#' Please see the example vignette \code{vignette("analysis", package = "BIGL")}
#' and the report "Lack of fit test for detecting synergy" included in the
#' \code{papers} folder for further details on the test statistics used:
#' \code{system.file("papers", "newStatistics.pdf", package = "BIGL")}
#'
#' @param data Dose-response dataframe.
#' @param fitResult Monotherapy (on-axis) model fit, e.g. produced by
#'   \code{\link{fitMarginals}}. It has to be a \code{"MarginalFit"} object or a
#'   list containing \code{df}, \code{sigma}, \code{coef},
#'   \code{shared_asymptote} and \code{method} elements for, respectively,
#'   marginal model degrees of freedom, residual standard deviation, named
#'   vector of coefficient estimates, logical value of whether shared asymptote
#'   is imposed and method for estimating marginal models during bootstrapping
#'   (see \code{\link{fitMarginals}}). If biological and power transformations
#'   were used in marginal model estimation, \code{fitResult} should contain
#'   \code{transforms} elements with these transformations. Alternatively, these
#'   can also be specified via \code{transforms} argument.
#' @param effect Name of the response column in the data ("effect")
#' @param d1 Name of the column with doses of the first compound ("d1")
#' @param d2 Name of the column with doses of the second compound ("d2")
#' @param transforms Transformation functions. If non-null, \code{transforms} is
#'   a list containing 5 elements, namely biological and power transformations
#'   along with their inverse functions and \code{compositeArgs} which is a list
#'   with argument values shared across the 4 functions. See vignette for more
#'   information.
#' @param statistic Which statistics should be computed. This argument can take
#'   one of the values from \code{c("none", "meanR", "maxR", "both")}.
#' @param cutoff Cut-off to use in maxR procedure for declaring non-additivity
#'   (default is 0.95).
#' @param B.CP Number of bootstrap iterations to use for CP matrix estimation
#' @param B.B Number of iterations to use in bootstrapping null distribution for
#'   either meanR or maxR statistics.
#' @param nested_bootstrap When statistics are calculated, if
#'   \code{nested_bootstrap = TRUE}, \code{CP} matrix is recalculated at each
#'   bootstrap iteration of \code{B.B} using \code{B.CP} iterations. Using such
#'   nested bootstrap may however significantly increase computational time. If
#'   \code{nested_bootstrap = FALSE}, \code{CP} bootstrapped data reuses
#'   \code{CP} matrix calculated from the original data.
#' @param error Type of error for resampling in the bootstrapping procedure.
#'   This argument will be passed to \code{\link{generateData}}. If \code{error
#'   = 4} (default), the error terms for generating distribution of the null
#'   will be resampled from the vector specified in \code{sampling}. If
#'   \code{error = 1}, normal errors are added. If \code{error = 2}, errors are
#'   sampled from a mixture of two normal distributions. If \code{error = 3},
#'   errors are generated from a rescaled chi-square distribution.
#' @param sampling_errors Sampling vector to resample errors from. Used only if
#'   \code{error} is 4 and is passed as argument to \code{\link{generateData}}.
#'   If \code{sampling_errors = NULL} (default), mean residuals at off-axis
#'   points between observed and predicted response are taken.
#' @param parallel Whether parallel computing should be used for bootstrap. This
#'   parameter can take either integer value to specify the number of threads to
#'   be used or logical \code{TRUE/FALSE}. If \code{parallel = TRUE}, then
#'   \code{max(1, detectCores()-1)} is set to be the number of threads. If
#'   \code{parallel = FALSE}, then a single thread is used and cluster object
#'   is not created.
#' @param CP Prediction covariance matrix. If not specified, it will be estimated
#'   by bootstrap using \code{B.CP} iterations.
#' @inheritParams generateData
#' @importFrom parallel makeCluster clusterSetRNGStream detectCores stopCluster
#' @return This function returns a \code{ResponseSurface} object with estimates
#'   of the predicted surface. \code{ResponseSurface} object is essentially a
#'   list with appropriately named elements.
#'
#'   Elements of the list include input data, monotherapy model coefficients and
#'   transformation functions, null model used to construct the surface as well
#'   as estimated CP matrix (see \code{\link{CPBootstrap}}), occupancy level at
#'   each dose combination according to the generalized Loewe model and
#'   \code{"offAxisTable"} element which contains observed and predicted effects
#'   as well as estimated z-scores for each dose combination.
#'
#'   If statistical testing was done, returned object contains \code{"meanR"}
#'   and \code{"maxR"} elements with output from \code{\link{meanR}} and
#'   \code{\link{maxR}} respectively.
#' @examples
#' \dontrun{
#'   data <- subset(directAntivirals, experiment == 4)
#'   ## Data should contain d1, d2 and effect columns
#'   transforms <- list("PowerT" = function(x, args) with(args, log(x)),
#'                      "InvPowerT" = function(y, args) with(args, exp(y)),
#'                      "BiolT" = function(x, args) with(args, N0 * exp(x * time.hours)),
#'                      "InvBiolT" = function(y, args) with(args, 1/time.hours * log(y/N0)),
#'                      "compositeArgs" = list(N0 = 1, time.hours = 72))
#'   fitResult <- fitMarginals(data, transforms)
#'   surf <- fitSurface(data, fitResult, statistic = "meanR")
#'   summary(surf)
#' }
#' @export
fitSurface <- function(data, fitResult,
                       transforms = fitResult$transforms,
                       null_model = c("loewe", "hsa", "bliss"),
                       effect = "effect", d1 = "d1", d2 = "d2",
                       statistic = c("none", "meanR", "maxR", "both"),
                       CP = NULL, B.CP = 50, B.B = NULL, nested_bootstrap = FALSE,
                       error = 4, sampling_errors = NULL, wild_bootstrap = FALSE,
                       cutoff = 0.95, parallel = TRUE) {

  ## Argument matching
  null_model <- match.arg(null_model)
  statistic <- match.arg(statistic)

  ## Verify column names of input dataframe
  if (!all(c(effect, d1, d2) %in% colnames(data)))
    stop("effect, d1 and d2 arguments must be column names of data")
  id <- match(c(effect, d1, d2), colnames(data))
  colnames(data)[id] <- c("effect", "d1", "d2")

  sigma0 <- fitResult$sigma
  df0 <- fitResult$df
  MSE0 <- sigma0^2
  coefFit <- fitResult$coef

  paramsEstimate <- list("data" = data,
                         "fitResult" = fitResult,
                         "coefFit" = coefFit,
                         "null_model" = null_model,
                         "transforms" = transforms)

  respS <- do.call(predictOffAxis, paramsEstimate)
  offaxisZScores <- with(respS$offaxisZTable,
                         (effect - predicted) / sigma0)
  offAxisTable <- cbind(respS$offaxisZTable, "z.score" = offaxisZScores)

  ## Setup parallel computation
  if ((is.logical(parallel) & parallel) | is.numeric(parallel)) {
    nCores <- ifelse(is.logical(parallel),
                     max(1, detectCores() - 1),
                     parallel)
    clusterObj <- makeCluster(nCores, outfile="")
    clusterSetRNGStream(clusterObj)
  } else {
    clusterObj <- NULL
  }


### Computation of MeanR/MaxR statistics
  ## Bootstrap sampling vector
  if (is.null(sampling_errors)) {
    ## Ensure errors are generated from transformed data if applicable
    dataT <- data[, c("d1", "d2", "effect")]
    if (!is.null(transforms)) {
      dataT$effect <- with(transforms,
                           PowerT(dataT$effect, compositeArgs))
    }

    mean_effects <- aggregate(effect ~ d1 + d2, dataT, mean)
    Total <- merge(dataT, mean_effects, by = c("d1", "d2"))
    colnames(Total) <- c("d1", "d2", "effect", "meaneffect")
    sampling_errors <- Total$effect - Total$meaneffect
  }

  ## NB: mean responses are taken
  Ymean <- aggregate(effect - predicted ~ d1 + d2,
                     respS$offaxisZTable, mean)
  R <- Ymean[["effect - predicted"]]
  reps <- aggregate(effect ~ d1 + d2, respS$offaxisZTable, length)$effect
  ## Covariance matrix of the predictions
  stopifnot(length(R) == length(reps))

  paramsBootstrap <- list("B.B" = B.B, "B.CP" = B.CP,
                          "error" = error, "sampling_errors" = sampling_errors,
                          "nested_bootstrap" = nested_bootstrap,
                          "wild_bootstrap" = wild_bootstrap,
                          "cutoff" = cutoff, "Ymean" = Ymean,
                          "reps" = reps, "R" = R,
                          "clusterObj" = clusterObj)

  ## If not provided, compute prediction covariance matrix by bootstrap
  if (is.null(CP)) CP <- do.call(CPBootstrap, c(paramsBootstrap, paramsEstimate))
  paramsBootstrap$CP <- CP

  retObj <- list("data" = data,
                 "fitResult" = fitResult,
                 "transforms" = transforms,
                 "null_model" = null_model,
                 "offAxisTable" = offAxisTable,
                 "occupancy" = respS$occupancy,
                 "CP" = CP)

  statObj <- NULL
  if (statistic %in% c("meanR", "both"))
    statObj <- c(statObj, list("meanR" = do.call(meanR, c(paramsBootstrap, paramsEstimate))))
  if (statistic %in% c("maxR", "both"))
    statObj <- c(statObj, list("maxR" = do.call(maxR, c(paramsBootstrap, paramsEstimate))))

  retObj <- c(retObj, statObj)
  if (!is.null(clusterObj)) stopCluster(clusterObj)

  class(retObj) <- append(class(retObj), "ResponseSurface")
  return(retObj)

}
