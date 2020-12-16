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
#'   will be resampled from the vector specified in \code{sampling_errors}. If
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
#' @param progressBar A boolean, should progress of bootstraps be shown?
#' @param method What assumption should be used for the variance of on- and 
#'   off-axis points. This argument can take one of the values from 
#'   \code{c("equal", "model", "unequal")}. With the value \code{"equal"} as the
#'   default. \code{"equal"} assumes that both on- and off-axis points have the 
#'   same variance, \code{"unequal"} estimates a different parameter for on- and 
#'   off-axis points and \code{"model"} predicts variance based on the average 
#'   effect of an off-axis point. If no transformations are used the 
#'   \code{"model"} method is recommended. If transformations are used, only the
#'   \code{"equal"} method can be chosen.
#' @param confInt a boolean, should confidence intervals be returned?
#' @param bootRS a boolean, should bootstrapped response surfaces be used in the
#'  calculation of the confidence intervals?
#' @param trans,invtrans the transformation function for the variance and its
#' inverse, possibly as strings
#' @param rescaleResids a boolean indicating whether to rescale residuals,
#' or else normality of the residuals is assumed.
#' @param newtonRaphson A boolean, should Newton-Raphson be used to find Loewe
#' response surfaces? May be faster but also less stable to switch on
#' @param asymptotes Number of asymptotes. It can be either \code{1}
#'   as in standard Loewe model or \code{2} as in generalized Loewe model.
#' @param bootmethod The resampling method to be used in the bootstraps. Defaults to the same as method
#' @inheritParams generateData
#' @importFrom parallel makeCluster clusterSetRNGStream detectCores stopCluster parLapply
#' @importFrom progress progress_bar
#' @importFrom stats aggregate
#' @return This function returns a \code{ResponseSurface} object with estimates
#'   of the predicted surface. \code{ResponseSurface} object is essentially a
#'   list with appropriately named elements.
#'
#'   Elements of the list include input data, monotherapy model coefficients and
#'   transformation functions, null model used to construct the surface as well
#'   as estimated CP matrix, occupancy level at
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
                       null_model = c("loewe", "hsa", "bliss", "loewe2"),
                       effect = "effect", d1 = "d1", d2 = "d2",
                       statistic = c("none", "meanR", "maxR", "both"),
                       CP = NULL, B.CP = 50, B.B = NULL, nested_bootstrap = FALSE,
                       error = 4, sampling_errors = NULL, wild_bootstrap = FALSE,
                       cutoff = 0.95, parallel = FALSE, progressBar = TRUE,
                       method = c("equal", "model", "unequal"), confInt = TRUE,
                       bootRS = TRUE, trans = "identity", rescaleResids = FALSE,
                       invtrans = switch(trans, "identity" = "identity", "log" = "exp"),
                       newtonRaphson = FALSE, asymptotes = 2, bootmethod = method) {

  ## Argument matching
  null_model <- match.arg(null_model)
  statistic <- match.arg(statistic)
  method <- match.arg(method)
  transFun = match.fun(trans); invTransFun = match.fun(invtrans)

  if (method %in% c("model", "unequal") && (!is.null(transforms) || !is.null(fitResult$transforms))) {
    stop("No transformations can be used when choosing the method 'model' or 'unequal'")
  }

  ## Verify column names of input dataframe
  if (!all(c("effect", "d1", "d2") %in% colnames(data)))
    stop("effect, d1 and d2 arguments must be column names of data")
  id <- match(c("effect", "d1", "d2"), colnames(data))
  colnames(data)[id] <- c("effect", "d1", "d2")
  data$d1d2 = apply(data[, c("d1", "d2")], 1, paste, collapse = "_")
  sigma0 <- fitResult$sigma
  df0 <- fitResult$df
  MSE0 <- sigma0^2

  #Off-axis data and predictions
  data_off = with(data, data[d1 & d2, , drop = FALSE])
  uniqueDoses <- with(data, list("d1" = sort(unique(d1)),
     "d2" = sort(unique(d2))))
  doseGrid <- expand.grid(uniqueDoses)
  offAxisFit = fitOffAxis(fitResult, null_model = null_model,
                               doseGrid = doseGrid, newtonRaphson = newtonRaphson)
  offAxisPred = predictOffAxis(fitResult, null_model = null_model,
                                          doseGrid = doseGrid, transforms = transforms,
                                          fit = offAxisFit)
  #Retrieve all off-axis points
  idOffDoseGrid = with(doseGrid, d1 & d2)
  doseGridOff = doseGrid[idOffDoseGrid,]
  d1d2off = apply(doseGridOff, 1, paste, collapse = "_")
  rownames(doseGridOff) = d1d2off
  idUnique = d1d2off[match(data_off$d1d2, d1d2off)]
  offAxisPredAll <- offAxisPred[idUnique]
  offaxisZTable <- cbind(data_off[, c("d1", "d2", "effect", "d1d2"), drop = FALSE],
                         "predicted" = offAxisPredAll)
  if(!is.null(transforms))
    offaxisZTable$effect <- with(transforms,
                      PowerT(offaxisZTable$effect, compositeArgs))
  offAxisTable <- cbind(offaxisZTable,
                        "z.score" = with(offaxisZTable, (effect - predicted) / sigma0))
  if(null_model == "loewe"){
   occupancy = offAxisFit$occupancy
    startvalues = offAxisFit$oc
  } else if(null_model == "loewe2"){
    startvalues = offAxisFit
    occupancy = NULL
  } else {
    occupancy = startvalues = NULL
  }
### Computation of MeanR/MaxR statistics
  ## Bootstrap sampling vector
  if (is.null(sampling_errors)) {
    ## Ensure errors are generated from transformed data if applicable
    dataT <- data[, c("d1", "d2", "effect", "d1d2")]
    if (!is.null(transforms)) {
      dataT$effect <- with(transforms,
                           PowerT(dataT$effect, compositeArgs))
    }
    mean_effects <- aggregate(data = dataT, effect ~ d1d2, mean)
    names(mean_effects) = c("d1d2", "meaneffect")
    Total <- merge(dataT, mean_effects, by = c("d1d2"))
    sampling_errors <- Total$effect - Total$meaneffect
  }

  ## NB: mean responses are taken
  R = with(offaxisZTable, tapply(effect-predicted, d1d2, mean))
  reps <- with(offaxisZTable, tapply(effect-predicted, d1d2, length))
  if (all(reps == 1) && method %in% c("model", "unequal")) {
    stop("Replicates are required when choosing the method 'model' or 'unequal'")
  }
  #Check predicted variances
  if(method == "model"){
    off_mean <- with(data_off, tapply(effect, d1d2, mean))
    off_var = with(data_off, tapply(effect, d1d2, var))
    Coef = lm.fit(transFun(off_var), x = cbind(1, off_mean))$coef
    #Don't allow negative variances
    if(Coef[2]<0){
      #warning("Variance was found to decrease with mean, check mean-variance trend!")
    }
    predVar = invTransFun(Coef[1] + Coef[2]*off_mean)
    if(any(predVar < 0)){
      #warning("Negative variances modelled on real data!\nCheck mean-variance trend with plotMeanVarFit and consider transforming the variance!")
    }
    model = c(Coef, "min" = min(off_var[off_var>0]),
              "max" = max(off_var)) #Store smallest observed variance
  } else model = NULL

  B = if(is.null(B.B)) B.CP else max(B.B, B.CP) #Number of bootstraps
  #If any bootstraps needed, do them first
  if(is.null(CP) && (is.null(B))){
    stop("No covariance matrix supplied, and no bootstraps required.\n Please provide either of both!")
  } else if(!is.null(B)){
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

      #Progess bar
      if(progressBar && !is.null(B)){
      pb <- progress_bar$new(format = "(bootstraps): [:bar]:percent",
                             total = B, width = 60)
      pb$tick(0)
      } else {pb = NULL}
     #Start bootstrapping
       paramsBootstrap <- list("data"  =data,
          "fitResult" = fitResult, "transforms" = transforms,
                              "null_model" = null_model,
          "error" = error, "sampling_errors" = sampling_errors,
                             "wild_bootstrap" = wild_bootstrap, "bootmethod" = bootmethod,
                              "method" = method, "doseGrid" = doseGrid,
          "startvalues" = startvalues, "pb" = pb, "progressBar" = progressBar,
          "model" = model, "means" = Total$meaneffect, "rescaleResids" = rescaleResids,
          "invTransFun" = invTransFun, "newtonRaphson" = newtonRaphson,
          "asymptotes" = asymptotes)

      bootStraps = if(is.null(clusterObj)) {
              lapply(integer(B), bootFun, args = paramsBootstrap)
          } else {
              parLapply(clusterObj, integer(B), bootFun, args = paramsBootstrap)
          }
  } else {bootStraps = clusterObj = NULL}
  ## If not provided, compute prediction covariance matrix by bootstrap
  if (is.null(CP)) CP <- getCP(bootStraps, null_model, transforms,
                               sigma0 = sigma0, doseGrid = doseGrid)
  CP = CP[names(R), names(R)]
  #Calculate test statistics
  paramsStatistics = list("bootStraps" = bootStraps, "CP" = CP, "cutoff" = cutoff,
                          "data_off" = data_off, "fitResult" = fitResult,
                          "null_model" = null_model, "transforms" = transforms,
                          "doseGrid" = doseGrid, "reps" = reps, "R" = R,
                          "idUnique" = idUnique, "B.B" = B.B,
                          "Total" = Total, "n1" = length(R), "method" = method,
                          "respS" = offAxisPredAll, "bootRS" = bootRS,
                          "doseGridOff" = doseGridOff[names(R),], "transFun" = transFun,
                          "invTransFun" = invTransFun, "model" = model, "rescaleResids" = rescaleResids)
  statObj <- NULL
  if (statistic %in% c("meanR", "both"))
      statObj <- c(statObj, list("meanR" = do.call(meanR, paramsStatistics)))
  if (statistic %in% c("maxR", "both"))
      statObj <- c(statObj, list("maxR" = do.call(maxR, paramsStatistics)))
  if(confInt && is.null(B.B) && is.null(B.CP)){
    warning("Confidence intervals only available with the bootstrap")
    confInt = FALSE
  }
  if(confInt)
      statObj <- c(statObj, list("confInt" = do.call(bootConfInt, paramsStatistics)))

  retObj <- c(list("data" = data, "fitResult" = fitResult,
                   "transforms" = transforms, "null_model" = null_model,
                   "method" = method, "offAxisTable" = offAxisTable, "asymptotes" = asymptotes,
                   "occupancy" = occupancy, "CP" = CP, "cutoff" = cutoff), statObj)
  if (!is.null(clusterObj)) stopCluster(clusterObj)
  # add compound names from marginal fit
  retObj$names <- fitResult$names
  class(retObj) <- append(class(retObj), "ResponseSurface")
  return(retObj)
}
