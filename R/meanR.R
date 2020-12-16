#' Compute meanR statistic for the estimated model
#'
#' \code{\link{meanR}} computes the meanR statistic for the provided model
#' and returns the computed F-statistic and the estimated p-value. p-value
#' can be calculated either by assuming an exact distribution or using
#' bootstrapping procedure. In the latter case, null distribution of
#' bootstrapped F-statistics is also returned.
#'
#' @param R Numeric vector containing mean deviation of predicted response
#'   surface from the observed one at each of the off-axis points. If missing,
#'   it will be calculated automatically from output of
#'   \code{\link{predictOffAxis}} function.
#' @param CP Matrix which is part of covariance matrix for the \code{R} argument
#' @param reps Numeric vector containing number of replicates for each off-axis
#'   dose combination. If missing, it will be calculated automatically from output
#'   of \code{\link{predictOffAxis}} function.
#' @param cl If parallel computations are desired, \code{cl} should be a cluster
#'   object created by \code{\link[parallel]{makeCluster}}. If parallel
#'   computing is active, progress reporting messages are not necessarily
#'   ordered as it should be expected.
#' @param bootStraps precomputed bootstrap objects
#' @param paramsBootstrap parameters for the nested bootstrap
#' @param idUnique unique combinations of on-axis points, a character vector
#' @param n1 the number of off-axis points
#' @param data_off data frame with off -axis information
#' @param transFun,invTransFun the transformation and inverse transformation functions for the variance
#' @param ... Further arguments that will be later passed to
#'   \code{\link{generateData}} function during bootstrapping
#' @inheritParams fitSurface
#' @importFrom progress progress_bar
#' @importFrom stats ecdf pf
#' @return This function returns a \code{meanR} object with estimates for the
#'   meanR statistical test. \code{meanR} object is essentially a list with
#'   appropriately named elements.
#'
#'   \code{meanR} object list includes notably the calculated F-statistic,
#'   p-value and degrees of freedom (\code{"n1"} and \code{"df0"} respectively)
#'   used to find the critical value of the F-distribution under the null.
#'
#'   If \code{\link{meanR}} test is run with bootstrapping, then p-value
#'   estimate is based on bootstrapped null distribution of test statistic and an
#'   additional element \code{"FDist"} (of class \code{"ecdf"}) is returned.
meanR <- function(data_off, fitResult, transforms = fitResult$transforms,
                  null_model = c("loewe", "hsa", "bliss", "loewe2"), R, CP, reps,
                  nested_bootstrap = FALSE, B.B = NULL, B.CP = NULL,
                  cl = NULL, method = c("equal", "model", "unequal"),
                  bootStraps, paramsBootstrap, idUnique, n1, transFun,
                  invTransFun, ...) {

    ## Argument matching
    null_model <- match.arg(null_model)
    method <- match.arg(method)
    df0 = fitResult$df

    if (all(reps == 1) && method %in% c("model", "unequal")) {
        stop("Replicates are required when choosing the method 'model' or 'unequal'")
    }
    FStat <- getMeanRF(data_off, fitResult, method, CP, reps, transforms, null_model,
                       R, n1, idUnique, transFun = transFun, invTransFun = invTransFun)
    if (is.null(B.B)) {
        ans <- list("FStat" = FStat,
                    "p.value" = pf(FStat, n1, df0, lower.tail = FALSE),
                    "n1" = n1, "df0" = df0)
        class(ans) <- append("meanR", class(ans))
        return(ans)
    }
    FStatb <- vapply(bootStraps, FUN.VALUE = FStat, function(x) {
        if(nested_bootstrap){
            paramsBootstrap <- list("data"  = x$data,
                                    "fitResult" = x$simFit, "transforms" = transforms,
                                    "null_model" = null_model)
            nestedBootstraps = lapply(integer(B.CP), bootFun, args = paramsBootstrap)
            CP = getCP(nestedBootstraps, null_model, transforms)
        }
        getMeanRF(data = x$data[x$data$d1 & x$data$d2,], fitResult = x$simFit, method = method, CP = CP,
                  reps = reps, transforms = transforms, null_model = null_model,
                  n1 = n1, idUnique = idUnique, respS = x$respS,
                  transFun = transFun, invTransFun = invTransFun)
    })
    pvalb <- mean(FStatb >= FStat)
    ans <- list("FStat" = FStat, "FDist" = ecdf(FStatb), "p.value" = pvalb,
                "n1" = n1, "df0" = df0)
    class(ans) <- append("meanR", class(ans))
    ans
}