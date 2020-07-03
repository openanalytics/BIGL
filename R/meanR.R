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
#' @param ... Further arguments that will be later passed to
#'   \code{\link{generateData}} function during bootstrapping
#' @inheritParams fitSurface
#' @importFrom parallel parSapply
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
#'   estimate is based on boostrapped null distribution of test statistic and an
#'   additional element \code{"FDist"} (of class \code{"ecdf"}) is returned.
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 2)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   CP <- CPBootstrap(data, fitResult, null_model = "loewe", B.CP = 5)
#'   meanR(data, fitResult, null_model = "loewe", CP = CP)
meanR <- function(data, fitResult, transforms = fitResult$transforms,
                  null_model = c("loewe", "hsa", "bliss", "loewe2"), R, CP, reps,
                  nested_bootstrap = FALSE, B.CP = NULL,
                  cl = NULL, method = c("equal", "model", "unequal"),
                  bootStraps, paramsBootstrap, ...) {

    ## Argument matching
    null_model <- match.arg(null_model)
    method <- match.arg(method)
    df0 = fitResult$df

    ## If not supplied, calculate these manually
    if (missing(R) | missing(reps)) {
        respS <- predictOffAxis(data = data, fitResult = fitResult,
                                transforms = transforms, null_model = null_model)
        if (missing(R)) R <- with(respS$offaxisZTable, tapply(effect - predicted, d1d2, mean))
        if (missing(reps)) reps <- with(respS$offaxisZTable,
                                        tapply(effect - predicted, d1d2, length))
    }
    if (all(reps == 1) && method %in% c("model", "unequal")) {
        stop("Replicates are required when choosing the method 'model' or 'unequal'")
    }
    n1 <- length(R)
    FStat <- getMeanRF(data, fitResult, method, CP, reps, transforms, null_model,
                       R, n1)
    if (is.null(bootStraps)) {
        ans <- list("FStat" = FStat,
                    "p.value" = pf(FStat, n1, df0, lower.tail = FALSE),
                    "n1" = n1, "df0" = df0)
        class(ans) <- append("meanR", class(ans))
        return(ans)
    }
    FStatb <- sapply(bootStraps, function(x) {
        if(nested_bootstrap){
            paramsBootstrap <- list("data"  =x$data,
                                    "fitResult" = x$simFit, "transforms" = transforms,
                                    "null_model" = null_model)
            nestedBootstraps = lapply(integer(B.CP), bootFun, args = paramsBootstrap)
            CP = getCP(nestedBootstraps, null_model, transforms)
        }
        getMeanRF(data = x$data, fitResult = x$simFit, method = method, CP = CP,
                  reps = reps, transforms = transforms, null_model = null_model,
                  n1 = n1)
    })
    pvalb <- mean(FStatb >= FStat)
    ans <- list("FStat" = FStat, "FDist" = ecdf(FStatb), "p.value" = pvalb,
                "n1" = n1, "df0" = df0)
    class(ans) <- append("meanR", class(ans))
    ans
}