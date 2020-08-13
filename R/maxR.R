#' Compute maxR statistic for each off-axis dose combination
#'
#' \code{\link{maxR}} computes maxR statistics for each off-axis dose
#' combination given the data provided. It provides a summary with results
#' indicating whether a given point is estimated to be synergetic or
#' antagonistic. These can be based either on normal approximation or a
#' fully bootstrapped distribution of the statistics.
#'
#' @param doseGridOff dose grid for off-axis points
#' @inheritParams fitSurface
#' @inheritParams meanR
#' @importFrom stats rt ecdf
#' @return This function returns a \code{maxR} object with estimates for the
#'   maxR statistical test. \code{maxR} object is essentially a list with
#'   appropriately named elements.
#'
#'   In particular, \code{maxR} object contains \code{"Ymean"} element which is
#'   a summary table of maxR test results for each dose combination. This table
#'   contains mean deviation from the predicted surface, normalized deviation
#'   (\code{"absR"}) as well as a statistical call whether this deviation is
#'   significant. Distributional information on which these calls are made can
#'   be retrieved from the attributes of the \code{"Ymean"} dataframe.
#'
#'   Also, \code{maxR} object contains \code{"Call"} element which indicates the
#'   general direction of the deviation of the observed surface from the null.
#'   This call is based on the strongest local deviation in the \code{"Ymean"}
#'   table. 4 values are available here: \code{"Syn"}, \code{"Ant"},
#'   \code{"None"}, \code{"Undefined"}. If one compound acts as an agonist while
#'   another one is an antagonist, then a deviation from the null is classified
#'   as \code{"Undefined"}. If both compounds act in the same direction, then a
#'   stronger than individual effect is classified as synergy while a weaker
#'   effect would be classified as antagonism.
maxR <- function(data_off, fitResult, transforms = fitResult$transforms,
                 null_model = c("loewe", "hsa", "bliss", "loewe2"), R,
                 CP, reps, nested_bootstrap = FALSE, B.B = NULL,
                 cutoff = 0.95, cl = NULL, B.CP = NULL,
                 method = c("equal", "model", "unequal"), bootStraps,
                 idUnique, n1, doseGridOff, transFun, invTransFun, ...) {
    ## Argument matching
    null_model <- match.arg(null_model)
    method <- match.arg(method)

    FStat <- getMaxRF(data_off, fitResult, method, CP, reps, transforms,
                      null_model, R, n1, transFun = transFun,
                      invTransFun = invTransFun)
    Ymean = data.frame(doseGridOff, R = FStat, absR = abs(FStat),
                       "effect - predicted" = R)
    df0 = fitResult$df

    ## Use normal approximation if B.B is not provided.
    ## Otherwise, bootstrap the procedure to find the distribution.
    if (is.null(B.B)) {

        ## MN: find distribution & overall call & points
        B <- 1e5  # iterations to find distribution of M under null
        sim1 <- abs(matrix(rt(B*n1, df = df0), ncol = n1, byrow = TRUE))

        M <- apply(sim1, 1, max)
        q <- quantile(M, cutoff)
        Rnull <- NULL

    } else {
        Rnull <- vapply(bootStraps, FUN.VALUE = FStat, function(x){
            if(nested_bootstrap){
            paramsBootstrap <- list("data" = x$data, "fitResult" = x$simFit,
                                    "transforms" = transforms,
                                    "null_model" = null_model)
            nestedBootstraps = lapply(integer(B.CP), bootFun, args = paramsBootstrap)
            CP = getCP(nestedBootstraps, null_model, transforms)
        }
        getMaxRF(data = x$data[x$data$d1 & x$data$d2,], fitResult = x$simFit,
                 method = method, CP = CP, reps = reps, transforms = transforms,
                 null_model = null_model, n1 = n1, idUnique = idUnique,
                 respS = x$respS, transFun = transFun, invTransFun = invTransFun)
        })
        M <- apply(abs(Rnull), 2, max)
        q <- quantile(M, cutoff)
    }
    coefFit = fitResult$coef
    eq1  <- coefFit["m1"] == coefFit["b"]
    eq2  <- coefFit["m2"] == coefFit["b"]
    inc1 <- coefFit["m1"] >= coefFit["b"]
    inc2 <- coefFit["m2"] >= coefFit["b"]
    dec1 <- coefFit["m1"] <= coefFit["b"]
    dec2 <- coefFit["m2"] <= coefFit["b"]

    call <- {
        if (max(Ymean$absR) > q) {
            invertCall <- Ymean$R[which.max(Ymean$absR)] < 0
            if (eq1 & eq2) "Undefined"
            else if (inc1 & inc2) c("Syn", "Ant")[1 + invertCall]
            else if (dec1 & dec2) c("Ant", "Syn")[1 + invertCall]
            else "Undefined"
        } else "None"
    }

    Ymean$sign <- Ymean$absR > q
    Ymean$call <- "None"

    ## MN: here make call based on the parameters
    Ymean$call[Ymean$sign] <- if (eq1 & eq2) {
        "Undefined"
    } else if (inc1 & inc2) {
        c("Syn", "Ant")[1+(Ymean$R[Ymean$sign] < 0)]
    } else if (dec1 & dec2) {
        c("Ant", "Syn")[1+(Ymean$R[Ymean$sign] < 0)]
    } else "Undefined"

    attr(Ymean, "df0") <- df0
    attr(Ymean, "cutoff") <- cutoff
    attr(Ymean, "q") <- q
    attr(Ymean, "distr") <- ecdf(M)

    ans <- list("Call" = call,
                "Ymean" = Ymean)
    class(ans) <- append(class(ans), "maxR")
    ans
}