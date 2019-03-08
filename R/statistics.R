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
                  nested_bootstrap = FALSE, B.B = NULL, B.CP = NULL,
                  cl = NULL,
                  method = c("equal", "model", "unequal"), ...) {

  ## Argument matching
  null_model <- match.arg(null_model)
  method <- match.arg(method)

  ## If not supplied, calculate these manually
  if (missing(R) | missing(reps)) {
    respS <- predictOffAxis(data = data, fitResult = fitResult,
                            transforms = transforms, null_model = null_model)
    Ymean <- aggregate(effect - predicted ~ d1 + d2,
                       respS$offaxisZTable, mean)
    if (missing(R)) R <- Ymean[["effect - predicted"]]
    if (missing(reps)) reps <- aggregate(effect ~ d1 + d2,
                                         data = respS$offaxisZTable, length)[["effect"]]
  }
  if (all(reps == 1) && method %in% c("model", "unequal")) {
    stop("Replicates are required when choosing the method 'model' or 'unequal'")
  }
  
  n1 <- length(R)
  MSE0 <- fitResult$sigma^2
  df0 <- fitResult$df

  dat_off  <- data[data$d1 != 0 & data$d2 != 0, ]
  off_var  <- aggregate(effect ~ d1 + d2, data = dat_off, var)[["effect"]]
  off_mean <- aggregate(effect ~ d1 + d2, data = dat_off, mean)[["effect"]]
  
  mse_off <- switch(method,
      "equal" = MSE0,
      "model" = {
        linmod <- lm(off_var ~ off_mean)
        predict(linmod)
      },
      "unequal" = mean(off_var)
  )

  A <- MSE0*CP + mse_off*diag(1/reps, nrow = n1)
  FStat <- as.numeric(t(R) %*% solve(A) %*% R / n1)
  
  if (is.null(B.B)) {
    ans <- list("FStat" = FStat,
                "p.value" = pf(FStat, n1, df0, lower.tail = FALSE),
                "n1" = n1, "df0" = df0)
    class(ans) <- append("meanR", class(ans))
    return(ans)
  }

  pb <- progress_bar$new(format = "(meanR): [:bar]:percent",
                         total = B.B, width = 60)
  pb$tick(0)

  iterFunction <- function(i) {

    ## Update progress bar
    pb$tick()

    out <- bootstrapData(data = data, fitResult = fitResult,
                         transforms = transforms, null_model = null_model, ...)

    MSE0b <- out$MSE0b
    Rb <- out$Rb
    n1b <- out$n1b
    repsb <- out$repsb
    fitResultb <- out$fitResult
    Predvarb <- out$Predvarb
    mse_offb <- out$mse_offb

    CPb <- CP
    if (nested_bootstrap)
      CPb <- CPBootstrap(data = data, fitResult = fitResultb,
                         transforms = transforms, null_model = null_model,
                         B.CP = B.CP, ...)
    
    mse_offb <- switch(method,
        "equal" = MSE0b,
        "model" = Predvarb,
        "unequal" = mse_offb
    )
      
    Ab <- MSE0b*CPb + mse_offb*diag(1/repsb, nrow = n1b)
    FStatb1 <- t(Rb) %*% solve(Ab) %*% Rb / n1b
    
    return(as.numeric(FStatb1))
  }

  ## Call parallel computation if needed
  if (is.null(cl)) {
    FStatb <- sapply(seq_len(B.B), iterFunction)
  } else {
    FStatb <- parSapply(cl, seq_len(B.B), iterFunction)
  }

  pvalb <- mean(FStatb >= FStat)

  ans <- list("FStat" = FStat,
              "FDist" = ecdf(FStatb),
              "p.value" = pvalb,
              "n1" = n1, "df0" = df0)
  class(ans) <- append("meanR", class(ans))
  ans
}

#' Compute maxR statistic for each off-axis dose combination
#'
#' \code{\link{maxR}} computes maxR statistics for each off-axis dose
#' combination given the data provided. It provides a summary with results
#' indicating whether a given point is estimated to be synergetic or
#' antagonistic. These can be based either on normal approximation or a
#' fully bootstrapped distribution of the statistics.
#'
#' @param Ymean Aggregate summary of off-axis predicted responses. In
#'   particular, it should contain \code{"effect - predicted"} column. If
#'   \code{Ymean} is missing, it will be calculated automatically from output of
#'   \code{\link{predictOffAxis}} function.
#' @inheritParams fitSurface
#' @inheritParams meanR
#' @importFrom parallel parSapply
#' @importFrom progress progress_bar
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
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 2)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   CP <- CPBootstrap(data, fitResult, null_model = "loewe", B.CP = 5)
#'   maxR(data, fitResult, null_model = "loewe", CP = CP)
maxR <- function(data, fitResult, transforms = fitResult$transforms,
                 null_model = c("loewe", "hsa", "bliss", "loewe2"), Ymean, CP, reps,
                 nested_bootstrap = FALSE, B.B = NULL, B.CP = NULL,
                 cutoff = 0.95, cl = NULL, 
                 method = c("equal", "model", "unequal"), ...) {

  ## Argument matching
  null_model <- match.arg(null_model)
  method <- match.arg(method)
  
  ## If not supplied, calculate these manually
  if (missing(reps) | missing(Ymean)) {
    respS <- predictOffAxis(data = data, fitResult = fitResult,
                            transforms = transforms, null_model = null_model)
    if (missing(Ymean)) Ymean <- aggregate(effect - predicted ~ d1 + d2,
                                           data = respS$offaxisZTable, mean)
    if (missing(reps)) reps <- aggregate(effect ~ d1 + d2,
                                         respS$offaxisZTable, length)[["effect"]]
  }
  
  if (all(reps == 1) && method %in% c("model", "unequal")) {
    stop("Replicates are required when choosing the method 'model' or 'unequal'")
  }
  
  MSE0 <- fitResult$sigma^2
  df0 <- fitResult$df
  coefFit <- fitResult$coef
  R <- Ymean[["effect - predicted"]]
  n1 <- length(R)

  dat_off  <- data[data$d1 != 0 & data$d2 != 0, ]
  off_var  <- aggregate(effect ~ d1 + d2, data = dat_off, var)[["effect"]]
  off_mean <- aggregate(effect ~ d1 + d2, data = dat_off, mean)[["effect"]]
  
  mse_off <- switch(method,
      "equal" = MSE0,
      "model" = {
        linmod <- lm(off_var ~ off_mean)
        Predvar <- predict(linmod)
        ifelse(Predvar < 0, 0.00001, Predvar)
      },
      "unequal" = mean(off_var)
  )
  
  A <- MSE0*CP + mse_off*diag(1/reps, nrow = n1)
  
  E <- eigen(A)
  V <- E$values
  Q <- E$vectors
  Amsq <- Q %*% diag(1/sqrt(V)) %*% t(Q)
  
  RStud <- t(R) %*% Amsq
  Ymean$R <- t(RStud)
  Ymean$absR <- abs(Ymean$R)

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

    pb <- progress_bar$new(format = "(maxR): [:bar]:percent",
                           total = B.B, width = 60)
    pb$tick(0)

    iterFunction <- function(i) {

      ## Update progress bar
      pb$tick()

      out <- bootstrapData(data = data, fitResult = fitResult,
                           transforms = transforms, null_model = null_model, ...)
      MSE0b <- out$MSE0b
      Rb <- out$Rb
      n1b <- out$n1b
      repsb <- out$repsb
      fitResultb <- out$fitResult
      Predvarb <- out$Predvarb
      mse_offb <- out$mse_offb

      CPb <- CP
      if (nested_bootstrap)
        CPb <- CPBootstrap(data = data, fitResult = fitResultb,
            transforms = transforms, null_model = null_model,
            B.CP = B.CP, ...)
      
      mse_offb <- switch(method,
          "equal" = MSE0b,
          "model" = ifelse(Predvarb < 0, 0.00001, Predvarb),
          "unequal" = mse_offb
      )
      
      Ab <- MSE0b*CPb + mse_offb*diag(1/repsb, nrow = n1b)

      Eb <- eigen(Ab)
      Vb <- Eb$values
      Qb <- Eb$vectors
      Amsqb <- Qb %*% diag(1/sqrt(Vb)) %*% t(Qb)

      RStudb <- t(Rb) %*% Amsqb
      return(RStudb)
    }

    ## Call parallel computation if needed
    if (is.null(cl)) {
      Rnull <- sapply(seq_len(B.B), iterFunction)
    } else {
      Rnull <- parSapply(cl, seq_len(B.B), iterFunction)
    }

    M <- apply(abs(Rnull), 2, max)
    q <- quantile(M, cutoff)

  }


  call <- {
    if (max(Ymean$absR) > q) {
      invertCall <- Ymean$R[which.max(Ymean$absR)] < 0
      inc1 <- coefFit["m1"] >= coefFit["b"]
      inc2 <- coefFit["m2"] >= coefFit["b"]
      if (inc1 & inc2) c("Syn", "Ant")[1 + invertCall]
      else if (!inc1 & !inc2) c("Ant", "Syn")[1 + invertCall]
      else "Undefined"
    } else "None"
  }

  Ymean$sign <- Ymean$absR > q
  Ymean$call <- "None"

  ## MN: here make call based on the parameters
  Ymean$call[Ymean$sign] <- c("Syn", "Ant")[1+(Ymean$R[Ymean$sign] < 0)]
  if (coefFit["b"] > coefFit["m1"]) {
    Ymean$call[Ymean$sign] <- c("Ant", "Syn")[1+(Ymean$R[Ymean$sign] < 0)]
  }

  attr(Ymean, "df0") <- df0
  attr(Ymean, "cutoff") <- cutoff
  attr(Ymean, "q") <- q
  attr(Ymean, "distr") <- ecdf(M)

  ans <- list("Call" = call,
              "Ymean" = Ymean)
  class(ans) <- append(class(ans), "maxR")
  ans

}
