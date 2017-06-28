#' Apply two-parameter Box-Cox transformation
#'
#' @param y Numeric vector
#' @param lambda Power parameter in power transform
#' @param alpha Shift paramater in 2-parameter power transform. Defaults to
#'   \code{0} which implies a 1-parameter Box-Cox transform.
#' @return Power-transformed data
boxcox.transformation <- function(y, lambda, alpha = 0){

    if (lambda == 0) {
        log(y + alpha)
    } else {
      ((y + alpha)^lambda - 1)/lambda
    }

}

#' Summarize data by factor
#'
#' @param value data to sumamrize
#' @param fac factor to summarize by
#' @importFrom stats median na.omit sd
get.summ.data <- function(value, fac){

  data <- data.frame("value" = value, "rep" = fac)
  aggs <- aggregate(data$value,
                    by = list(data$rep),
                    FUN = function(x) c("median" = median(x, na.rm = TRUE),
                                        "mean" = mean(x, na.rm = TRUE),
                                        "sd" = sd(x, na.rm = TRUE),
                                        "N" = length(x)))

  aggs <- data.frame("Rep" = aggs[["Group.1"]], aggs$x)
  data <- na.omit(data)	   	        # Remove observations with n_i < 2
  return(aggs)
}

#' Return absolute t-value, used in optimization call in
#' \code{\link{optim.boxcox}}
#'
#' @param value data
#' @param fac factor
#' @param lambda box-cox parameter
#' @param zero.add2 2nd box-cox parameter
#' @importFrom robustbase lmrob
#' @importFrom graphics plot abline
get.abs_tval <- function(value, fac, lambda, zero.add2 = 0){

  value <- boxcox.transformation(value, lambda, zero.add2)
  data <- get.summ.data(value, fac)

  if (all(data$mean > 0))
    absmin <- 0
  else
    absmin <- abs(min(data$mean)) + 1e-05 * diff(range(data$mean))

  ## We check whether variance and mean in power-transformed data
  ## appear to be related
  m <- lmrob(log10(data$sd + zero.add2) ~ log10(data$mean + absmin))
  abs.tval <- abs(summary(m)$coefficients[2,3])

  abs.tval
}

#' Find optimal Box-Cox transformation parameters
#'
#' @param value Response variable in the data, e.g. \code{"effect"} column
#' @param fac Factor indicating groups of replicates, e.g.
#'   \code{interaction(d1,d2)}
#' @param shift Whether to use 2-parameter Box-Cox transformation. Input may be
#'   \code{TRUE/FALSE} or a numeric value indicating the shift parameter to use.
#'   If \code{FALSE}, shift parameter is set to zero.
#' @importFrom stats optim optimize
#' @return Numeric vector with power and shift parameter in that order.
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   optim.boxcox(data$effect, interaction(data$d1, data$d2))
#' @export
optim.boxcox <- function(value, fac, shift = FALSE){

  ## If shift parameter is optimized, absmin would be its
  ## lower bound in order to guarantee positivity of the
  ## numerator.
  if (all(value > 0))
    absmin <- 0
  else
    absmin <- abs(min(value)) + 1e-05 * diff(range(value))

  if (isTRUE(shift)) {
    ## Optimize over two parameters
    fn2 <- function(p) get.abs_tval(value, fac, p[1], p[2])

    optim(c(0.5, 0), fn2, method = "L-BFGS-B",
          lower = c(0.1, absmin), upper = c(0.9, Inf))$par
  } else {
    ## Optimize over a single parameter
    if (!is.numeric(shift)) shift <- 0
    fn1 <- function(p) get.abs_tval(value, fac, p, shift)

    c(optimize(fn1, c(0.1, 0.9))[["minimum"]], shift)
  }
}

#' Return a list with transformation functions
#'
#' This function takes in response data from a dose-response model and attempts
#' to find an optimal Box-Cox power transform based on
#' \code{\link{optim.boxcox}} function. It then returns a list of transformation
#' functions which contains this power transform and its inverse which can be
#' subsequently used in \code{\link{fitMarginals}} and \code{\link{fitSurface}}.
#'
#' Additionally, returned list contains biological transform and its inverse
#' based on a simple exponential growth model, especially useful when response
#' data is provided in cell counts. User can additionally provide arguments for
#' these biological transforms where \code{N0} stands for initial cell count and
#' \code{time.hours} indicates number in hours after which response data was
#' measured.
#'
#' @param shift If \code{TRUE} or is a numeric value, then a two-parameter
#'   Box-Cox transformation is assumed. This parameter will be passed on to
#'   \code{\link{optim.boxcox}} function.
#' @param args List with elements that are added to the list of transformation
#'   function and which can be used by these functions. In particular, this
#'   list should be of type \code{args = list("N0" = 1, "time.hours" = 1)} where
#'   \code{N0} and \code{time.hours} are arguments used for the biological
#'   transform.
#' @inheritParams fitSurface
#' @details \code{\link{getTransformations}} relies on
#'   \code{\link{optim.boxcox}} to obtain the optimal Box-Cox transformation
#'   parameters. However, \code{\link{optim.boxcox}} optimizes for the power
#'   parameter only within the interval (0.1, 0.9). Hence, if obtained power
#'   parameter is close to 0.1, then a logarithmic transformation is applied
#'   instead.
#' @return This function returns a list with transformation functions. These
#'   include power transformation (\code{"PowerT"}) and its inverse
#'   (\code{"InvPowerT"}) as well as biological transformation (\code{"BiolT"})
#'   and its inverse (\code{"InvBiolT"}).
#'
#'   Power transformation is a 1-parameter Box-Cox transformation. If
#'   \code{shift = TRUE}, then power transformation is a 2-parameter Box-Cox
#'   transformation. Optimal values for power and shift operators are selected
#'   by means of \code{\link{optim.boxcox}} function.
#'
#'   Biological transformation \code{y = N0 * exp(x * t)} where \code{N0} is the
#'   initial cell count and \code{t} is the incubation time. If response/effect
#'   variable (\code{y}) is given in terms of cell counts, biological
#'   transformation ensures that modelisation is done for the growth rate
#'   instead (\code{x}).
#'
#'   Returned list also contains \code{"compositeArgs"} elements shared by all
#'   the transformation functions. These arguments include initial cell count
#'   (\code{"N0"}) and incubation time (\code{"time.hours"}).
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   getTransformations(data)
#' @export
getTransformations <- function(data, shift = FALSE,
                               args = list("N0" = 1, "time.hours" = 1)) {

  if (any(data$effect <= 0))
    stop("Values below zero are not allowed by power transform.")

  bcPars <- optim.boxcox(data$effect,
                         interaction(data$d1, data$d2), shift)
  if (shift == FALSE) shift <- 0

  if (shift == FALSE) {
    if (abs(bcPars[1] - 0.1) < 1e-2) {
      PowerT <- function(x) log(x)
      InvPowerT <- function(y) exp(y)
    } else {
      PowerT <- function(x) (x^bcPars[1] - 1) / bcPars[1]
      InvPowerT <- function(y) (y*bcPars[1] + 1)^(1/bcPars[1])
    }
  } else {
    if (abs(bcPars[1] - 0.1) < 1e-2) {
      PowerT <- function(x) log(x + bcPars[2])
      InvPowerT <- function(y) exp(y) - bcPars[2]
    } else {
      PowerT <- function(x) ((x + bcPars[2])^bcPars[1] - 1) / bcPars[1]
      InvPowerT <- function(y) (y*bcPars[1] + 1)^(1/bcPars[1]) - bcPars[2]
    }
  }

  transforms <- list(
    "PowerT" = function(x, args) with(args, PowerT(x)),
    "InvPowerT" = function(x, args) with(args, InvPowerT(x)),
    "BiolT" = function(x, args) with(args, N0*exp(x*time.hours)),
    "InvBiolT" = function(y, args) with(args, log(y/N0) / time.hours),
    "compositeArgs" = args
  )

  return(transforms)
}
