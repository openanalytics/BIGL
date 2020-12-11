#' Compute combined predicted response from drug doses according to standard or
#' generalized Loewe model.
#'
#' @param doseInput Dose-response dataframe containing \code{"d1"} and
#'   \code{"d2"} columns
#' @param parmInput Numeric vector or list with appropriately named
#'   parameter inputs. Typically, it will be coefficients from a
#'   \code{MarginalFit} object.
#' @inheritParams fitSurface
#' @param newtonRaphson a boolean, is Newton raphson used for finding the
#' response surface? May be faster but also less stable
#' @param ... Further arguments that are currently unused
#' @inheritParams simulateNull
#' @importFrom nleqslv nleqslv
#' @importFrom stats uniroot
generalizedLoewe <- function (doseInput, parmInput, asymptotes = 2,
                              startvalues = NULL, newtonRaphson = FALSE, ...) {
  parmInput[c("h1", "h2")] = abs(parmInput[c("h1", "h2")])
  ## Need good accuracy here: solve for -logit(o)
  solver <- function(dose, par){
    fun0 <- function(x){
      res <- exp(dose[1] + x/par["h1"]) + exp(dose[2] + x/par["h2"]) - 1
      if(is.finite(res)) res else 1
    }
    gr0 = function(x){
      exp(dose[1] + x/par["h1"])/par["h1"] +
        exp(dose[2] + x/par["h2"])/par["h2"]
    }
    if(newtonRaphson)
      nleqslv(fn = fun0, x = dose[3], jac = gr0)$x
    else
      uniroot(fun0, c(-5000, 5000), tol = .Machine$double.eps)$root
  }

  ## Remove observations where both drugs are dosed at zero
  allZero <- !rowSums(doseInput != 0)
  dose <- doseInput[!allZero,, drop = FALSE]
  logDoseMinE = log(dose) - rep(parmInput[c("e1", "e2")], each = nrow(dose))

  ## In case of a single asymptote, use an artificial one for the second drug
  ## equal to the first one.
  if (asymptotes == 1) {
    parm <- c(parmInput[1:4], parmInput[4], parmInput[5:6])
    names(parm)[4:5] <- c("m1", "m2")
  } else {
    parm <- parmInput
  }

  increasing <- parm["m1"] >= parm["b"] && parm["m2"] >= parm["b"]
  decreasing <- parm["m1"] <= parm["b"] && parm["m2"] <= parm["b"]
  
  ## If agonist and antagonist, give a warning
  if (!(increasing || decreasing)) {
    warning("Marginal curves are diverging. The synergy/antagonism calls may be reversed")
  }

  logDoseMinE = cbind(logDoseMinE,
                      "start" = if(is.null(startvalues)) integer(nrow(logDoseMinE)) else
                        startvalues)
  ## For each combination of Compound 1 and Compound 2, find transformed
  ## occupancy rate, i.e. -logit(o*), which satisfies Loewe model equation.
  oc <- apply(logDoseMinE, 1, solver, parmInput)
  if (all(is.na(oc))) stop("genLoewe: no roots found between starting parameters")
  if (any(is.na(oc))) warning("genLoewe: some roots not found")

  LogOccupancy <- function(d, e, o, h) log(d) - e + abs(1/h) * o
  xf <- with(as.list(parm), {
    logO1 <- LogOccupancy(dose[["d1"]], e1, oc, h1)
    logO2 <- LogOccupancy(dose[["d2"]], e2, oc, h2)
    rvInt <- (m1 - b) * exp(logO1) + (m2 - b) * exp(logO2)
    rv <- b + rvInt/(exp(oc) + 1)
    rv
  })

  ## Set baseline response for observations where both doses are zero.
  ## Otherwise, use the estimate above.
  if (any(allZero)) {
    rv <- rep(NA, length(allZero))
    rv[!allZero] <- xf
    rv[allZero] <- parm["b"]
    xf <- rv
  }

  return(list("response" = xf, "oc" = oc,
              "occupancy" = cbind(dose, "occupancy" = 1/(exp(oc)+1))))
}
