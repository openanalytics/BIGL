# Inverse 4-parameter log-logistic (response-dose) function
invL4 <- function(y, h, b, m) {
  ((y-b)/(m-y))^(1/h)
}
invL4deriv <- function(y, h, b, m) {
  ((y-b)/(m-y))^(-1/h-1)*(1/h*(1/(m-y)+(y-b)/(m-y)^2))
}

# calculation of d/D in the additivity/synergy equation
doseRatio <- function(response, d, h, b, m, expe, lower, upper) {
  # In case of different asymptotes, response can be outside of [lower, upper]
  # range, this needs to be handled separately (otherwise formula gives NaN)
  if(!d){
    return(0)
  } else if(response <= lower){
    return(if(b < m) Inf else 0)
  } else if (response >= upper){
    return(if(b < m) 0 else Inf)
  } else {
    return(d/(invL4(response, h, b, m)*expe))
  }
}
doseRatioGr = function(response, d, h, b, m, expe, lower, upper){
  if(!d || response < lower || response > upper){
    return(0)
  } else {
    return(d/expe* (-invL4deriv(response, h, b, m)))
  }
}
#' Alternative Loewe generalization
#'
#' @inheritParams generalizedLoewe
harbronLoewe <- function (doseInput, parmInput, asymptotes = 2,
                          startvalues = NULL, newtonRaphson = FALSE, ...) {
  parmInput[c("h1", "h2")] = abs(parmInput[c("h1", "h2")])
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
  
  ## If agonist and antagonist, give an error
  if (!(increasing || decreasing)) {
    stop("Alternative Loewe generalization does not work for diverging marginal curves.")
  }
  expE1 = exp(parm["e1"]);expE2 = exp(parm["e2"])
  lower1 <- min(parm[c("b", "m1")]);upper1 <- max(parm[c("b", "m1")])
  lower2 <- min(parm[c("b", "m2")]);upper2 <- max(parm[c("b", "m2")])
  lower = min(lower1, lower2); upper = max(upper1, upper2)
  solver <- function(dose, par) {
    fun0 <- function(y) {
      res <- doseRatio(y, dose[1], par["h1"], par["b"], par["m1"], expE1, lower1, upper1) +
          doseRatio(y, dose[2], par["h2"], par["b"], par["m2"], expE2, lower2, upper2) - 1
      if(is.finite(res)) res else 1
    }
    gr0 = function(y){
      doseRatioGr(y, dose[1], par["h1"], par["b"], par["m1"], expE1, lower1, upper1) +
        doseRatioGr(y, dose[2], par["h2"], par["b"], par["m2"], expE2, lower2, upper2)
    }
    if(newtonRaphson){
    out = nleqslv(fn = fun0, x = dose[3], jac = gr0,
            control = list(ftol = .Machine$double.eps))$x
    out = if(out <= lower) {lower} else if(out >= upper) {upper} else {out}
    return(out)
    } else {
      parRange <- range(par[c("b", "m1", "m2")])
      if (length(unique(parRange)) == 1) # special case of 2 flat profiles
        return(par["b"])
      uniroot(fun0, parRange, tol = .Machine$double.eps)$root
    }
  }

  ## Remove observations where both drugs are dosed at zero
  allZero <- !rowSums(doseInput != 0)
  dose <- doseInput[!allZero, , drop = FALSE]
  parRange <- range(parm[c("b", "m1", "m2")])
  dose = cbind(dose, if(is.null(startvalues)) {
    mean(parm[c("b", "m1", "m2")])
    } else { startvalues[!allZero]})
  res = if (length(unique(parRange)) == 1){
    rep(parm["b"], nrow(dose))
  } else {# special case of 2 flat profiles
     apply(dose, 1, solver, parm)
  }

  if (all(is.na(res))) stop("Alternative Loewe generalization: no roots found between starting parameters")
  if (anyNA(res)) warning("Alternative Loewe generalization: some roots not found")

  ## Set baseline response for observations where both doses are zero.
  ## Otherwise, use the estimate above.
  if (any(allZero)) {
    rv <- rep(NA, length(allZero))
    rv[!allZero] <- res
    rv[allZero] <- parm["b"]
    res <- rv
  }

  return(res)
}
