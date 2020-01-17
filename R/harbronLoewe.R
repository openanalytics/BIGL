# Inverse 4-parameter log-logistic (response-dose) function 
invL4 <- function(y, h, b, m, logEC50) {
  
  exp(logEC50) * ((y-b)/(m-y))^(1.0/abs(h))
  
}

# calculation of d/D in the additivity/synergy equation
doseRatio <- function(response, d, h, b, m, e) {
  
  stopifnot(length(response) == length(d))
  D <- rep(0, length(d))
  
  lower <- min(b, m)
  upper <- max(b, m)
  
  # In case of different asymptotes, response can be outside of [lower, upper] 
  # range, this needs to be handled separately (otherwise formula gives NaN)
  validIdx <- response > lower & response < upper
  
  D[validIdx] <- invL4(response[validIdx], h, b, m, e)
  D[response <= lower] <- if (b < m) 0 else Inf
  D[response >= upper] <- if (b < m) Inf else 0
  
  dPart <- ifelse(d == 0, 0, d/D)
  dPart
  
}


#' Alternative Loewe generalization
#' 
#' @inheritParams generalizedLoewe
harbronLoewe <- function (doseInput, parmInput, asymptotes = 2, ...) {

  stopifnot(asymptotes %in% c(1, 2))
  stopifnot(identical(colnames(doseInput), c("d1", "d2")))

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
    
  solver <- function(dose, par) {
    dose <- as.numeric(dose)
    
    fun0 <- function(y) {
      res <- doseRatio(y, dose[1], par["h1"], par["b"], par["m1"], par["e1"]) + 
          doseRatio(y, dose[2], par["h2"], par["b"], par["m2"], par["e2"]) - 1
      res
      ifelse(is.finite(res), res, 1)
    }
    uniroot(fun0, range(par[c("b", "m1", "m2")]), tol = .Machine$double.eps)$root
    
  }
  
  ## Remove observations where both drugs are dosed at zero
  allZero <- !rowSums(doseInput != 0)
  dose <- doseInput[!allZero, , drop = FALSE]
  
  res <- apply(dose, 1, solver, parmInput)

  if (all(is.na(res))) stop("Alternative Loewe generalization: no roots found between starting parameters")
  if (any(is.na(res))) warning("Alternative Loewe generalization: some roots not found")

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
