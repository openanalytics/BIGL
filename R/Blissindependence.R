#' Bliss Independence Model
#' 
#' This function returns fractional response levels for when these are based on
#' Bliss Independence Model.
#' 
#' @inheritParams generalizedLoewe

Blissindependence <- function(doseInput, parmInput, ...) {
  
  pars <- parmInput

  increasing <- pars["m1"] >= pars["b"] && pars["m2"] >= pars["b"]
  decreasing <- pars["m1"] <= pars["b"] && pars["m2"] <= pars["b"]
  
  ## If agonist and antagonist, give an error
  if (!(increasing || decreasing)) {
    stop("Bliss independence does not work for diverging marginal curves.")
  }  
  
  # Calculate prediction mono and rescale to max upper for percentage
  maxRange <- max(abs(pars["m1"]-pars["b"]), abs(pars["m2"]-pars["b"]))
  pred1 <- L4(doseInput[["d1"]], b = pars["h1"], logEC50 = pars["e1"], L = 0, U = 1) * 
      abs(pars["m1"]-pars["b"]) / maxRange
  pred2 <- L4(doseInput[["d2"]], b = pars["h2"], logEC50 = pars["e2"], L = 0, U = 1) *
      abs(pars["m2"]-pars["b"]) / maxRange
  
  # Bliss independence combination - prediction value in percentage
  bliss <- pred1 + pred2 - pred1*pred2
  
  # rescale
  direction <- 1*increasing - 1*decreasing
  pars["b"] + direction * maxRange * bliss
 
}
