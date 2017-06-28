#' 4-parameter logistic dose-response function
#'
#' @param dose Dose level
#' @param b Hill's coefficient (slope of the curve)
#' @param L Baseline effect (at zero dose)
#' @param U Asymptote effect (at infinite dose)
#' @param logEC50 Point of inflection (in logarithmic terms)
L4 <- function(dose, b, L, U, logEC50) {

  denum <- 1 + (exp(logEC50) / dose)^(abs(b))
  return(L + (U - L) / denum)

}
