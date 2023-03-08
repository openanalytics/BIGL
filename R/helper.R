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


#' R color to RGB (red/green/blue) conversion.
#' @param cname vector of any of the three kinds of R color specifications, i.e., either a color name (as listed by \code{\link{colors}}()), a hexadecimal string of the form "#rrggbb" or "#rrggbbaa" (see \code{\link{rgb}}), or a positive integer i meaning \code{\link{palette}}()[i].
#' @param alpha logical value indicating whether the alpha channel (opacity) values should be returned.
#' @importFrom grDevices rgb col2rgb
col2hex <- function (cname, alpha = FALSE) 
{
  colMat <- col2rgb(cname)
  rgb(red = colMat[1, ]/255, green = colMat[2, ]/255, blue = colMat[3, ]/255)
}
