#' Compute off-axis predictions
#'
#' Given a dataframe with dose-response data, this function uses coefficient
#' estimates from the marginal (on-axis) monotherapy model to compute the
#' expected values of response at off-axis dose combinations using a provided
#' null model.
#'
#' @param doseGrid A dose grid with unique combination of doses
#' @param fit a pre-calculated off-axis fit
#' @param ... Further arguments passed on to the Loewe fitters
#' @inheritParams fitSurface
#' @return This functions returns a named vector with predicted off-axis points
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   transforms <- getTransformations(data)
#'   fitResult <- fitMarginals(data, transforms)
#'     uniqueDoses <- with(data, list("d1" = sort(unique(data$d1)),
#'     "d2" = sort(unique(data$d2))))
#'     doseGrid <- expand.grid(uniqueDoses)
#'   predictOffAxis(fitResult, null_model = "hsa", doseGrid = doseGrid)
predictOffAxis <- function(doseGrid, fitResult, transforms = fitResult$transforms,
                           null_model = c("loewe", "hsa", "bliss", "loewe2"),
                           fit = NULL,...) {
  nm = match.arg(null_model)
  if(is.null(fit)){
    fit = fitOffAxis(doseGrid, fitResult, nm,  ...)
  }
  out = (if(nm %in% c("loewe")) fit$response else fit)[doseGrid$d1 & doseGrid$d2]
  names(out) = getd1d2(doseGrid[doseGrid$d1 & doseGrid$d2, ])
  if (!is.null(transforms)) {
      CompositeT <- with(transforms,
                         function(y, args) PowerT(BiolT(y, args), args))
      out <- with(transforms, CompositeT(out, compositeArgs))
  }
  return(out)
}
fitOffAxis = function(doseGrid, fitResult,
                      null_model = c("loewe", "hsa", "bliss", "loewe2"),
                      ...){
  nm = match.arg(null_model)
  switch(nm,
               "loewe" = generalizedLoewe(doseGrid, fitResult$coef, ...),
               "hsa" = hsa(doseGrid, fitResult$coef),
               "bliss" = Blissindependence(doseGrid, fitResult$coef),
               "loewe2" = harbronLoewe(doseGrid, fitResult$coef,...)
  )
}
