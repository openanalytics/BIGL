#' Compute off-axis predictions
#'
#' Given a dataframe with dose-response data, this function uses coefficient
#' estimates from the marginal (on-axis) monotherapy model to compute the
#' expected values of response at off-axis dose combinations using a provided
#' null model.
#'
#' @param ... Further arguments that are currently unused
#' @inheritParams fitSurface
#' @return This functions returns a vector with predicted off-axis points
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
predictOffAxis <- function( doseGrid, fitResult, transforms = fitResult$transforms,
                           null_model = c("loewe", "hsa", "bliss", "loewe2"), startvalues = NULL,
                            fit = NULL,...) {
  nm = match.arg(null_model)
  if(is.null(fit)){
    fit = fitOffAxis(doseGrid, fitResult, nm, startvalues)
  }
  out = (if(nm %in% c("loewe")) fit$response else fit)[doseGrid$d1 & doseGrid$d2]
  names(out) = apply(doseGrid[doseGrid$d1 & doseGrid$d2, ], 1, paste, collapse = "_")
  if (!is.null(transforms)) {
      CompositeT <- with(transforms,
                         function(y, args) PowerT(BiolT(y, args), args))
      out <- with(transforms, CompositeT(out, compositeArgs))
  }
  return(out)
}
fitOffAxis = function(doseGrid, fitResult,
                      null_model = c("loewe", "hsa", "bliss", "loewe2"),
                      startvalues = NULL){
  nm = match.arg(null_model)
  switch(nm,
               "loewe" = generalizedLoewe(doseGrid, fitResult$coef,
                                          startvalues = startvalues),
               "hsa" = hsa(doseGrid, fitResult$coef),
               "bliss" = Blissindependence(doseGrid, fitResult$coef),
               "loewe2" = harbronLoewe(doseGrid, fitResult$coef,
                                       startvalues = startvalues)
  )
}
