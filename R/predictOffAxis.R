#' Compute off-axis predictions
#'
#' Given a dataframe with dose-response data, this function uses coefficient
#' estimates from the marginal (on-axis) monotherapy model to compute the
#' expected values of response at off-axis dose combinations using a provided
#' null model.
#'
#' @param ... Further arguments that are currently unused
#' @inheritParams fitSurface
#' @return This functions returns a list with 3 elements.
#'
#'   \code{"offaxisZTable"} is a dataframe containing dose levels, observed
#'   effects and effects predicted according to the specified null model. This
#'   dataframe also contains replicates, if there are any.
#'
#'   \code{"predSurface"} are the predicted effects (without replicates)
#'   according to the specified null model. These effects are arranged in a
#'   matrix form so that each direction of the matrix rightward or downward
#'   corresponds to increasing dose of one of the compounds.
#'
#'   \code{"occupancy"} contains occupancy levels at each dose
#'   combination as (always) computed by generalized Loewe model.
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   transforms <- getTransformations(data)
#'   fitResult <- fitMarginals(data, transforms)
#'   predictOffAxis(data, fitResult, null_model = "hsa")
predictOffAxis <- function(data, fitResult,
                           transforms = fitResult$transforms,
                           null_model = c("loewe", "hsa", "bliss", "loewe2"), ...) {

  ## Argument matching
  null_model <- match.arg(null_model)

  uniqueDoses <- with(data, list("d1" = sort(unique(data$d1)),
                                 "d2" = sort(unique(data$d2))))
  doseGrid <- expand.grid(uniqueDoses)

  occupancy <- NULL
  
  predSurface <- array(dim = lengths(uniqueDoses))
  
  predSurface[] <- switch(null_model,
                          "loewe" = {
                            fitLoewe <- generalizedLoewe(doseGrid, fitResult$coef)
                            occupancy <- fitLoewe$occupancy
                            fitLoewe$response
                          },
                          "hsa" = hsa(doseGrid, fitResult$coef),
                          "bliss" = Blissindependence(doseGrid, fitResult$coef),
                          "loewe2" = harbronLoewe(doseGrid, fitResult$coef)
                      )

  if (!is.null(transforms)) {
    CompositeT <- with(transforms,
                       function(y, args) PowerT(BiolT(y, args), args))
    data$effect <- with(transforms,
                        PowerT(data$effect, compositeArgs))
    predSurface <- with(transforms,
                        CompositeT(predSurface, compositeArgs))
  }

  dataOffAxis <- with(data, data[d1 & d2, , drop = FALSE])
  predOffAxis <- predSurface[cbind(match(dataOffAxis$d1, uniqueDoses$d1),
                                   match(dataOffAxis$d2, uniqueDoses$d2))]

  offaxisZTable <- cbind(dataOffAxis[, c("d1", "d2", "effect"), drop = FALSE],
                         "predicted" = predOffAxis)

  return(list("offaxisZTable" = offaxisZTable,
              "predSurface" = predSurface,
              "occupancy" = occupancy))
}
