#' Functions for scaling, and rescaling residuals
#' @details Residuals are calculated with respect to the average observation on
#' the off-axis point, so replicates are required!
#' @param rawResids A vector of raw residuals
#' @param means a vector of means
#' @param model a variance model, see \code{\link{fitSurface}}
scaleResids = function(rawResids, means, model){
    predVar = predictVar(means, model)
    rawResids/sqrt(predVar)
}
backscaleResids = function(scaledResids, means, model){
    predVar = predictVar(means, model)
    scaleResids*sqrt(predVar)
}
predictVar(means, model){
    predVar = model[1] + model[2]*means
}