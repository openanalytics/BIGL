#' Functions for scaling, and rescaling residuals
#' @details Residuals are calculated with respect to the average observation on
#' the off-axis point, so replicates are required!
#' @param sampling_errors A vector of raw  residuals
#' @param means a vector of means
#' @param model coefficients of a fitted mean-variance trend
scaleResids = function(sampling_errors, means, model){
    predVar = predictVar(means, model)
    sampling_errors/sqrt(predVar)
}
#' Backscale residuals
#' @param scaledResids scaled residuals
#' @inheritParams scaleResids
backscaleResids = function(scaledResids, means, model){
    predVar = predictVar(means, model)
    scaledResids*sqrt(predVar)
}
#'Predict variance
#' @inheritParams scaleResids
predictVar = function(means, model){
    predVar = model[1] + model[2]*means
    predVar[predVar<=0] = model[3] #Correct for negative variances
    predVar
}
#' Add residuals by adding to mean effects
#' @param method a variance method, see \code{\link{fitSurface}}
#' @inheritParams scaleResids
addResids = function(means, sampling_errors, method, model){
    means + sampleResids(means, sampling_errors, method, model)
}
#' Sample residuals according to a new model
#' @inheritParams addResids
#' @return sampled residuals
sampleResids = function(means, sampling_errors, method, model){
    if(method %in% c("equal", "unequal")){
        return(sample(sampling_errors, replace = TRUE))
    } else if(method == "model"){
        scaledResids = scaleResids(sampling_errors, means, model)
        sampledResids = sample(scaledResids, replace = TRUE)
        return(backscaleResids(sampledResids, means, model))
    } else{}
}