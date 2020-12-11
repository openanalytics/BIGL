#' Functions for scaling, and rescaling residuals. May lead to unstable behaviour in practice
#' @details Residuals are calculated with respect to the average observation on
#' the off-axis point, so replicates are required!
#' @param sampling_errors A vector of raw  residuals
#' @param ... passed on to predictVar
scaleResids = function(sampling_errors, ...){
    predVar = predictVar(...)
    sampling_errors/sqrt(predVar)
}
#' Backscale residuals
#' @param scaledResids scaled residuals
#' @inheritParams scaleResids
backscaleResids = function(scaledResids, ...){
    predVar = predictVar(...)
    scaledResids*sqrt(predVar)
}
#'Predict variance
#' @param means a vector of means
#' @inheritParams generateData
predictVar = function(means, model, invTransFun){
    predVar = invTransFun(model[1] + model[2]*means)
    predVar[predVar<=0] = model["min"] #Correct for negative variances
    predVar[predVar > model["max"]] = model["max"] #Upper bound
    predVar
}
#' Add residuals by adding to mean effects
#' @inheritParams scaleResids
#' @inheritParams predictVar
addResids = function(means, ...){
    means + sampleResids(means, ...)
}
#' Sample residuals according to a new model
#' @inheritParams fitSurface
#' @inheritParams predictVar
#' @inheritParams scaleResids
#' @return sampled residuals
sampleResids = function(means, sampling_errors, method, rescaleResids,...){
    if(method %in% c("equal", "unequal")){
        return(sample(sampling_errors, size = length(means), replace = TRUE))
    } else if(method == "model"){
        resids = if(rescaleResids){
            scaledResids = scaleResids(sampling_errors, means, ...)
            sampledResids = sample(scaledResids, replace = TRUE)
            backscaleResids(sampledResids, means, ...)
        } else{
            rnorm(length(means), sd = sqrt(predictVar(means, ...)))
        }
        return(resids)
    } else{}
}