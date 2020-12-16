#' Helper functions for the test statistics
#' @param idUnique id of unique off axis points
#' @param data the datasets
#' @param respS the evaluated response surface
#'@inheritParams predictOffAxis
getR = function(data, idUnique, transforms, respS){
    if(!is.null(transforms)){
        data$effect <- with(transforms, PowerT(data$effect, compositeArgs))
    }
    tapply(data$effect - respS[idUnique], data$d1d2, mean)
}
#'@inheritParams predictOffAxis
#'@inheritParams getR
getMeanRF = function(data, fitResult, method, CP, reps, transforms, null_model,
                     R, n1, idUnique, respS, transFun, invTransFun){
    if(missing(R)){
        R = getR(data = data, idUnique = idUnique, transforms = transforms,
                 respS = respS)
    }
    A <- getA(data, fitResult, method, CP, reps, n1, transFun, invTransFun)
    FStat <- max(0, as.numeric(crossprod(R, solve(A)) %*% R / n1))
    return(FStat)
}
getMaxRF = function(data, fitResult, method, CP, reps, transforms, null_model,
                    R, n1, idUnique, respS, transFun, invTransFun){
    if(missing(R)){
        R = getR(data = data, idUnique = idUnique, transforms = transforms,
                 respS = respS)
    }
    A <- getA(data, fitResult, method, CP, reps, n1, transFun, invTransFun)
    E <- eigen(A)
    V <- E$values
    Q <- E$vectors
    Amsq <- Q %*% tcrossprod(diag(1/sqrt(V)), Q)
    RStud <- crossprod(R, Amsq)
    return(as.numeric(RStud))
}
getA = function(dat_off, fitResult, method, CP, reps, n1, transFun, invTransFun){
    MSE0 <- fitResult$sigma^2
    mse_off <- switch(method,
                      "equal" = MSE0, "model" = c(modelVar(dat_off, transFun, invTransFun)),
                      "unequal" = mean(with(dat_off, tapply(effect, d1d2, var)))
    )
    A <- MSE0*CP + mse_off*diag(1/reps, nrow = n1)
    return(A)
}