#' Helper functions for the test statistics
#'@inheritParams predictOffAxis
getR = function(data, fitResult, transforms, null_model){
    respS <- predictOffAxis(data = data, fitResult = fitResult,
                            transforms = transforms, null_model = null_model)
    R <- with(respS$offaxisZTable, tapply(effect - predicted, d1d2, mean))
}
#'@inheritParams predictOffAxis
getMeanRF = function(data, fitResult, method, CP, reps, transforms, null_model,
                     R, n1){
    if(missing(R)){
        R = getR(data = data, fitResult = fitResult, transforms = transforms,
                 null_model = null_model)
    }
    A <- getA(data, fitResult, method, CP, reps, transforms, null_model, R, n1)
    FStat <- max(0, as.numeric(t(R) %*% solve(A) %*% R / n1))
    return(FStat)
}
getMaxRF = function(data, fitResult, method, CP, reps, transforms, null_model, R, n1){
    A <- getA(data, fitResult, method, CP, reps, transforms, null_model, R, n1)
    E <- eigen(A)
    V <- E$values
    Q <- E$vectors
    Amsq <- Q %*% diag(1/sqrt(V)) %*% t(Q)
    RStud <- t(R) %*% Amsq
    return(as.numeric(RStud))
}
getA = function(data, fitResult, method, CP, reps, transforms, null_model, R, n1){
    MSE0 <- fitResult$sigma^2
    dat_off  <- data[data$d1 & data$d2, ]
    mse_off <- switch(method,
                      "equal" = MSE0,
                      "model" = modelVar(dat_off),
                      "unequal" = mean(with(dat_off, tapply(effect, d1d2, var)))
    )
    A <- MSE0*CP + mse_off*diag(1/reps, nrow = n1)
    return(A)
}