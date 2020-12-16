#' Calculate model variance, assuming variance increases linearly with mean
#' @importFrom stats lm.fit
#' @param dat_off off-axis points data
#' @param transFun,invTransFun the transformation and inverse transformation functions for the variance
#' @return the predicted model variance
modelVar = function(dat_off, transFun, invTransFun){
    off_mean <- with(dat_off, tapply(effect, d1d2, mean))
    off_var = with(dat_off, tapply(effect, d1d2, var))
    linmod <- c(lm.fit(transFun(off_var), x = cbind(1, off_mean))$coef,
                "min" = min(off_var), "max" = max(off_var))
    Var = predictVar(off_mean, linmod, invTransFun)
    Var
}