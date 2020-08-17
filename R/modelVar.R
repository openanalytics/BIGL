#' Calculate model variance, assuming variance increases linearly with mean
#' @importFrom stats lm.fit
#' @param dat_off off-axis points data
#' @param transFun,invTransFun the transformation and inverse transformation functions for the variance
#' @return the predicted model variance
modelVar = function(dat_off, transFun, invTransFun){
    off_mean <- with(dat_off, tapply(effect, d1d2, mean))
    off_var = with(dat_off, tapply(effect, d1d2, var))
    linmod <- lm.fit(transFun(off_var), x = cbind(1, off_mean))
    Var = invTransFun(linmod$coef[1] + linmod$coef[2]*off_mean)
    Var[Var<=0] = min(off_var[off_var>0]) #Correct negative variances
    Var
}