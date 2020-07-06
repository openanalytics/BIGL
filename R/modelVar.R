#' Calculate model variance, assuming variance increases linearly with mean
#' @importFrom stats lm.fit
#' @param dat_off off-axis points data
#' @return the predicted model variance
modelVar = function(dat_off){
    off_mean <- with(dat_off, tapply(effect, d1d2, mean))
    off_var = with(dat_off, tapply(effect, d1d2, var))
    linmod <- lm.fit(off_var, x = cbind(1, off_mean))
    linmod$coef[1] + linmod$coef[2]*off_mean
}