#' Make a mean-variance plot
#' @param data a dataset or matrix with d1, d2 and effect column
#' @param trans,invtrans the transformation function for the variance and its inverse, possibly as strings
#' @param ... passed on to plot()
#' @return Plots the mean-variance trend
#' @export
#' @importFrom graphics lines
#' @details This is a crucial graphical check for deciding on the
plotMeanVarFit = function(data, trans = "identity",
                          invtrans = switch(trans, "identity" = "identity", "log" = "exp"),...){
    transFun = match.fun(trans); invtransFun = match.fun(invtrans)
    if(!all(c("d1", "d2", "effect") %in% colnames(data)))
        stop("Data must contain d1, d2 and effect columns!")
    dat_off = data[data$d1 & data$d2,]
    dat_off$d1d2 = getd1d2(dat_off)
    off_mean <- with(dat_off, tapply(effect, d1d2, mean))
    off_var = with(dat_off, tapply(effect, d1d2, var))
    predVar = modelVar(dat_off, transFun, invtransFun)
    if(any(predVar<0))
        warning("Negative variances modelled!\n")
    plot(off_mean, off_var, log = switch(trans, "identity" = "", "log" = "y", ""),
         ylab = "Variance", xlab ="Mean",...)
    lines(off_mean, predVar)
}