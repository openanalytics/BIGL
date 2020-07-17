#' Obtain confidence intervals for the raw effect sizes
bootConfInt = function(Total, idUnique, bootStraps,
                       transforms, respS, B.B, method,
                       CP, reps, n1, cutoff, R, fitResult,
                       bootRS, data_off,...){
    Total = Total[Total$d1 & Total$d2,]
    sampling_errors <- Total$effect - Total$meaneffect
    bootEffectSizes = vapply(bootStraps, FUN.VALUE = c(R), function(bb){
        #Do use bootstrapped response surface for complete mimicry of variability
        dat_off_resam = within(Total, {
            effect = meaneffect + sample(sampling_errors)
        })
        getR(data = dat_off_resam, idUnique = dat_off_resam$d1d2,
             transforms = NULL, respS = if(bootRS) bb$respS else respS)
    })
    A = getA(data_off, fitResult, method, CP, reps, n1)
    bootEffectSizesStand = abs((bootEffectSizes-c(R))/sqrt(diag(A)))
    maxEffectSizes = apply(bootEffectSizesStand, 2, max)
    effectSizeQuant = quantile(maxEffectSizes, cutoff)
    confInt = c(R) + outer(effectSizeQuant*sqrt(diag(A)),
                           c("lower" = -1, "upper" = 1))
    rownames(confInt) = rownames(bootEffectSizesStand)
    return(confInt)
}