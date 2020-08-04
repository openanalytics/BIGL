#' Obtain confidence intervals for the raw effect sizes on every off-axis point and overall
#'
#' @param Total data frame with all effects and mean effects
#' @inheritParams fitSurface
#' @inheritParams meanR
#' @param posEffect a boolean, are effects restricted to be positive
#' @param respS the observed response surface
#' @return A list with components
#' \item{offAxis}{The off-axis bootstrapped confidence intervals}
#' \item{single}{A mean effect and percentile and studentized boostrap intervals}
bootConfInt = function(Total, idUnique, bootStraps,
                       transforms, respS, B.B, method,
                       CP, reps, n1, cutoff, R, fitResult,
                       bootRS, data_off, posEffect = all(Total$effect >= 0), ...){
    Total = Total[Total$d1 & Total$d2,]
    sampling_errors <- Total$effect - Total$meaneffect
    A = getA(data_off, fitResult, method, CP, reps, n1)
    bootEffectSizesList = lapply(bootStraps, function(bb){
        #Do use bootstrapped response surface for complete mimicry of variability
        dat_off_resam = within(Total, {
            effect = meaneffect + sample(sampling_errors)
            if(posEffect){
                effect = abs(effect)
            }
        })
        #Transforms have occurred in Total already
        bootR = getR(data = dat_off_resam, idUnique = dat_off_resam$d1d2,
             transforms = NULL, respS = if(bootRS) bb$respS else respS)
        bootA = getA(dat_off_resam, bb$simFit, method, CP, reps, n1)
        list("R" = bootR, "A" = bootA)
    })
    bootEffectSizes = vapply(bootEffectSizesList, FUN.VALUE = c(R), function(x) x$R)
    bootAs = vapply(bootEffectSizesList, FUN.VALUE = diag(A), function(x) sqrt(diag(x$A)))
    #Off axis confidence interval
    bootEffectSizesStand = abs(bootEffectSizes-c(R))/bootAs
    maxEffectSizes = apply(bootEffectSizesStand, 2, max)
    effectSizeQuant = quantile(maxEffectSizes, cutoff)
    confInt = c(R) + outer(effectSizeQuant*sqrt(diag(A)),
                           c("lower" = -1, "upper" = 1))
    rownames(confInt) = rownames(bootEffectSizesStand)
    #Single measure of effect size
    singleMeasure = mean(R)
    bootR = colMeans(bootEffectSizes)
    bootRstand = (bootR-singleMeasure)/vapply(bootEffectSizesList, FUN.VALUE = double(1), function(x) mean(x$A))
    percentileCI = quantile(bootR, c("lower" = (1-cutoff)/2,
                                     "upper" = (1+cutoff)/2))
    sdA = mean(A)
    studentizedCI = singleMeasure + sdA*
                            quantile(bootRstand, c("lower" = (1-cutoff)/2,
                                                   "upper" = (1+cutoff)/2))
    return(list("offAxis" = confInt,
                "single" = list("meanEffect" = singleMeasure,
                                "confIntMeanEffect" =
                                    rbind("percentile" = percentileCI,
                                 "studentized" = studentizedCI))))
}