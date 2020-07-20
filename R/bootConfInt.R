#' Obtain confidence intervals for the raw effect sizes on every off-axis point and overall
#'
#' @param Total data frame with all effects and mean effects
#' @inheritParams fitSurface
#' @inheritParams meanR
#' @return A list with components
#' \item{offAxis}{The off-axis bootstapped confidence intervals}
#' \item{single}{A mean effect and percentile and studentized boostrap intervals}
bootConfInt = function(Total, idUnique, bootStraps,
                       transforms, respS, B.B, method,
                       CP, reps, n1, cutoff, R, fitResult,
                       bootRS, data_off, posEffect = all(Total$effect >= 0), ...){
    Total = Total[Total$d1 & Total$d2,]
    sampling_errors <- Total$effect - Total$meaneffect
    bootEffectSizes = vapply(bootStraps, FUN.VALUE = c(R), function(bb){
        #Do use bootstrapped response surface for complete mimicry of variability
        dat_off_resam = within(Total, {
            effect = meaneffect + sample(sampling_errors)
            if(posEffect){
                effect = abs(effect)
            }
        })
        getR(data = dat_off_resam, idUnique = dat_off_resam$d1d2,
             transforms = transforms, respS = if(bootRS) bb$respS else respS)
    })
    A = getA(data_off, fitResult, method, CP, reps, n1)
    #Off axis confidence interval
    bootEffectSizesStand = abs((bootEffectSizes-c(R))/sqrt(diag(A)))
    maxEffectSizes = apply(bootEffectSizesStand, 2, max)
    effectSizeQuant = quantile(maxEffectSizes, cutoff)
    confInt = c(R) + outer(effectSizeQuant*sqrt(diag(A)),
                           c("lower" = -1, "upper" = 1))
    rownames(confInt) = rownames(bootEffectSizesStand)
    #Single measure of effect size
    singleMeasure = mean(R)
    bootR = colMeans(bootEffectSizes)
    percentileCI = quantile(bootR, c("lower" = (1-cutoff)/2,
                                     "upper" = (1+cutoff)/2))
    studentizedCI = singleMeasure + effectSizeQuant*sd(bootR)*
                                 c("lower" = -1, "upper" = 1)
    return(list("offAxis" = confInt,
                "single" = list("meanEffect" = singleMeasure,
                                "confIntMeanEffect" =
                                    rbind("percentile" = percentileCI,
                                 "studentized" = studentizedCI))))
}