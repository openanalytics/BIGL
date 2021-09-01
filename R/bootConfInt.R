#' Obtain confidence intervals for the raw effect sizes on every off-axis point and overall
#'
#' @param Total data frame with all effects and mean effects
#' @inheritParams fitSurface
#' @inheritParams meanR
#' @inheritParams generateData
#' @param posEffect a boolean, are effects restricted to be positive
#' @param respS the observed response surface
#' @return A list with components
#' \item{offAxis}{The off-axis bootstrapped confidence intervals}
#' \item{single}{A mean effect and percentile and studentized boostrap intervals}
bootConfInt = function(Total, idUnique, bootStraps,
    transforms, respS, B.B, method,
    CP, reps, n1, cutoff, R, fitResult,
    bootRS, data_off, posEffect = all(Total$effect >= 0),
    transFun, invTransFun, model, rescaleResids, ...) {
  Total <- Total[Total$d1 & Total$d2, ]
  sampling_errors <- Total$effect - Total$meaneffect
  A <- getA(data_off, fitResult, method, CP, reps, n1, transFun = transFun,
      invTransFun = invTransFun)
  bootEffectSizesList <- lapply(bootStraps, function(bb) {
        #Do use bootstrapped response surface for complete mimicry of variability
        dat_off_resam <- within(Total, {
              effect <- addResids(Total$meaneffect, sampling_errors, method,
                  rescaleResids, model, invTransFun)
              #Sample with replacement
              if (posEffect) {
                effect <- abs(effect)
              }
            })
        #Transforms have occurred in Total already
        bootR <- getR(data = dat_off_resam, idUnique = dat_off_resam$d1d2,
            transforms = NULL, respS = if(bootRS) bb$respS else respS)
        bootA <- getA(dat_off_resam, bb$simFit, method, CP, reps, n1, transFun = transFun,
            invTransFun = invTransFun)
        list("R" = bootR, "A" = bootA)
      })
  bootEffectSizes <- vapply(bootEffectSizesList, FUN.VALUE = c(R), function(x) x$R)
  bootAs <- vapply(bootEffectSizesList, FUN.VALUE = diag(A), function(x) sqrt(diag(x$A)))
  #Off axis confidence interval
  bootEffectSizesStand <- abs(bootEffectSizes-c(R))/bootAs
  maxEffectSizes <- apply(bootEffectSizesStand, 2, max)
  effectSizeQuant <- quantile(maxEffectSizes, cutoff, na.rm = TRUE)
  confInt <- c(R) + outer(effectSizeQuant*sqrt(diag(A)),
      c("lower" = -1, "upper" = 1))
  rownames(confInt) <- rownames(bootEffectSizesStand)
  
  coefFit <- fitResult$coef
  eq  <- coefFit["m1"] == coefFit["b"] && coefFit["m2"] == coefFit["b"]
  inc <- coefFit["m1"] >= coefFit["b"] && coefFit["m2"] >= coefFit["b"]
  dec <- coefFit["m1"] <= coefFit["b"] && coefFit["m2"] <= coefFit["b"]
  
  call <- rep("None", length(R))
  call[confInt[, "lower"] >= 0] <- if (eq) {
        "Undefined"
      } else if (inc) {
        "Syn"
      } else if (dec) {
        "Ant"
      } else 
        "Undefined"
  call[confInt[, "upper"] <= 0] <- if (eq) {
        "Undefined"
      } else if (inc) {
        "Ant"
      } else if (dec) {
        "Syn"
      } else 
        "Undefined"
  confInt <- data.frame("estimate" = R, confInt, "call" = call)
  
  #Single measure of effect size
  singleMeasure <- mean(R)
  bootR <- colMeans(bootEffectSizes)
  bootRstand <- (bootR-singleMeasure)/vapply(bootEffectSizesList, 
      FUN.VALUE = double(1), function(x) mean(x$A))
  sdA <- mean(A)
  studentizedCI <- singleMeasure + sdA*quantile(bootRstand, 
      c((1-cutoff)/2, (1+cutoff)/2), na.rm = TRUE)
  names(studentizedCI) <- c("lower", "upper")
  overallCall <- if (eq || any(is.na(studentizedCI))) {
        "Undefined"
      } else {
        if (studentizedCI["lower"] > 0) {
          if (inc) {
            "Syn"
          } else if (dec) {
            "Ant"
          } else
            "Undefined"
        } else if (studentizedCI["upper"] < 0) {
          if (inc) {
            "Ant"
          } else if (dec) {
            "Syn"
          } else 
            "Undefined"
        } else {
          "None"
        }
      }
  
  ans <- list("offAxis" = confInt,
      "single" = list("meanEffect" = singleMeasure,
          "confIntMeanEffect" = studentizedCI,
          "Call" = overallCall),
      "cutoff" = cutoff)
  class(ans) <- append("BIGLconfInt", class(ans))
  return(ans)
}