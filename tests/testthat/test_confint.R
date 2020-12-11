context("confidence intervals")

test_that('confidence intervals', {

    ## Try computation and check that all values are coherent
    expect_silent(confInt <- rs$confInt)
    #Check plot
    expect_silent(plotConfInt(rs))

    # Off axis confidence intervals
    expect_true(all(with(confInt, offAxis[,"lower"] < offAxis[,"upper"])))
    expect_true(all(with(confInt, offAxis[,"estimate"] < offAxis[,"upper"])))
    expect_true(all(with(confInt, offAxis[,"lower"] < offAxis[,"estimate"])))

    # Single measure confidence intervals
    expect_true(confInt$single$confIntMeanEffect[1] <=
                        confInt$single$confIntMeanEffect[2])
    #Only confidence interval with bootstrap
    CP = diag(with(data, sum(d1&d2)))
    colnames(CP) = rownames(CP) = getd1d2(data)[with(data, d1&d2)]
    expect_warning(rsci <- fitSurface(data, fit, transforms = transforms,
                                           B.CP = NULL, B.B = NULL, parallel = FALSE,
                                           statistic = "both", CP = CP))
    expect_null(rsci$confInt)
})
