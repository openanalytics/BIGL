context("confidence intervals")

test_that('confidence intervals', {

    ## Try computation and check that all values are coherent
    expect_silent(confInt <- rs$confInt)

    # Off axis confidence intervals
    expect_type(confInt$offAxis, "double")
    expect_true(all(with(confInt, offAxis[,1] < offAxis[,2])))

    # Single measure confidence intervals
    expect_true(all(confInt$single$confIntMeanEffect[,1] <
                        confInt$single$confIntMeanEffect[, 2]))
    expect_true(all(confInt$single$confIntMeanEffect[,1] <
                        confInt$single$meanEffect))
    expect_true(all(confInt$single$confIntMeanEffect[,2] >
                    confInt$single$meanEffect))
})