context("Plotting the mean-variance trend")

test_that('plotting the mean-variance trend', {

    ## Tsome basic checkcs for the plotting function
    expect_silent(plotMeanVarFit(data))
    expect_silent(plotMeanVarFit(data, trans = "log"))
})
