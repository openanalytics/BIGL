
test_that('synergy plot by component', {

    skip_on_cran()

    # Check plot for every null model
    null_model <- c("loewe", "loewe2", "bliss", "hsa")
    expect_no_error({
      rs_list <- lapply(null_model, function(x_model) {
        fitSurface(data, fit, null_model = x_model, B.CP = 50, statistic = "none", parallel = FALSE)
      })
      names(rs_list) <- null_model
    })

    # Check plots objects inherit the right classes
    expect_warning(
      p_bycomp <- synergy_plot_bycomp(
        ls = rs_list, xlab = NULL, plotBy = "Drug A", color = TRUE, ylab = "Response"
      ), "Unrecognized name"
    )

    expect_silent(
        p_bycomp <- synergy_plot_bycomp(
            ls = rs_list, xlab = NULL, plotBy = "Compound 2", color = TRUE, ylab = "Response"
        )
    )
    expect_equal(sort(class(p_bycomp)), sort(c("gg", "ggplot")))   
})


