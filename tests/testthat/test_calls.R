context("maxR calls for synergy/antagonism")

test_that("calls", {

  ## Decreasing marginal curves
  sims <- genData(parsDec, 2, -1)
  rs <- fitSurface(sims$data, sims$pars, statistic = "maxR",
                   B.CP = 2, B.B = NULL, parallel = FALSE)
  expect_equal(rs$maxR$Call, "Syn")
  expect_named(rs$occupancy, c("d1", "d2", "occupancy"))
  
  sims <- genData(parsDec, 2, 1)
  rs <- fitSurface(sims$data, sims$pars, statistic = "maxR",
                   B.CP = 2, B.B = NULL, parallel = FALSE)
  expect_equal(rs$maxR$Call, "Ant")

  ## Increasing marginal curves
  sims <- genData(parsInc, 2, 1)
  rs <- fitSurface(sims$data, sims$pars, statistic = "maxR",
                   B.CP = 2, B.B = NULL, parallel = FALSE)
  expect_equal(rs$maxR$Call, "Syn")

  sims <- genData(parsInc, 2, -1)
  rs <- fitSurface(sims$data, sims$pars, statistic = "maxR",
                   B.CP = 2, B.B = NULL, parallel = FALSE)
  expect_equal(rs$maxR$Call, "Ant")

  ## Diverging marginal curves
  expect_warning(sims <- genData(parsDiv, 2, 1), "Marginal curves are diverging")
  expect_warning(rs <- fitSurface(sims$data, sims$pars, statistic = "maxR",
                   B.CP = 2, B.B = NULL, parallel = FALSE), "Marginal curves are diverging")
  expect_equal(rs$maxR$Call, "Undefined")

  expect_warning(sims <- genData(parsDiv, 2, -1), "Marginal curves are diverging")
  expect_warning(rs <- fitSurface(sims$data, sims$pars, statistic = "maxR",
                   B.CP = 2, B.B = NULL, parallel = FALSE), "Marginal curves are diverging")
  expect_equal(rs$maxR$Call, "Undefined")

})

test_that("calls with no replicates", {
      
      ## Decreasing marginal curves
      sims <- genData(parsDec, 1, -1)
      rs <- fitSurface(sims$data, sims$pars, statistic = "both",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      
      expect_error(rs <- fitSurface(sims$data, sims$pars, statistic = "both",
              B.CP = 2, B.B = NULL, parallel = FALSE, method = "model"),
          "Replicates are required")
      
      expect_error(rs <- fitSurface(sims$data, sims$pars, statistic = "both",
              B.CP = 2, B.B = NULL, parallel = FALSE, method = "unequal"),
          "Replicates are required")
      
    })

test_that("HSA calls", {
      
      ## Decreasing marginal curves
      sims <- genData(parsDec, 2, -1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "hsa", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      expect_null(rs$occupancy) # no occupancy
      
      sims <- genData(parsDec, 2, 1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "hsa", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Ant")
      
      ## Increasing marginal curves
      sims <- genData(parsInc, 2, 1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "hsa", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      
      sims <- genData(parsInc, 2, -1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "hsa", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Ant")
            
    })

test_that("Bliss calls", {
      
      ## Decreasing marginal curves
      sims <- genData(parsDec, 2, -1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "bliss", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      expect_null(rs$occupancy) # no occupancy
      
      sims <- genData(parsDec, 2, 1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "bliss", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Ant")
      
      ## Increasing marginal curves
      sims <- genData(parsInc, 2, 1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "bliss", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      
      sims <- genData(parsInc, 2, -1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "bliss", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Ant")
      
      ## Diverging marginal curves
      expect_warning(sims <- genData(parsDiv, 2, 1), "Marginal curves are diverging")
      expect_error(
          rs <- fitSurface(sims$data, sims$pars, null_model = "bliss", statistic = "maxR",
              B.CP = 2, B.B = NULL, parallel = FALSE),
          "Bliss independence does not work for diverging marginal curves")
      
      expect_warning(sims <- genData(parsDiv, 2, -1), "Marginal curves are diverging")
      expect_error(
          rs <- fitSurface(sims$data, sims$pars, null_model = "bliss", statistic = "maxR",
              B.CP = 2, B.B = NULL, parallel = FALSE),
          "Bliss independence does not work for diverging marginal curves")
      
    })

test_that("Alternative Loewe calls", {
      
      ## Decreasing marginal curves
      sims <- genData(parsDec, 2, -1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "loewe2", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      expect_null(rs$occupancy) # no occupancy
      
      sims <- genData(parsDec, 2, 1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "loewe2", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Ant")
      
      ## Increasing marginal curves
      sims <- genData(parsInc, 2, 1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "loewe2", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Syn")
      
      sims <- genData(parsInc, 2, -1)
      rs <- fitSurface(sims$data, sims$pars, null_model = "loewe2", statistic = "maxR",
          B.CP = 2, B.B = NULL, parallel = FALSE)
      expect_equal(rs$maxR$Call, "Ant")
      
      ## Diverging marginal curves
      expect_warning(sims <- genData(parsDiv, 2, 1), "Marginal curves are diverging")
      expect_error(
          rs <- fitSurface(sims$data, sims$pars, null_model = "loewe2", statistic = "maxR",
              B.CP = 2, B.B = NULL, parallel = FALSE),
          "Alternative Loewe generalization does not work for diverging marginal curves.")
      
      expect_warning(sims <- genData(parsDiv, 2, -1), "Marginal curves are diverging")
      expect_error(
          rs <- fitSurface(sims$data, sims$pars, null_model = "loewe2", statistic = "maxR",
              B.CP = 2, B.B = NULL, parallel = FALSE),
          "Alternative Loewe generalization does not work for diverging marginal curves.")
      
    })
