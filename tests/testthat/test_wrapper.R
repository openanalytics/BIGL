context("Response surface evaluation")

test_that('statistics 1', {

  expect_true(inherits(rs, "ResponseSurface"))
  expect_true(inherits(rs$meanR, "meanR"))
  expect_true(inherits(rs$maxR, "maxR"))

  expect_identical(rs$transforms, transforms)
  expect_identical(rs$meanR,
		  			fitSurface(data, fit, transforms, "loewe", R, CP = rs$CP, reps, statistic = "meanR")$meanR)
  
  ## Check equality without the attributes - FIXME: not working
#  skip_on_cran()
#  expect_equal(
#      rs$maxR,
#      maxR(data, fit, null_model = "loewe", B.CP = NULL,
#          B.B = NULL, cutoff = 0.95, Ymean = offAxisTable,
#          CP = rs$CP, reps = reps), 
#      check.attributes = FALSE, 
#      tolerance = 1e-4
#  )

})

test_that('statistics 2', {

  ## Check for correct classes and outputs if both are bootstrapped
  rs <- fitSurface(data, fit, statistic = "both",
                   B.CP = 2, B.B = 2, parallel = FALSE)

  expect_true(inherits(rs, "ResponseSurface"))
  expect_true(inherits(rs$meanR, "meanR"))
  expect_true(inherits(rs$maxR, "maxR"))
  expect_identical(rs$transforms, transforms)

  expect_silent(summary(rs))

})

test_that('root finders', {
  expect_silent(fitSurface(data, fit, transforms = transforms,
                           B.CP = 2, B.B = NULL, parallel = FALSE,
                           statistic = "both", newtonRaphson = FALSE))
  expect_silent(fitSurface(data, fit, transforms = transforms,
                           B.CP = 2, B.B = NULL, parallel = FALSE,
                           statistic = "both", newtonRaphson = TRUE))
  expect_silent(fitSurface(data, fit, transforms = transforms,
                           null_model = "loewe2",
                           B.CP = 2, B.B = NULL, parallel = FALSE,
                           statistic = "both", newtonRaphson = FALSE))
  expect_silent(fitSurface(data, fit, transforms = transforms,
                           null_model = "loewe2",
                           B.CP = 2, B.B = NULL, parallel = FALSE,
                           statistic = "both", newtonRaphson = TRUE))
})

test_that("plots", {
#      contour(rs)
#      contour(rs, xlab = "x-label", ylab = "y-label")
      
      expect_silent(isobologram(rs))
      
      skip_on_cran()  # rgl fails on mac on CRAN
      
      expect_silent({
            plot(rs, "maxR")
            close3d()
          })

      expect_silent({
        plot(rs, "confInt")
        close3d()
      })

      expect_silent({
        plot(rs, "z-score")
        close3d()
      })
    })

test_that("Error and warnings", {
  expect_error(fitSurface(data, fit, transforms = transforms,
                          B.CP = NULL, B.B = NULL, parallel = FALSE,
                          statistic = "both"))
})
