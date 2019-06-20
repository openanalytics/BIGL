context("Response surface evaluation")

test_that('statistics 1', {

  expect_true(inherits(rs, "ResponseSurface"))
  expect_true(inherits(rs$meanR, "meanR"))
  expect_true(inherits(rs$maxR, "maxR"))

  expect_identical(rs$transforms, transforms)
  expect_identical(rs$meanR,
                   meanR(data, fit, transforms, "loewe", R, rs$CP, reps))

  
  ## Check equality without the attributes
  skip_on_cran()
  expect_equal(
      rs$maxR,
      maxR(data, fit, null_model = "loewe", B.CP = NULL,
          B.B = NULL, cutoff = 0.95, Ymean = offAxisTable,
          CP = rs$CP, reps = reps), 
      check.attributes = FALSE, 
      tolerance = 1e-4
  )

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

test_that("plots", {
      contour(rs)
      contour(rs, xlab = "x-label", ylab = "y-label")
      
      expect_silent(isobologram(rs))
      
      expect_silent({
            plot(rs, "maxR")
            rgl.close()
          })
    })
