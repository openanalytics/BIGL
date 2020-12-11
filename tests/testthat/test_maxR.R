context("maxR")

test_that('statistics', {

  expect_silent(MaxR <- rs$maxR)

  expect_true(inherits(MaxR, "maxR"))
  expect_true(exists("q", where = attributes(MaxR$Ymean)) &
              is.numeric(attr(MaxR$Ymean, "q")))

  expect_true(exists("df0", where = attributes(MaxR$Ymean)) &
              is.numeric(attr(MaxR$Ymean, "df0")))

  expect_true(exists("cutoff", where = attributes(MaxR$Ymean)) &
              is.numeric(attr(MaxR$Ymean, "cutoff")))

  expect_true(exists("distr", where = attributes(MaxR$Ymean)) &
              inherits(attr(MaxR$Ymean, "distr"), "ecdf"))

  expect_silent(summary(MaxR))
      
})
