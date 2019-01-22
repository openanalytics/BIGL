#' Isobologram of the response surface predicted by the null model
#'
#' If transformation functions are used, then the isobologram response levels
#' will be plotted on the transformed scale.
#'
#' @param x Output of \code{\link{fitSurface}}
#' @param grid.len Number of concentrations to plot for each compound in the
#'   contour plot. An evenly spaced grid of doses will be generated for each
#'   compound given its respective observed minimum and maximum doses. Note that
#'   \code{grid.len^2} computations will be needed later so this number should
#'   stay reasonably low.
#' @param logScale If \code{logScale = TRUE}, then grid of doses is evenly
#'   spaced in the logarithmic scale.
#' @param ... Further parameters that are not used at this moment.
#' @import ggplot2
#' @export
isobologram <- function(x, grid.len = 100, logScale = TRUE, ...) {

  ## Generate evenly spaced grid either on a linear or a log-scale
  genSeq <- function(doses) {

    if (logScale) {
      ## Log-scale removed zero dose
      doses <- setdiff(doses, 0)
      seq.range <- log(range(doses))
      c(0, exp(seq(seq.range[1], seq.range[2], length.out = grid.len - 1)))
    } else {
      ## Linear scale
      seq.range <- range(doses)
      seq(seq.range[1], seq.range[2], length.out = grid.len)
    }

  }

  coefs <- coef(x$fitResult)

  ## Generate a grid of doses for Compound 1 and predict the response
  doses1 <- genSeq(unique(x$data$d1))
  resp1 <- L4(doses1, coefs["h1"], coefs["b"], coefs["m1"], coefs["e1"])
  ## Generate a grid of doses for Compound 2 and predict the response
  doses2 <- genSeq(unique(x$data$d2))
  resp2 <- L4(doses2, coefs["h2"], coefs["b"], coefs["m2"], coefs["e2"])
  ## Combine both compounds and their marginal predictions
  data <- rbind(data.frame("d1" = doses1, "d2" = 0, "effect" = resp1),
                data.frame("d1" = 0, "d2" = doses2, "effect" = resp2))

  ## Based on marginal data, generate null model predictions
  predSurface <- predictOffAxis(data, x$fitResult,
                                null_model = x$null_model)$predSurface

  melt.surface <- data.frame("d1" = rep(doses1, length(doses2)),
                             "d2" = rep(doses2, each = length(doses1)),
                             "effect" = as.numeric(predSurface))

  labnames <- c("Response", 
      if (!is.null(x$names)) x$names else c("Compound 1", "Compound 2"))
  if (!is.null(attr(x$data, "orig.colnames"))) {
    labnames <- unlist(attr(x$data, "orig.colnames"))
  }

  p <- ggplot(melt.surface,
              aes_string(x = "d1", y = "d2",
                         z = "effect", fill = "effect")) +
    theme_bw() +
    geom_tile() +
    labs(x = labnames[2], y = labnames[3]) +
    scale_fill_gradientn(labnames[1],
                         colours = c("steelblue", "lightsteelblue", "lightblue",
                                     "floralwhite", "beige", "khaki",
                                     "orange1", "tomato3", "red")) +
    geom_contour(bins = 7, col = "black", size = 0.2)

  if (logScale) {
    p <- p +
      scale_x_log10(breaks = unique(x$data$d1)) +
      scale_y_log10(breaks = unique(x$data$d1))
  } else {
    p <- p +
      scale_x_continuous(breaks = unique(x$data$d1)) +
      scale_y_continuous(breaks = unique(x$data$d1))
  }

  p

}


