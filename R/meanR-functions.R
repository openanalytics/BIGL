#' Summary of meanR object
#'
#' @param object Output from \code{\link{meanR}}
#' @param ... Further arguments
#' @export
summary.meanR <- function(object, ...) {

  ans <- list()
  ans$bootstrapped <- !is.null(object$FDist)
  ans$dfs <- c("n1" = object$n1, "df0" = object$df0)
  ans$res <- data.frame("F" = object$FStat,
                        "p-value" = object$p.value,
                        check.names = FALSE)

  class(ans) <- append("summary.meanR", class(ans))
  ans
}

#' Print summary of meanR object
#'
#' @param x Summary of meanR object
#' @inheritParams summary.meanR
#' @export
print.summary.meanR <- function(x, ...) {

  type <- if (x$bootstrapped) "Bootstrapped" else "Exact"

  if (x$res$`p-value` < 2e-16)
    x$res$`p-value` <- "< 2e-16"
  else
    x$res$`p-value` <- paste("=", round(x$res$`p-value`, 4))

  cat(paste(type, "meanR test (H0 = no synergy/antagonism):\n"))
  cat("\tF(", x$dfs["n1"], ",", x$dfs["df0"], ") = ",
      round(x$res$F, 4), " (p-value ",
      x$res$`p-value`, ")\n", sep="")

}

#' Plot bootstrapped cumulative distribution function of meanR null distribution
#'
#' @param x Output from \code{\link{meanR}}
#' @inheritParams summary.meanR
#' @export
plot.meanR <- function(x, ...) {

  if (!is.null(x$FDist))
    plot(x$FDist)
  else
    stop("CDF plotting only available for bootstrapped meanR.")

}
