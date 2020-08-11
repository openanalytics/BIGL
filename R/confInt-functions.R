#' Summary of confidence intervals object
#'
#' @param object Output from \code{\link{bootConfInt}}
#' @param ... Further arguments
#' @export
summary.BIGLconfInt <- function(object, ...) {

ans <- list()
ans$estimate = object$single$meanEffect
ans$sigLevel = paste0(round(object$cutoff*100), "%")
ans$singleCI = object$single$confIntMeanEffect
ans$call = object$single$Call

ans$confInt = object$offAxis[object$offAxis$call %in% c("Syn", "Ant"),]
ans$confInt[, c("estimate", "lower", "upper")] = round(ans$confInt[, c("estimate", "lower", "upper")], 4)
ans$totals <- data.frame("Syn" = sum(object$offAxis$call == "Syn"),
                         "Ant" = sum(object$offAxis$call == "Ant"),
                         "Total" = nrow(object$offAxis))
rownames(ans$totals) = NULL

class(ans) <- append("summary.BIGLconfInt", class(ans))
ans
}

#' Print summary of BIGLconfInt object
#'
#' @param x Summary of BIGLconfInt object
#' @inheritParams summary.BIGLconfInt
#' @export
print.summary.BIGLconfInt <- function(x, ...) {

    #Overall
    cat("Overall effect\n")
    cat(sep = "", "Estimated mean departure from null response surface with ",
        x$sigLevel, " confidence interval: ", round(x$estimate, 4), " [", round(x$singleCI[1], 4), ", ", round(x$singleCI[2], 4), "]\n")
    cat("Evidence for effects in data:", x$call, "\n\n")

    #Pointwise
    cat("Pointwise effects\n")
    print(x$confInt)
    cat("\nPointwise", x$sigLevel, "confidence intervals summary:\n")
    print(x$totals)
    cat("\n")
}

#' Plot confidence intervals in a contour plot
#'
#' @param x Output from \code{\link{BIGLconfInt}}
#' @inheritParams summary.BIGLconfInt
#' @export
plot.BIGLconfInt <- function(x, ...) {

    if (!is.null(x$FDist))
        plot(x$FDist)
    else
        stop("CDF plotting only available for bootstrapped BIGLconfInt.")

}