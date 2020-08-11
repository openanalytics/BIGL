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
#' @param x off axis confidence intervals, a data frame
#' @param ... additional arguments, currently ignored
#' @importFrom stats setNames
#' @export
#' @note after contour in the \code{drugCombo} package
plot.BIGLconfInt <- function(x, ...) {
        x = x$offAxis
        # prepare fill legend
        synCalls <- c("additive", "antagonism", "synergy")
        x$synLabel <- factor(x$call, labels = synCalls, levels = c("Add", "Ant", "Syn"))

        legendColors <- c("white", "pink", "lightblue")
        names(legendColors) = synCalls
        # subset to only the colors that are present in the data
        legendColors <- legendColors[names(legendColors) %in% as.character(unique(x$synLabel))]

        # text to show
        x$label <- sprintf("%.2f\n(%.2f, %.2f)", x$estimate, x$lower, x$upper)

        # show doses on equidistant grid
        d1d2 = rownames(x)
        d1d2split = sapply(d1d2, function(y) strsplit(y, split = "_")[[1]])
        x$d1 <- as.numeric(d1d2split[1,])
        x$d1 = factor(x$d1, levels = sort(unique(x$d1)),
                      labels = sort(unique(x$d1)), ordered = TRUE)
        x$d2 <- as.numeric(d1d2split[2,])
        x$d2 = factor(x$d2, levels = sort(unique(x$d2)),
                      labels = sort(unique(x$d2)), ordered = TRUE)

        p <- ggplot(data = x, aes_string(x = "d1", y = "d2")) +
            geom_tile(aes_string(fill = "synLabel"), color = "grey") +
            geom_text(aes_string(label = "label"), show.legend = FALSE, size = 3) +
            # invisible points, used only for labels
            geom_point(aes_string(color = "synLabel"), alpha = 0) +
            # round dose labels to digits
            scale_x_discrete(labels = format(as.numeric(levels(x$d1)), digits = 4)) +
            scale_y_discrete(labels = format(as.numeric(levels(x$d2)), digits = 4)) +
            scale_fill_manual(values = legendColors,
                              guide = FALSE) +
            scale_color_manual( # for a nicer legend
                values = setNames(seq_len(3), nm = synCalls),
                guide = guide_legend(title = "call:",
                                     override.aes = list(alpha = 1, shape = 22, size = 8, color = "grey",
                                                         fill = legendColors))
            ) +
            theme_minimal() +
            theme(
                panel.grid.major = element_blank(),
                legend.position = "bottom",
                axis.text.x = element_text(angle = 45, hjust = 1)
            )
        p
}
#' Plot confidence intervals from BIGL object in a contour plot
#'
#' @param BIGLobj Output from \code{\link{fitSurface}}
#' @param ... passed on to \code{\link{plot.BIGLconfInt}}
#' @export
plotConfInt <- function(BIGLobj, ...) {
    plot(BIGLobj$confInt, ...)
}