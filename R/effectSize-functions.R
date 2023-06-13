#' Plot of effect-size object
#'
#' @param x Object of class \code{effect-size}.
#' @param colorPalette Vector of color values
#' @param logScale logScale
#' @param zTransform zTransform
#' @param digitsFunc Function to be applied to numeric values like doses. This expects a single parameter.
#' @param digits Numeric value indicating the number of digits used for numeric values. Whether \code{digitsFunc} is provided, this will be ignored.
#' @param ... Further arguments that are passed to \code{\link{format}} function
#'   for formatting of axis labels
#' @inheritParams graphics::title
#' @importFrom graphics axis filled.contour points title
#' @importFrom grDevices extendrange rgb
#' @importFrom ggplot2 ggplot scale_colour_stepsn
#' @importFrom data.table rbindlist
#' @export
`plot.effect-size` <- function(
  x,
  main = "Contour plot for effect size",
  xlab = "Dose (Compound 1)", ylab = "Dose (Compound 2)",
  colorPalette,
  logScale = TRUE, 
  zTransform = function(z) { z },
  digits,
  digitsFunc,
  ...
) {
  
  labels <- names(colorPalette)
  if (is.null(labels)) stop("Names for the vector `colorPalette` are mandatory")
  
  if (missing(digitsFunc)) {
    if (!missing(digits)) {
      digitsFunc <- function(x) round(x, digits = digits)
    } else {
      digitsFunc <- function(x) { x }
    }
  }
  
  if ("maxR" %in% names(x)) {
    synOut <- x$maxR$Ymean
    names(synOut)[names(synOut) == "call"] <- "synCall"
    
    effectOut <- x$confInt$offAxis
    names(effectOut)[names(effectOut) == "call"] <- "effectCall"
    effectOut$d1 <- as.numeric(gsub("(.+)_.+", "\\1", rownames(effectOut)))
    effectOut$d2 <- as.numeric(gsub(".+_(.+)", "\\1", rownames(effectOut)))
    x <- merge(synOut, effectOut, by = c("d1","d2"))
  } else {
    x <- x$offAxis
    names(x)[names(x) == "call"] <- "effectCall"
    #show doses on equidistant grid
    d1d2 <- rownames(x)
    d1d2split <- sapply(d1d2, function(y) strsplit(y, split = "_")[[1]])
    x$d1 <- as.numeric(d1d2split[1,])
    x$d2 <- as.numeric(d1d2split[2,])
  }
  uniqueDoses <- with(x, list("d1" = sort(unique(d1)),
                              "d2" = sort(unique(d2))))
  doseGrid <- expand.grid(uniqueDoses)
  
  log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
  transformF <- if (logScale) log10T else function(z) z

  breaks <- 1:3
  colourVec  <- colorPalette #colorRampPalette(colorPalette)(length(breaks) - 1)
  
  breaksInfo <- data.frame(
    breaks  = breaks,
    label   = labels,
    colour  = colorPalette,
    effect  = labels
  )  
  # This will be used to plot a continuous & invisible variable to get a continuous color bar legend
  breaksInfo$x <- seq(min(x$d1), max(x$d1), length.out = nrow(breaksInfo))
  breaksInfo$y <- seq(min(x$d2), max(x$d2), length.out = nrow(breaksInfo))
  
  if (nrow(x) < length(unique(x$d1))*length(unique(x$d2))) {
    # Then we need to add the missing combinations in order to allow geom_contour_filled to generate a contour
    x2 <- expand.grid(d1 = x$d1, d2 = x$d2)
    x2$effectCallNum <- which(labels == "None")
    x2$effectCall <- "None"
    x2$estimate <- 1e-8
    x <- rbindlist(list(x, x2), fill = TRUE)
    x <- x[!duplicated(x[, c("d1", "d2")]), ]
  }
  
  adjFactor <- 10
  x$effectCallNum <- as.numeric(factor(x$effectCall, levels = labels))
  # Fix glitch of blue line around the plot when all values are "None"
  if (all(x$effectCallNum == which(labels == "None"))) x$effectCallNum[1] <- x$effectCallNum[1] + 0.0001
  x$effectCallNum <- x$effectCallNum/adjFactor - 1/(adjFactor*2)
  x$d1_t <- transformF(x$d1)
  x$d2_t <- transformF(x$d2)
  
  # When only one point is Ant/Syn, geom_contour_filled won't be able to display any coloured polygon, so we are 
  # adding an artificial value next to it
  x <- rbind(x, x[x$effectCall %in% c("Syn", "Ant"),])
  p <- ggplot(data = x, aes(x = .data$d1_t, y = .data$d2_t)) + 
    geom_contour_filled(
      aes(z = .data$effectCallNum, colour = .data$effectCallNum), 
      breaks = c(0, breaks/adjFactor)
    ) + 
    geom_point(
      colour = rgb(0, 0, 0, 0.3),
      size = abs(x$estimate)/max(abs(x$estimate))*4 ## Size proportional to the effect size (normalized to be from 0 to 1)
    ) +
    scale_fill_manual("Call:", values = as.character(colourVec), labels = labels, drop = FALSE) + # values = colorPalette
    scale_x_continuous(
      breaks = unique(x$d1_t), 
      labels = digitsFunc(unique(x$d1)), 
    ) + 
    scale_y_continuous(
      breaks = unique(x$d2_t), 
      labels = digitsFunc(unique(x$d2))
    ) + 
    theme(
      panel.background = element_rect(
        fill = "white"
      ),
      axis.ticks.length = unit(0.1, "cm"),
      axis.ticks   = element_line(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
      legend.position = "bottom"
    ) +
    guides(colour = "none") +
    labs(title = main) + xlab(xlab) + ylab(ylab) 
  
  p
  
}
