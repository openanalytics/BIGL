#' Plot response surface
#'
#' Plot the 3-dimensional response surface predicted by one of the null
#' models. This plot allows for a visual comparison between the null
#' model prediction and observed points. This function is mainly used
#' as the workhorse of \code{\link{plot.ResponseSurface}} method.
#'
#' Title for the plot and legend are drawn as bitmaps and do not rotate with the
#' rest of the plot. Since they are bitmaps, they do not scale properly, hence
#' resizing window will result in unappealing visuals. For them to look
#' properly, it suffices to set the appropriate RGL window size and rerun the
#' plotting command.
#'
#' @param predSurface Vector of all predicted responses based on
#'   \code{expand.grid(uniqueDoses)}. If not supplied, it will be computed
#'   with \code{\link{predictOffAxis}} function.
#' @param null_model If \code{predSurface} is not supplied, it is computed using
#'   one of the available null models, i.e. \code{"loewe"} or \code{"hsa"}. See
#'   also \code{\link{fitSurface}}.
#' @param breaks Numeric vector with numerical breaks. To be used in conjunction
#'   with \code{colorPalette} argument.
#' @param colorPalette Vector of color names for surface
#' @param colorBy This parameter determines values on which coloring is based
#'   for the 3-dimensional surface. If matrix or a data frame with \code{d1} and
#'   {d2} columns is supplied, dose combinations from \code{colorBy} will be
#'   matched automatically to the appropriate dose combinations in \code{data}.
#'   Unmatched dose combinations will be set to 0. This is especially useful for
#'   plotting results for off-axis estimates only, e.g. off-axis Z-scores or
#'   maxR test statistics. If \code{colorBy = "colors"}, surface will be colored
#'   using colors in \code{colorPalette} argument.
#' @param colorPoints Colors for off-axis and on-axis points. Character vector
#'   of length four with colors for 1) off-axis points; 2) on-axis points of the
#'   first drug (i.e. second drug is dosed at zero); 3) on-axis points of the
#'   second drug; 4) on-axis points where both drugs are dosed at zero.
#' @param radius Radius of spheres. If missing, an educated guess based on
#'   number of digits in average effect will be made.
#' @param logScale Draw doses on log-scale (setting zeroes to be finite constant)
#' @param zTransform Optional transformation function for z-axis. By default,
#'   identity function is used.
#' @param main Fixed non-moving title for the 3D plot
#' @param legend Whether legend should be added
#' @param add Add the predicted response surface to an existing plot. Will not
#'   draw any points, just the surface. Must be called after another call to
#'   \code{\link{plotResponseSurface}}.
#' @param xat x-axis ticks: "pretty", "actual" or a numeric vector
#' @param yat y-axis ticks: "pretty", "actual" or a numeric vector
#' @param colorfun If replicates in \code{colorBy} variable are present, these
#'   will be aggregated using \code{colorfun} function. This can also be a
#'   custom function returning a scalar.
#' @param plotfun If replicates for dose combinations in \code{data} are
#'   available, points can be aggregated using \code{plotfun} function.
#'   Typically, it will be \code{\link{mean}}, \code{\link[stats]{median}},
#'   \code{\link{min}} or \code{\link{max}} but a custom-defined function
#'   returning a scalar from a vector is also possible.
#' @param ... Further arguments to format axis labels
#' @return Plot is shown on a rgl device.
#' @import rgl
#' @importFrom stats median predict quantile
#' @importFrom grDevices axisTicks colorRampPalette colors terrain.colors
#' @importFrom graphics plot.new
#' @inheritParams fitSurface
#' @export
#' @examples
#' \dontrun{
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   data_mean <- aggregate(effect ~ d1 + d2, data = data[, c("d1", "d2", "effect")],
#'                          FUN = mean)
#'
#'   ## Construct the surface from marginal fit estimates based on HSA
#'   ## model and color it by mean effect level
#'   plotResponseSurface(data, fitResult, null_model = "hsa",
#'                       colorBy = data_mean, breaks = 10^(c(0, 3, 4, 6)),
#'                       colorPalette = c("grey", "blue", "green"))
#'
#'   ## Response surface based on Loewe additivity model and colored with
#'   ## rainbow colors. Legend will not be displayed in any case.
#'   plotResponseSurface(data, fitResult, null_model = "loewe",
#'                       colorBy = "colors", colorPalette = rainbow(6))
#' }
plotResponseSurface <- function(data, fitResult = NULL,
                                transforms = fitResult$transforms,
                                predSurface = NULL, null_model = c("loewe", "hsa", "bliss"),
                                colorPalette = c("blue", "grey70", "red"),
                                colorBy = "none",
                                colorPoints = c("black", "sandybrown", "brown", "white"),
                                breaks = c(-Inf, 0, Inf),
                                radius = NULL,
                                logScale = TRUE,
                                colorfun = median,
                                zTransform = function(x) x,
                                add = FALSE,
                                main = "", legend = TRUE,
                                xat = "actual", yat = "actual",
                                plotfun = NULL, ...) {

  ## Argument matching
  null_model <- match.arg(null_model)

  if (missing(fitResult) & missing(predSurface))
    stop("Marginals fit result or predicted surface need to be supplied.")

  if (is.character(colorBy) & all(colorBy %in% colors())) {
    colorPalette <- colorBy
    colorBy <- "colors"
  }

  ## Calculate extra arguments
  uniqueDoses <- with(data, list("d1" = sort(unique(d1)),
                                 "d2" = sort(unique(d2))))
  doseGrid <- expand.grid(uniqueDoses)
  logT <- function(z) log(z + 0.5 * min(z[z != 0]))
  log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
  ## Transform function for doses
  transformF <- if (logScale) log10T else function(z) z

  zGrid <- predSurface
  ## If marginal fit information is provided, response surface can be
  ## automatically calculated.
  if (is.null(predSurface)) {
    respSurface <- predictOffAxis(data, fitResult,
                                  null_model = null_model,
                                  transforms = transforms)
    if (!is.null(transforms)) {
      predSurface <- with(transforms,
                          InvPowerT(respSurface$predSurface, compositeArgs))
    } else {
      predSurface <- respSurface$predSurface
    }
    zGrid <- predSurface
  }

  ## Rougly, radius is proportional to the number of digits of mean response
  if (missing(radius)) {
    avgEffect <- abs(mean(zTransform(zGrid)))
    radius <- max(10^(log10(avgEffect) - 1.5), 0.05)
  }

  ## If colorVec is a matrix with d1 and d2 columns, reorder it so that it
  ## matches the data ordering.
  if (inherits(colorBy, c("matrix", "data.frame"))) {
    stopifnot(c("d1", "d2") %in% colnames(colorBy))
    colorVec <- colorBy
    colorBy <- "asis"

    coloredBy <- colorVec
    cols <- colnames(coloredBy)
    pCols <- which(!(cols %in% c("d1", "d2")))
    if (length(pCols) > 1) pCols <- min[pCols]

    if (any(duplicated(coloredBy[, c("d1", "d2")]))) {
      coloredBy <- aggregate(coloredBy[, pCols],
                             by = coloredBy[, c("d1", "d2")], FUN = colorfun)
    }

    cols <- colnames(coloredBy)
    pCols <- which(!(cols %in% c("d1", "d2")))
    if (length(pCols) > 1) pCols <- min[pCols]

    colorVec <- rep(NA, nrow(doseGrid))
    for (i in 1L:nrow(doseGrid)) {
      ind <- match(paste(doseGrid[i,], collapse=";"),
                   apply(coloredBy[, c("d1", "d2")], 1,
                         function(x) paste(x, collapse=";")))
      if (!is.na(ind)) colorVec[i] <- coloredBy[[pCols]][ind]
    }
  }

  dataOffAxis <- with(data, data[d1 & d2, , drop = FALSE])
  predOffAxis <- predSurface[cbind(match(dataOffAxis$d1, uniqueDoses$d1),
                                   match(dataOffAxis$d2, uniqueDoses$d2))]
  if (nrow(dataOffAxis) == 0) {
    warning("No off-axis observations were found. Surface won't be custom colored..")
    colorBy <- "none"
  }

  ## This function computes a grid of colors for the 3d surface. Values
  ## that exceed specified boundary in absolute value take the boundary
  ## color.
  surfaceColors <- colorRampPalette(colorPalette)(length(breaks) - 1)

  surfaceColor <- function(response) {
    ff <- cut(response, breaks = breaks, include.lowest = TRUE)
    zcol <- surfaceColors[ff]
    return(zcol)
  }

  getLabels <- function(response) {
    ff <- cut(response, breaks = breaks, include.lowest = TRUE)
    labels <- gsub(",", ", ", levels(ff))
    return(labels)
  }

  ## Generate colors for the surface plot
  if (colorBy == "asis") {
    colorVec[is.na(colorVec)] <- 0
    zcol <- surfaceColor(colorVec)
    labels <- getLabels(colorVec)

  } else if (colorBy == "colors") {
    ## Use specified colors and recycle if necessary
    zcol <- rep(colorPalette, length(zGrid))

  } else {
    ## Generate colors from terrain.colors by looking at range of zGrid
    zGridFloor <- floor(100 * zGrid)
    col <- terrain.colors(diff(range(zGridFloor, na.rm = TRUE)))
    zcol <- col[zGridFloor - min(zGridFloor, na.rm = TRUE) + 1]
  }

  ##
  ## 3-dimensional surface plotting
  ##
  labnames <- c("Response", "Compound 1", "Compound 2")
  if (!is.null(attr(data, "orig.colnames")))
    labnames <- attr(data, "orig.colnames")

  if (!add) {
    ## Plot with no axes/labels
    if (!is.null(plotfun))
      data <- aggregate(effect ~ d1 + d2, data, FUN = plotfun)[, names(data)]

    plot3d(transformF(data$d1), transformF(data$d2),
           zTransform(data$effect), xlab = "", ylab = "", zlab = "",
           box = FALSE, axes = FALSE)

    ## Get x and y ticks
    ## TODO: use scales:log_breaks()(range(x))
    if (!is.numeric(xat)) {
      xat <- match.arg(xat, c("pretty", "actual"))
      if (xat == "pretty") {
        xlab <-  axisTicks(range(transformF(uniqueDoses$d1)),
                           log = logScale, nint = 3)
        if (logScale && length(xlab) > 4)
          xlab <- xlab[!(log10(xlab) %% 1)]
      } else xlab <- uniqueDoses$d1
      xat <- transformF(xlab)
    } else {
      xlab <- xat
      xat <- transformF(xat)
    }
    if (!is.numeric(yat)) {
      yat <- match.arg(yat, c("pretty", "actual"))
      if (yat == "pretty") {
        ylab <-  axisTicks(range(transformF(uniqueDoses$d2)),
                           log = logScale, nint = 3)
        if (logScale && length(ylab) > 4)
          ylab <- ylab[!(log10(ylab) %% 1)]
      } else ylab <- uniqueDoses$d2
      yat <- transformF(ylab)
    } else {
      ylab <- yat
      yat <- transformF(yat)
    }
    ## add points
    with(data, spheres3d(transformF(d1), transformF(d2),
                         zTransform(effect), radius = radius, add = TRUE,
                         col = colorPoints[1 + 1 * (d2 == 0) + 2 * (d1 == 0)]))
  }

  planes3d(0, 0, 1, zTransform(0),
           col = "grey", lit = FALSE, alpha = 0.3)
  persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2),
          zTransform(zGrid), add = TRUE, col = zcol, alpha = 0.6,
          lit = FALSE, aspect = FALSE, lwd = 2)

  bgplot3d({
    plot.new()
    title(main = main)

    if (legend && is.character(labels))
      legend("topleft", legend = labels, pch = 16, cex = 1.5, col = surfaceColors)
  })


  if (!add){
    ## Add box and annotation
    ## Must come after persp3d, otherwise only get wireframe
    rgl.bbox(xlen = 0, ylen = 0, zlen = 0, # no tickmarks
             expand = 1.03,
             color = "#000000", front = "lines", back = "cull")
    axis3d(edge = "x--", at = xat, labels = format(xlab, ...))
    axis3d(edge = "y--", at = yat, labels = format(ylab, ...))
    axis3d(edge = "z--")
    mtext3d(labnames[2], edge = "x--", line = 2)
    mtext3d(labnames[3], edge = "y--", line = 2)
    mtext3d(labnames[1], edge = "z--", line = 2)
  }

  persp3d(transformF(uniqueDoses$d1), transformF(uniqueDoses$d2),
          zTransform(zGrid), add = TRUE, col = zcol, alpha = 0.6,
          lit = FALSE, aspect = FALSE, lwd = 2)

}
