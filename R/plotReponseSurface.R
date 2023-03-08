#' Plot response surface
#'
#' Plot the 3-dimensional response surface predicted by one of the null
#' models. This plot allows for a visual comparison between the null
#' model prediction and observed points. This function is mainly used
#' as the workhorse of \code{\link{plot.ResponseSurface}} method.
#'
#' @param data Object "data" from the output of \code{\link{fitSurface}}
#' @param fitResult Object "fitResult" from the output of \code{\link{fitSurface}}
#' @param transforms Object "transforms" from the output of \code{\link{fitSurface}}
#' @param predSurface Vector of all predicted responses based on
#'   \code{expand.grid(uniqueDoses)}. If not supplied, it will be computed
#'   with \code{\link{predictOffAxis}} function.
#' @param null_model If \code{predSurface} is not supplied, it is computed using
#'   one of the available null models, i.e. \code{"loewe"}, \code{"hsa"}, 
#'   \code{"bliss"} and \code{"loewe2"}. See also \code{\link{fitSurface}}.
#' @param breaks Numeric vector with numerical breaks. To be used in conjunction
#'   with \code{colorPalette} argument. If named, the labels will be displayed in the legend
#' @param colorPalette Vector of color names for surface
#' @param colorPaletteNA Color used in the matrix of colours when the combination of doses doesn't exist (NA)
#' @param colorBy This parameter determines values on which coloring is based
#'   for the 3-dimensional surface. If matrix or a data frame with \code{d1} and
#'   {d2} columns is supplied, dose combinations from \code{colorBy} will be
#'   matched automatically to the appropriate dose combinations in \code{data}.
#'   Unmatched dose combinations will be set to 0. This is especially useful for
#'   plotting results for off-axis estimates only, e.g. off-axis Z-scores or
#'   maxR test statistics. If \code{colorBy = "colors"}, surface will be colored
#'   using colors in \code{colorPalette} argument.
#' @param addPoints Boolean whether the dose points should be included
#' @param colorPoints Colors for off-axis and on-axis points. Character vector
#'   of length four with colors for 1) off-axis points; 2) on-axis points of the
#'   first drug (i.e. second drug is dosed at zero); 3) on-axis points of the
#'   second drug; 4) on-axis points where both drugs are dosed at zero.
#' @param radius Size of spheres (default is 4)
#' @param logScale Draw doses on log-scale (setting zeroes to be finite constant)
#' @param zTransform Optional transformation function for z-axis. By default,
#'   identity function is used.
#' @param main Fixed non-moving title for the 3D plot
#' @param legend Whether legend should be added (default FALSE)
#' @param add (deprecated) Add the predicted response surface to an existing plot. Will not
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
#' @param gradient Boolean indicating whether colours should be interpolated between breaks (default TRUE). 
#'   If FALSE, \code{colorPalette} must contain length(breaks)-1 colours
#' @param width Width in pixels (optional, defaults to 800px).
#' @param height Height in pixels (optional, defaults to 800px).
#' @param title String title (default "")
#' @param digitsFunc Function to be applied to the axis values
#' @param ... Further arguments to format axis labels
#' @return Plotly plot
#' @importFrom stats median predict quantile
#' @importFrom grDevices axisTicks colorRampPalette colors terrain.colors
#' @importFrom graphics plot.new
#' @importFrom plotly plot_ly add_surface add_markers layout config
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
#'   ## rainbow colors.
#'   plotResponseSurface(data, fitResult, null_model = "loewe", breaks = c(-Inf, 0, Inf),
#'                       colorBy = "colors", colorPalette = rainbow(6))
#' }
plotResponseSurface <- function(data, fitResult = NULL, 
                                transforms = fitResult$transforms, 
                                predSurface = NULL, null_model = c("loewe", "hsa", "bliss", "loewe2"),
                                colorPalette   = c("blue", "grey70", "red"), 
                                colorPaletteNA = "grey70",
                                colorBy = "none", 
                                addPoints = TRUE,
                                colorPoints = c("black", "sandybrown", "brown", "white"), 
                                breaks, # = c(-Inf, 0, Inf) 
                                radius = 4, 
                                logScale = TRUE, 
                                colorfun = median, 
                                zTransform = function(x) x, 
                                add = FALSE, 
                                main = "", 
                                legend = FALSE, 
                                xat = "actual", yat = "actual", 
                                plotfun = NULL,
                                gradient = TRUE,
                                width = 800, height = 800, 
                                title = "", 
                                digitsFunc = function(x) {x}, ...
) {
  ## Argument matching
  null_model <- match.arg(null_model)

  if (missing(fitResult) & missing(predSurface)) 
    stop("Marginals fit result or predicted surface need to be supplied.")

  if (is.character(colorBy) & all(colorBy %in% colors())) {
    colorPalette <- colorBy
    colorBy <- "colors"
  }

  ## Calculate extra arguments
  uniqueDoses <- with(data, list("d1" = sort(unique(d1)), "d2" = sort(unique(d2))))
  doseGrid <- expand.grid(uniqueDoses)
  logT <- function(z) log(z + 0.5 * min(z[z != 0]))
  log10T <- function(z) log10(z + 0.5 * min(z[z != 0]))
  ## Transform function for doses
  transformF <- if (logScale) log10T else function(z) z

  zGrid <- predSurface
  ## If marginal fit information is provided, response surface can be
  ## automatically calculated.
  if (is.null(predSurface)) {
    respSurface <- predictResponseSurface(doseGrid, fitResult,
                                  null_model = null_model,
                                  transforms = transforms)
    if (!is.null(transforms)) {
      predSurface <- with(transforms,
                          InvPowerT(respSurface, compositeArgs))
    } else {
      predSurface <- respSurface
    }
    zGrid <- predSurface
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
    if (length(pCols) > 1) pCols <- min(pCols)

    if (any(duplicated(coloredBy[, c("d1", "d2")]))) {
      coloredBy <- aggregate(coloredBy[, pCols],
                             by = coloredBy[, c("d1", "d2")], FUN = colorfun)
    }

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

  surfaceColors <- colorRampPalette(colorPalette)(length(breaks) - 1)
  
  getFF = function(response){
    if(is.numeric(response)) {
      cut(response, breaks = breaks, include.lowest = TRUE)
    } else if(is.factor(response)){
      response
    } else if(is.character(response)){
      factor(response, levels = c("Syn", "None", "Ant"), 
             labels = c("Syn", "None", "Ant"),
             ordered = TRUE)
    }
  }
  surfaceColor <- function(response) {
    ff <- getFF(response)
    zcol <- surfaceColors[ff]
    return(zcol)
  }

  getLabels <- function(response) {
    ff <- getFF(response)
    labels <- gsub(",", ", ", levels(ff))
    return(labels)
  }

  ## Generate colors for the surface plot
  if (colorBy == "asis") {
    if(is.numeric(colorVec))
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
  labnames <- c("Response", if (!is.null(fitResult$names))
            fitResult$names else c("Compound 1", "Compound 2"))
  if (!is.null(attr(data, "orig.colnames")))
    labnames <- attr(data, "orig.colnames")
  

  ## Plot with no axes/labels
  if (!is.null(plotfun))
    data <- aggregate(effect ~ d1 + d2, data, FUN = plotfun)[, names(data)]
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

  MatrixForColor <- matrix(as.numeric(factor(zcol, levels = unique(surfaceColors))), nrow = nrow(zGrid), ncol = ncol(zGrid))
  
  if (any(is.na(MatrixForColor))) {
    if (missing(colorPaletteNA)) colorPaletteNA <- "grey70"
    if (!colorPaletteNA %in% colorPaletteNA) stop("Please indicate a colour for `colorPaletteNA` that is present in the `colorPalette` list")
    MatrixForColor[is.na(MatrixForColor)] <- which(colorPalette %in% colorPaletteNA)[1]
  }
  
  # Layout setting (x, y and z axis)
  axx <- list(
    backgroundcolor = "rgb(250, 250, 250)",
    gridcolor = "rgb(150, 150, 150)",
    showbackground = TRUE,
    ticketmode = 'array',
    ticktext = digitsFunc(as.numeric(xlab)),
    tickvals = xat,
    title = fitResult$names[1]
  )
  
  axy <- list(
    backgroundcolor = "rgb(250, 250, 250)",
    gridcolor = "rgb(150, 150, 150)",
    showbackground = TRUE,
    ticketmode = 'array',
    ticktext = digitsFunc(as.numeric(ylab)),
    tickvals = yat,
    title = fitResult$names[2]
  )
  
  axz <- list(
    backgroundcolor = "rgb(250, 250, 250)",
    gridcolor = "rgb(150, 150, 150)",
    showbackground = TRUE,
    title = "Response"
  )
  
  data$color <- colorPoints[1 + 1 * (data$d2 == 0) + 2 * (data$d1 == 0)]
  data$color <- factor(data$color, levels = colorPoints)


  # Plot the main surface
  p <- plotly::plot_ly(
    height = height, width = width, colors = unique(colorPalette), 
    type = "surface", showlegend = FALSE
  )
  
  if (gradient) {
    p <- plotly::add_surface(
      p,
      x = transformF(uniqueDoses$d1),
      y = transformF(uniqueDoses$d2),
      #NOTE: Plotly requires to transpose the zGrid matrix (and consequently the MatrixForColor)
      z = zTransform(t(as.matrix(zGrid))),
      opacity = 0.8,
      surfacecolor = t(MatrixForColor),
      cauto = FALSE,
      colors = unique(surfaceColors),
      cmin = 1,
      cmax = length(unique(surfaceColors)),
      text = "",
      hoverinfo = 'text',
      showlegend = FALSE, showscale = FALSE
    ) 
  } else {
    
    colorPaletteHex <- col2hex(unique(colorPalette))
    colorVecAlt <- colorPaletteHex[sort(unique(as.numeric(MatrixForColor)))]
    
    mz <- seq(0,1, length.out = length(colorVecAlt)+1)
    if (length(mz) > 2)
      mz <- c(mz[1], rep(mz[2:(length(mz)-1)], each = 2), mz[length(mz)])
    custom_colors <- data.frame(
      z   =  mz,
      col = rep(colorVecAlt, each = 2)
    )
    
    p <- plotly::add_surface(
      p,
      x = transformF(uniqueDoses$d1),
      y = transformF(uniqueDoses$d2),
      #NOTE: Plotly requires to transpose the zGrid matrix (and consequently the MatrixForColor)
      z = zTransform(t(as.matrix(zGrid))),
      # zauto = FALSE,
      # opacity = 0.8,
      surfacecolor = t(MatrixForColor),
      colorscale = custom_colors,
      # cauto = FALSE,
      text = "",
      hoverinfo = 'text',
      showlegend = FALSE, showscale = FALSE
    )
    
  }
    
  if (legend) {
    # Categorical legend (only if colorPalette was named)
    if (!is.null(names(colorPalette))) {
      uColorPalette <- colorPalette[!duplicated(colorPalette)]
      for (xcolor in names(uColorPalette)) {
        p <- add_markers(
          p, x = -900, y = -900, color = I(uColorPalette[as.character(xcolor)]), size = 10,
          legendgroup = "call", name = xcolor, opacity = 1, 
          showlegend = TRUE, marker = list(symbol = "square")
        )
      }

      p <- layout(
        p,
        xaxis = list(
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE
        ),
        yaxis = list(
          showgrid = FALSE,
          zeroline = FALSE,
          showticklabels = FALSE
        ),
        legend = list(
          title = "Call:",
          itemsizing = 'constant',
          font = list(size = 10),
          orientation = "h",   # show entries horizontally
          xanchor = "center",  # use center of legend as anchor
          x = 0.5,
          y = 0
        )
      )
    } else {
      #TODO Add legend for continuous variables
    }
  }
  
  if (addPoints) {
    # Legend names
    pointNames <- c(
      paste0("Combination \n", paste(fitResult$names, collapse = " and ")), 
      fitResult$names, 
      "Dose 0 for both compounds"
    )
    names(pointNames) <- colorPoints
    
    ## add points (spheres)
    for (xcolor in levels(data$color)) {
      df <- data.frame(
        d1_transf = transformF(data$d1)[data$color == xcolor],
        d1        = data$d1[data$color == xcolor],
        d2_transf = transformF(data$d2)[data$color == xcolor],
        d2        = data$d2[data$color == xcolor],
        effect    = zTransform(data$effect)[data$color == xcolor]
      )
      p <- plotly::add_markers(
        p = p,
        opacity = 1,
        x = df$d1_transf,
        y = df$d2_transf,
        z = df$effect,
        marker = list(color = xcolor, size = radius,  line = list(color = "#333", width = 1)),
        showlegend = legend,
        text = paste0("d1: ", digitsFunc(df$d1), "\nd2: ", digitsFunc(df$d2), "\neffect: ", digitsFunc(df$effect)),
        hoverinfo = "text",
        name = pointNames[xcolor],
        legendgroup = "points"
      ) 
    }

    # set legend layout
    p <- plotly::layout(
      p = p,
      title = title,
      scene = list(
        xaxis = axx, yaxis = axy, zaxis = axz,
        aspectratio = list(x = 1, y = 1, z = 1.1)
      ),
      legend = list(
        itemsizing = 'constant',
        font = list(size = 10),
        orientation = "h",   # show entries horizontally
        xanchor = "center",  # use center of legend as anchor
        x = 0.5,
        y = 0
      )
    )
  }
  
  
  p <- layout(
    p,
    scene = list(
      xaxis = axx, yaxis = axy, zaxis = axz,
      aspectratio = list(x = 1, y = 1, z = 1.1)
    ),
    title = title
  )
  
  p <- config(
    p,
    displaylogo = FALSE,
    modeBarButtonsToRemove = c("orbitRotation", "resetCameraLastSave3d", "hoverClosest3d")
  )

  p

}



