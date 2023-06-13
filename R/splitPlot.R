globalVariables(c("loewe", "loewe2", "hsa", "bliss", "predicted", "d2", "effect"))

#' Plot 2D cross section of response surface
#' @param ls list of results objects obtained from \code{\link{fitSurface}}. Names of list objects 
#' expected to be one of the null model options i.e. loewe, loewe2, hsa, bliss
#' @param xlab label for x-axis
#' @param ylab label for y-axis
#' @param color plot lines in colour? Defaults to FALSE
#' @param plotBy compound name to be used for order of plotting. If plotBy = "Compound 1" then plots are split by 
#' concentrations in Compound 1 and concentrations in Compound 2 are shown on the x-axis.
#' @author Mohammed Ibrahim
#' @examples
#' \dontrun{
#'   data <- subset(directAntivirals, experiment == 1)
#'   transforms <- list("PowerT" = function(x, args) with(args, log(x)),
#'                      "InvPowerT" = function(y, args) with(args, exp(y)),
#'                      "BiolT" = function(x, args) with(args, N0 * exp(x * time.hours)),
#'                      "InvBiolT" = function(y, args) with(args, 1/time.hours * log(y/N0)),
#'                      "compositeArgs" = list(N0 = 1, time.hours = 72))
#'   fitResult <- fitMarginals(data, transforms)
#'   nullModels <- c("loewe", "loewe2", "bliss", "hsa")
#'   rs_list <- Map(fitSurface, null_model = nullModels, MoreArgs = list(
#'       data = data, fitResult = fitResult, B.CP = 50, statistic = "none"))
#'   synergy_plot_bycomp(ls = rs_list, plotBy = "Compound 1", color = TRUE)
#'   synergy_plot_bycomp(ls = rs_list, plotBy = "Compound 2", color = TRUE)
#' } 
#' @export
synergy_plot_bycomp <- function(ls, xlab = NULL, ylab = NULL, color = FALSE, plotBy = NULL) {
	
  ls <- Filter(function(x) { !inherits(x, "try-error") && !is.null(x) }, ls)
	
	nmes <- names(ls)
	tmp <- lapply(nmes, function(nn) {
	  res <- ls[[nn]]$offAxisTable[, c("d1", "d2", "effect", "predicted")]
	  colnames(res)[colnames(res) == "predicted"] <- nn
	  res
	})
	
	plot_df <- tmp[[1]]
	
	if (length(tmp) > 1) {
		for (i in 2:length(tmp)) {
		  plot_df <- merge(plot_df, tmp[[i]], by = c("d1", "d2", "effect"))
		}
	}
	
	if (is.null(plotBy)) {
		plotBy <- ls[[1]]$names[1]
  } else {
    if (!plotBy %in% ls[[1]]$names) {
      warning("Unrecognized name in `plotBy`, ", toString(ls[[1]]$names[[1]]), " will be used.")
      plotBy <- ls[[1]]$names[1]
    } else if (plotBy == ls[[1]]$names[2]) {
      names(plot_df)[1:2] <- names(plot_df)[2:1]
      if (is.null(xlab)) {
        xlab <- paste(ls[[1]]$names[1], "Concentration")
      }
    }
  }

  if (is.null(xlab)) {
    xlab <- paste(ls[[1]]$names[2], "Concentration")
  }
  
  
	facet_d1 <- unique(plot_df$d1)
	facet_d1 <- facet_d1[order(facet_d1)]
	
	plot_df$d1 <- factor(plot_df$d1, levels = facet_d1, labels = formatC(facet_d1, format = "fg", digits = 2))
	
	allLabels <- c("Bliss" = "bliss", "HSA" = "hsa", "Generalized Loewe" = "loewe", "Alternative Loewe" = "loewe2")
	allValues <- c("Bliss" = 1, "HSA" = 2, "Generalized Loewe" = 3, "Alternative Loewe" = 4)
	
	labels <- allLabels[which(allLabels %in% nmes)]
	values <- allValues[which(allLabels %in% nmes)]
	
	breaks <- unique(plot_df$d2)
	breakLabels <- trimws(formatC(breaks, format = "fg", digits = 1))
	
	if (color) {
		
		color_vec <- c(
				"loewe2" = "gold3",         #"Alternative Loewe"
				"bliss"  = "green3",        #"Bliss"
				"loewe"  = "coral2",        #"Generalized Loewe"
				"hsa"    = "cornflowerblue" #"HSA"
		)
		color_vec <- color_vec[names(color_vec) %in% labels]
		labels    <- labels[match(names(color_vec), labels)]
		
		p <- ggplot(plot_df, aes(x = d2, y = effect)) +
				geom_point() +
				facet_wrap(~d1) +
				scale_x_continuous(trans = "log10", breaks = breaks, labels = breakLabels) +
				#scale_color_manual(values = values, labels = names(labels)) +
				scale_color_manual(values = as.character(color_vec[labels]), labels = names(labels)) +
				theme_bw() +
				theme(
				  panel.grid = element_blank(), 
				  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
				  legend.position = "top", legend.justification = c(0.5,0),
					legend.title = element_text(hjust = 0.5), 
					legend.background = element_rect(linewidth = 0.25, color = "black")) +
				labs(x = xlab, y = ylab, color = "Model") +
				guides(color = guide_legend(nrow = 1, byrow = T))
		
		if ("bliss" %in% labels) {
			p <- p + geom_line(aes(y = bliss, color = "Bliss"))
		}
		if ("hsa" %in% labels) {
			p <- p + geom_line(aes(y = hsa, color = "HSA"))
		}
		if ("loewe" %in% labels) {
			p <- p + geom_line(aes(y = loewe, color = "Generalized Loewe"))
		}
		if ("loewe2" %in% labels) {
			p <- p + geom_line(aes(y = loewe2, color = "Alternative Loewe"))
		}
		
	} else {
		
		line_vec <- c(
				"Alternative Loewe" = 1, #"loewe2"
				"Bliss"             = 2, #"bliss"
				"Generalized Loewe" = 3, #"loewe"
				"HSA"               = 4  #"hsa"
		)
		line_vec <- line_vec[names(line_vec) %in% names(labels)]
		labels    <- labels[match(names(line_vec), names(labels))]
		
		p <- ggplot(plot_df, aes(x = d2, y = effect)) +
				geom_point() +
				facet_wrap(~d1) +
				scale_x_continuous(trans = "log10", breaks = breaks, labels = breakLabels) +
				scale_linetype_manual(values = line_vec, labels = names(labels)) +
				theme_bw() +
				theme(
				  panel.grid = element_blank(), 
				  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
				  legend.position = "top", legend.justification = c(0.5,0),
					legend.title = element_text(hjust = 0.5), 
					legend.background = element_rect(linewidth = 0.25, color = "black")
				) +
				labs(x = xlab, y = ylab, linetype = "Model") +
				guides(linetype = guide_legend(nrow = 1, byrow = T))
		
		if ("bliss" %in% labels) {
			p <- p + geom_line(aes(y = bliss, linetype = "Bliss"))
		}
		if ("hsa" %in% labels) {
			p <- p + geom_line(aes(y = hsa, linetype = "HSA"))
		}
		if ("loewe" %in% labels) {
			p <- p + geom_line(aes(y = loewe, linetype = "Generalized Loewe"))
		}
		if ("loewe2" %in% labels) {
			p <- p + geom_line(aes(y = loewe2, linetype = "Alternative Loewe"))
		}
	}
	
	p
	
}