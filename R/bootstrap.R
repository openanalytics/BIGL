#' Generate data from parameters of marginal monotherapy model
#'
#' This function is used to generate data for bootstrapping of the null
#' distribution for various estimates. Optional arguments such as specific
#' choice of sampling vector or corrections for heteroskedasticity can be
#' specified in the function arguments.
#'
#' @param pars Coefficients of the marginal model along with their appropriate
#'   naming scheme. These will typically be estimated using
#'   \code{\link{fitMarginals}}. Futhermore, \code{pars} can simply be a
#'   \code{MarginalFit} object and \code{transforms} object will be
#'   automatically extracted.
#' @param sigma Standard deviation to use for randomly generated error terms. This
#'   argument is unused if \code{error = 4} so that sampling error vector is
#'   provided.
#' @param data Data frame with dose columns \code{("d1", "d2")} to generate the
#'   effect for. Only \code{"d1"} and \code{"d2"} columns of the dose-response
#'   dataframe should be passed to this argument. \code{"effect"} column should
#'   not be passed and if it is, the column will be replaced by simulated data.
#' @param null_model Specified null model for the expected response surface.
#'   Currently, allowed options are \code{"loewe"} for generalized Loewe model,
#'   \code{"hsa"} for Highest Single Agent model, \code{"bliss"} for Bliss additivity,
#'   and \code{"loewe2"} for the alternative Loewe generalization.
#' @param error Type of error for resampling. \code{error = 1} (Default) adds
#'   normal errors to the simulated effects, \code{error = 2} adds errors sampled
#'   from a mixture of two normal distributions, \code{error = 3} generates errors
#'   from a rescaled chi-square distribution. \code{error = 4} will use bootstrap.
#'   Choosing this option, the error terms will be resampled from the vector
#'   specified in \code{sampling_errors}.
#' @param sampling_errors Sampling vector to resample errors from. Used only if
#'   \code{error = 4}.
#' @param wild_bootstrap Whether special bootstrap to correct for
#'   heteroskedasticity should be used. If \code{wild_bootstrap = TRUE}, errors
#'   are generated from \code{sampling_errors} multiplied by a random variable
#'   following Rademacher distribution. Argument is used only if \code{error = 4}.
#' @param model The mean-variance model
#' @param means The vector of mean values of the response surface, for variance modelling
#' @param invTransFun the inverse transformation function, back to the variance domain
#' @param ... Further arguments
#' @inheritParams fitSurface
#' @inheritParams predictOffAxis
#' @importFrom stats lm.fit rnorm rchisq rbinom
#' @return Dose-response dataframe with generated data including \code{"effect"}
#'   as well as \code{"d1"} and \code{"d2"} columns.
#' @export
#' @examples
#'   coefs <- c("h1" = 1, "h2" = 1.5, "b" = 0,
#'              "m1" = 1, "m2" = 2, "e1" = 0.5, "e2" = 0.1)
#'
#'   ## Dose levels are set to be integers from 0 to 10
#'   generateData(coefs, sigma = 1)
#'
#'   ## Dose levels are taken from existing dataset with d1 and d2 columns
#'   data <- subset(directAntivirals, experiment == 1)
#'   generateData(data = data[, c("d1", "d2")], pars = coefs, sigma = 1)
generateData <- function(pars, sigma, data = NULL,
		transforms = NULL,
		null_model = c("loewe", "hsa", "bliss", "loewe2"),
		error = 1, sampling_errors = NULL, means = NULL,
		model = NULL, method = "equal", wild_bootstrap = FALSE,
		rescaleResids, invTransFun, newtonRaphson = FALSE, bootmethod = method, ...) {
	
	if(bootmethod == "model") bootmethod <- "unequal"
	
	## Argument matching
	null_model <- match.arg(null_model)
	
	if (inherits(pars, "MarginalFit")) {
		transforms <- pars$transforms
		pars <- pars$coef
	}
	
	if (is.null(data)) data <- expand.grid("d1" = rep(0:10, each = 2),
				"d2" = rep(0:10, each = 2))
	
	if ("effect" %in% colnames(data)) {
		warning("effect column is unneeded for generateData() function and will be dropped.")
		data <- data[, c("d1", "d2")]
	}
	
	## Use identity transformation if no transform functions are supplied
	if (is.null(transforms)) {
		idF <- function(z, ...) z
		transforms <- list("PowerT" = idF,
				"InvPowerT" = idF,
				"BiolT" = idF,
				"compositeArgs" = NULL)
	}
	
	ySim <- switch(null_model,
			"loewe" = generalizedLoewe(data, pars, asymptotes = 2,
					newtonRaphson =  newtonRaphson)$response,
			"hsa" = hsa(data[, c("d1", "d2")], pars),
			"bliss" = Blissindependence(data[, c("d1", "d2")], pars),
			"loewe2" = harbronLoewe(data[, c("d1", "d2")], pars,
					newtonRaphson = newtonRaphson))
	ySim <- with(transforms, PowerT(BiolT(ySim, compositeArgs), compositeArgs))
	charEr = as.character(error)
	if(charEr %in% c("1", "2", "3")){
		errors = switch(charEr,
				## Normal
				"1" = {rnorm(length(ySim), 0, sigma)},
				## Two normals
				"2" = {ru <- sample(seq_len(2), replace = TRUE, size = length(ySim))
					mus <- c(-sigma, sigma)
					sigmas <- c(sigma/2, sigma/3)
					rnorm(length(ySim), mus[ru], sigmas[ru])},
				## Distribution with right-tail outliers
				"3" = {sigma*(rchisq(length(ySim), df = 4)-4)/8 })
	} else if(charEr == "4"){
		if (wild_bootstrap) {
			## Use Rademacher distribution to account for heteroskedasticity
			errors = sampling_errors*(2*rbinom(length(ySim), size = 1, prob = 0.5)-1)
		} else {
			if(bootmethod == "equal"){
				errors = sampleResids(means = ySim, sampling_errors = sampling_errors,
						method = "equal", rescaleResids = FALSE)
			} else {
				idd1d2 = with(data, d1&d2)
				errors = integer(length(ySim))
				#On-axis points
				errors[!idd1d2] = sampleResids(means = ySim[!idd1d2], sampling_errors = sampling_errors[!idd1d2],
						method = "equal", rescaleResids = FALSE)
				#Off-axis points
				errors[idd1d2] = sampleResids(means = ySim[idd1d2], sampling_errors = sampling_errors[idd1d2],
						method = bootmethod, rescaleResids = rescaleResids,
						model = model, invTransFun = invTransFun)
			}
		}
	} else {
		stop("Unavailable error type.")
	}
	ySim <- with(transforms, InvPowerT(ySim + errors, compositeArgs))
	if(all(data$effect>0)){
		ySim = abs(ySim)
	}
	return(data.frame("effect" = ySim, data))
}

#' Estimate CP matrix from bootstraps
#'
#' This function is generally called from within \code{\link{fitSurface}}.
#'
#' @param bootStraps the bootstraps carried out already
#' @param sigma0 standard deviation of the null model on the real data
#' @param doseGrid a grid of dose combinations
#' @inheritParams fitSurface
#' @inheritParams generateData
#' @importFrom stats lm.fit var
#' @return Estimated CP matrix
getCP = function(bootStraps, null_model, transforms, sigma0, doseGrid){
	pred <- vapply(bootStraps,
			FUN.VALUE = bootStraps[[1]]$respS,
			function(b) {b$respS/sigma0})
	var(t(pred))
}
#' Simulate data from a given null model and monotherapy coefficients
#'
#' @param ... Further parameters that will be passed to
#'   \code{\link{generateData}}
#' @param doseGrid A grid of dose combinations
#' @param startvalues Starting values for the non-linear equation,
#'   from the observed data
#' @inheritParams fitSurface
#' @inheritParams generateData
#' @return List with \code{data} element containing simulated data and
#'   \code{fitResult} element containing marginal fit on the simulated data.
#' @export
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   fitResult <- fitMarginals(data)
#'   simDat <- simulateNull(data, fitResult, expand.grid(d1 = data$d1, d2 = data$d2),
#'   null_model = "hsa")
simulateNull <- function(data, fitResult, doseGrid,
		transforms = fitResult$transforms, startvalues,
		null_model = c("loewe", "hsa", "bliss", "loewe2"), ...) {
	## Argument matching
	null_model <- match.arg(null_model)
	
	method <- fitResult$method
	coefFit0 <- fitResult$coef
	sigma0 <- fitResult$sigma
	model <- fitResult$model
	control <- {
		if (method %in% c("nls", "nlslm"))
			list("maxiter" = 200)
	}
	## Parameter estimates may at times return an error due to non-convergence. If
	## necessary, repeat the step until it functions properly and 1000 times at
	## most.
	counter <- 0
	initPars <- coefFit0
	repeat {
		simData <- generateData(pars = coefFit0, sigma = sigma0,
				data = data[, c("d1", "d2")],
				transforms = transforms,
				null_model = null_model, ...)
		## In cases where added errors put the response into negative domain, revert
		## it back to the positive one. Usually, values of such observations tend to
		## be quite small.
		simData$effect <- abs(simData$effect)
		simData$d1d2 = data$d1d2
		
		## construct a list of arguments, including ... passed to original
		## `fitMarginals` call (saved as `extraArgs`)
		paramsMarginal <- list("data" = simData, "method" = method,
				"start" = initPars, "model" = model, "transforms" = transforms,
				"control" = control)
		if (!is.null(fitResult$extraArgs) && is.list(fitResult$extraArgs))
			# use `modifyList` here, since `control` could be user-defined
			paramsMarginal <- modifyList(paramsMarginal, fitResult$extraArgs)
		simFit <- try({
					do.call(fitMarginals, paramsMarginal)
				}, silent = TRUE)
		counter <- counter + 1
		initPars <- NULL
		if (counter > 1000)
			stop(paste("Data simulation process failed. ",
							"Check that transformation functions correspond ",
							"to the marginal model."))
		if (!inherits(simFit, "try-error")) break
	}
	#Also precalculate response surface, quite computation intensive
	respS <- predictOffAxis(fitResult = simFit,
			transforms = transforms, startvalues = startvalues,
			doseGrid = doseGrid, null_model = null_model, ...)
	return(list("data" = simData, "simFit" = simFit, "respS" = respS))
}
bootFun = function(i, args) {
	if(args$progressBar) args$pb$tick()
	do.call(simulateNull, args)
}#Wrapper with index
