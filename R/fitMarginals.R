#' Fit two 4-parameter log-logistic functions for a synergy experiment
#'
#' This function uses dose-response data for two compounds and estimates
#' coefficients for monotherapy models of both of these compounds such that they
#' share a common baseline. Currently, these coefficients are estimated by
#' default using a non-linear least squares approximation. Although entire
#' dose-response data can be provided, estimation will subset the part of data
#' where at least one of the compounds is dosed at zero, i.e. on-axis data.
#'
#' Model formula is specified as \code{effect ~ fn(h1, h2, ...)} where \code{fn}
#' is a hard-coded function which fits two 4-parameter log-logistic functions
#' simultaneously so that the baseline can be shared. If transformation
#' functions are provided, \code{fn} is consequently adjusted to account for
#' them.
#'
#' @param data Dose-response dataframe. Marginal data will be extracted from
#'   it automatically.
#' @param start Starting parameter values. If not specified, they will be
#'   obtained from \code{\link{initialMarginal}}.
#' @param constraints List of constraint matrix and vector which will be passed
#'   to \code{\link{constructFormula}}. If \code{constraints = NULL}, no
#'   constraints on parameter estimation will be imposed.
#' @param fixed This arguments provides a user-friendly alternative to impose a
#'   fixed value for marginal parameters. It must be a named vector with names
#'   contained in \code{c("h1", "h2", "b", "m1", "m2", "e1", "e2")}. For
#'   example, \code{fixed = c("m1" = 1, "h1" = 1)} will automatically generate
#'   appropriate constraint matrix and vector to set the maximal response and
#'   the Hill coefficient of the first compound to 1. If both \code{constraints}
#'   and \code{fixed} arguments are passed, then only \code{fixed} will be used.
#' @param method Which estimation method should be used to obtain the estimates.
#'   If \code{method = "nls"}, simple non-linear least squares
#'   \code{\link[stats]{nls}} will be used. If \code{method = "nlslm"}
#'   Levenberg-Marquardt non-linear least squares
#'   \code{\link[minpack.lm]{nlsLM}} is used instead (default). If \code{method
#'   = "optim"}, residual sum of squares will be minimized using general purpose
#'   optimization based on Nelder-Mean algorithm in \code{\link[stats]{optim}}.
#'   This method can be noticeably slower than the non-linear least squares
#'   methods.
#' @param ... Further arguments that are passed to the optimizer function, 
#' such as \code{lower} or \code{upper} (for the "nlslm" method), or 
#' \code{control}.
#' @inheritParams fitSurface
#' @importFrom methods hasArg
#' @importFrom minpack.lm nlsLM
#' @importFrom stats nls
#' @importFrom utils modifyList
#' @return This function returns a \code{MarginalFit} object with monotherapy
#'   coefficient estimates and diverse information regarding monotherapy
#'   estimation. \code{MarginalFit} object is essentially a list with
#'   appropriately named elements.
#'
#'   Among these list elements, \code{"coef"} is a named vector with parameter
#'   estimates. \code{h1} and \code{h2} are Hill's slope coefficients for each
#'   of the compounds, \code{m1} and \code{m2} are their maximal response levels
#'   whereas \code{b} is the shared baseline. Lastly, \code{e1} and \code{e2}
#'   are log-transformed EC50 values.
#'
#'   \code{"sigma"} is standard deviation of residuals for the estimated
#'   monotherapy model and \code{"df"} is the degrees of freedom for the
#'   residuals. \code{"vcov"} is the variance-covariance matrix of the estimated
#'   parameters.
#'
#'   Return object also contains information regarding data, biological and
#'   power transformations used in this estimation as well as model construct
#'   and method of estimation.
#' @examples
#'   data <- subset(directAntivirals, experiment == 1)
#'   ## Data must contain d1, d2 and effect columns
#'   transforms <- getTransformations(data)
#'   fitMarginals(data, transforms)
#' @export
fitMarginals <- function(data, transforms = NULL, start = NULL,
                         constraints = NULL, fixed = NULL,
                         method = c("nlslm", "nls", "optim"), ...) {

  method <- match.arg(method)

  ## Verify column names of input dataframe
  if (!all(c("effect", "d1", "d2") %in% colnames(data)))
    stop("effect, d1 and d2 arguments must be column names of data")

  ## Keep only marginal data
  data <- data[(abs(data$d1) < .Machine$double.eps |
                abs(data$d2) < .Machine$double.eps), ]

  fitArgs <- list("data" = data,
                  "transforms" = transforms)

  ## If no starting parameters are provided, make an initial guess
  if (is.null(start))
    start <- initialMarginal(data, transforms)

  ## If model argument is provided in the ellipsis, re-use it instead of
  ## recalculating everything based on constraints.
  if (!hasArg("model") || is.null(list(...)$model)) {

    if (!is.null(fixed)) {

      if (!is.null(constraints))
        warning("Both `fixed` and `constraints` parameters were specified. Only `fixed` will be used.")

      fixedInd <- match(names(fixed), c("h1", "h2", "b", "m1", "m2", "e1", "e2"))
      if (any(is.na(fixedInd))) {
        stop(paste0("Please use parameter names: h1, h2, b, m1, m2, e1 or m2, ",
                    "or construct the formula using constructFormula()."))
      }

      constraints <- list("matrix" = diag(1, 7)[fixedInd,],
                          "vector" = fixed)
    }

    ## vars_x <- colnames(data)[!grepl("effect", colnames(data))]
    fitArgs$model <- {
      if (is.null(constraints))
        ## No constraints provided
        constructFormula()
      else
        ## With constraints
        constructFormula(constraints$matrix, constraints$vector)
    }
  } else {
    fitArgs$model <- list(...)$model
  }

  ## Pass ... along, except for 'model' :-/ 
  ## FIXME by renaming 'model' argument here and in marginalNLS
  extraArgs <- list(...)
  extraArgs <- extraArgs[names(extraArgs) != "model"]
  fitArgs <- modifyList(fitArgs, extraArgs, keep.null = TRUE)
  
  ## Subset only free parameters
  fitArgs$start <- start[fitArgs$model$free]

  fitResult <- switch(method,
                      "nlslm" = { fitArgs$nlsfn <- nlsLM; do.call(marginalNLS, fitArgs) },
                      "nls" = { fitArgs$nlsfn <- nls; do.call(marginalNLS, fitArgs) },
                      "optim" = { do.call(marginalOptim, fitArgs) })

  fitResult$method <- method
  class(fitResult) <- append("MarginalFit", class(fitResult))

  fitResult
}


#' Construct a model formula from parameter constraint matrix
#'
#' For parameter names defined in \code{naming} vector, formula is constructed
#' so that \code{consMatrix \%*\% naming = consVector} is satisfied. Constraint
#' coefficients are normalized and convert into fractions.
#'
#' @param consMatrix Constraint matrix
#' @param consVector Constraint vector
#' @param naming Parameter names
#' @param extraVars Non-parameter variables used in the formula and function
#'   evaluation. These will be appended to the formula.
#' @param formulaArgs Character vector of length two. First element indicates
#'   name for the response variable. Second element indicates name of the
#'   function.
#' @importFrom MASS fractions
#' @return This function returns a model construct appropriate for
#'   \code{\link{fitMarginals}} function. It also separates variables into those
#'   that are free and those which are constrained.
#' @export
#' @examples
#'   constM <- rbind(c(0, 0, 1, 0, 0, 0, 0),
#'                   c(0, 0, 0, -1, 1, 0, 0))
#'   constV <- c(0.9, 0)
#'   constructFormula(constM, constV)
constructFormula <- function(consMatrix = NULL, consVector = NULL,
                             naming = c("h1", "h2", "b", "m1", "m2", "e1", "e2"),
                             extraVars = c("d1", "d2"),
                             formulaArgs = c("effect", "fn")) {

  ## If no constraints are provided, return a default formula
  if (is.null(consMatrix) | is.null(consVector)) {
    model <- list("formula" = paste0(formulaArgs[1], " ~ ", formulaArgs[2], "(",
                                     paste(c(naming, extraVars), collapse = ", "), ")"),
                  "free" = naming,
                  "vars" = extraVars,
                  "order" = naming)
    return(model)
  }

  if (inherits(consMatrix, "numeric"))
    consMatrix <- matrix(consMatrix, nrow = 1, ncol = length(consMatrix))

  if (ncol(consMatrix) < length(naming))
    stop("Constraint matrix does not have enough columns.")
  if (nrow(consMatrix) != length(consVector))
    stop("Number of constraints in the matrix and vector do not match.")

  consIndex <- apply(consMatrix, 1, function(x) max(which(x != 0)))
  consFactors <- apply(consMatrix, 1, function(x) 1 / x[max(which(x != 0))])

  if (anyDuplicated(consIndex))
    stop(naming[consIndex[duplicated(consIndex)]], " cannot be constrained twice.")

  ## Normalize constraint constants
  normMatrix <- fractions(diag(consFactors, length(consFactors)) %*% consMatrix)
  normVector <- fractions(consFactors * consVector)

  namingOrig <- naming

  for (i in seq_along(consIndex)) {

    consVar <- consIndex[i]
    naming[consIndex[i]] <- paste0(naming[consVar], " = ", normVector[i])

    ## Check if value is constrained in terms of other variables
    otherVar <- setdiff(which(normMatrix[i,] != 0), consVar)
    if (length(otherVar) > 0) {
      naming[consVar] <- paste0(naming[consVar], " - ",
                                paste(paste(normMatrix[i, otherVar],
                                            naming[otherVar], sep = "*"),
                                      collapse = " + "))
    }
  }

  ## Make formula look nicer by doing some substitution
  naming <- gsub("- -", "+ ", naming)
  naming <- gsub("= 0 (\\+|\\-)", "=", naming)
  naming <- gsub(" \\+ 1[:space:]*\\*", " +", naming)
  naming <- gsub(" - 1[:space:]*\\*", " -", naming)
  naming <- gsub("= 1\\*", "= ", naming)

  list("formula" = paste0(formulaArgs[1], " ~ ", formulaArgs[2], "(",
                          paste(c(naming, extraVars), collapse = ", "), ")"),
       "free" = namingOrig[-consIndex],
       "nonfree" = namingOrig[consIndex],
       "order" = namingOrig,
       "vars" = extraVars,
       "constraints" = list("matrix" = normMatrix,
                            "vector" = normVector))

}


#' Fit two 4-parameter log-logistic functions with non-linear least squares
#'
#' This function does not automatically extract marginal data and requires
#' model input obtained from \code{\link{constructFormula}}.
#'
#' @param model List with model parameters. Typically, this is an output from
#'   \code{\link{constructFormula}}.
#' @param nlsfn Non-linear least-squares optimizer function
#' @importFrom methods formalArgs
#' @importFrom stats as.formula coef df.residual vcov
#' @inheritParams fitMarginals
marginalNLS <- function(data, transforms = NULL, start, model,
                        nlsfn = nls, ...) {

  dataU <- data

  if (is.null(transforms)) {
    PowerT <- function(x) x
    BiolT <- function(x) x
  } else {
    PowerT <- function(x)
      transforms$PowerT(x, transforms$compositeArgs)
    BiolT <- function(x)
      transforms$BiolT(x, transforms$compositeArgs)
    ## Power-transform response
    data$effect <- PowerT(data$effect)
  }

  fn <- function(h1, h2, b, m1, m2, e1, e2, d1, d2) {
    PowerT(BiolT(ifelse(d2 == 0,
                        L4(d1, h1, b, m1, e1),
                        L4(d2, h2, b, m2, e2))))
  }

  ## Starting parameters have to be in a list
  start <- as.list(start)
  nlsArgs <- list("formula" = as.formula(model$formula),
                  "start" = start, "data" = data)

  ## Pick out args for NLS function from all ellipsis arguments
  extraArgs <- as.list(substitute(list(...)))[-1L]
  nlsArgs <- c(nlsArgs, extraArgs[names(extraArgs) %in% formalArgs(nlsfn)])

  fit <- do.call(nlsfn, args = nlsArgs)

  ## If there were any constraints, coefficients need to be re-arranged
  if (is.null(model$constraints)) {
    coefs <- coef(fit)
  } else {
    ## Rearrange parameters
    coefs <- rep(0, length(model$order))
    names(coefs) <- model$order
    coefs[match(names(coef(fit)), names(coefs))] <- coef(fit)
    coefs[model$nonfree] <- model$constraints$vector - model$constraints$matrix %*% coefs
  }

  list("coef" = coefs,
       "sigma" = summary(fit)$sigma,
       "df" = df.residual(fit),
       "data" = dataU,
       "transforms" = transforms,
       "vcov" = vcov(fit),
       "model" = model,
       "shared_asymptote" = as.logical(coefs["m1"] == coefs["m2"]),
       "extraArgs" = extraArgs)
}

#' Fit two 4-parameter log-logistic functions with common baseline
#'
#' This function is an alternative to non-linear least squares and
#' provides optimization framework with \code{\link{optim}} function.
#' It is however noticeably slower than NLS methods and can be especially
#' time consuming in large datasets, in particular if bootstrap statistics
#' are calculated.
#'
#'
#' @param ... Further parameters passed to \code{\link[stats]{optim}} function
#' @inheritParams fitMarginals
#' @inheritParams marginalNLS
#' @return Variance-covariance matrix which is returned by \code{\link{optim}}
#'   is based on the fact that minimization of sum-of-squared residuals leads
#'   essentially to a maximum likelihood estimator and so variance-covariance
#'   matrix can be estimated using inverse Hessian evaluated at the optimal
#'   parameters. In some cases, so obtained variance-covariance matrix might not
#'   be positive-definite which probably means that estimates are unstable
#'   because of either a poor choice of initial values or poor properties of the
#'   data itself.
#' @importFrom numDeriv hessian
#' @importFrom stats optim
marginalOptim <- function(data, transforms = NULL, start, model, ...) {

  dataU <- data

  if (is.null(transforms)) {
    PowerT <- function(x) x
    BiolT <- function(x) x
  } else {
    PowerT <- function(x)
      transforms$PowerT(x, transforms$compositeArgs)
    BiolT <- function(x)
      transforms$BiolT(x, transforms$compositeArgs)
    ## Power-transform response
    data$effect <- PowerT(data$effect)
  }

  fn <- function(pars) {
    pred <- apply(data[, c("d1", "d2")], 1, function(x) {
      PowerT(BiolT(ifelse(x[2] == 0,
                          L4(x[1], pars["h1"], pars["b"], pars["m1"], pars["e1"]),
                          L4(x[2], pars["h2"], pars["b"], pars["m2"], pars["e2"]))))
    })

    data$effect - pred
  }

  ## Starting parameters have to be in a list
  start <- as.list(start)

  opt_fn <- function(pars) {

    ## Rearrange parameters in case of constraints
    if (!is.null(model$constraints)) {
      pars_new <- rep(0, length(model$order))
      names(pars_new) <- model$order
      pars_new[match(names(pars), names(pars_new))] <- pars
      pars_new[model$nonfree] <- model$constraints$vector - model$constraints$matrix %*% pars_new
      pars <- pars_new
    }

    sum(fn(pars)^2) / 2
  }

  fit <- optim(start, opt_fn, ...)

  ## If there were any constraints, coefficients need to be re-arranged
  if (!is.null(model$constraints)) {
    ## Rearrange parameters
    coefs <- rep(0, length(model$order))
    names(coefs) <- model$order
    coefs[match(model$free, names(coefs))] <- fit$par
    coefs[model$nonfree] <- model$constraints$vector - model$constraints$matrix %*% coefs
  } else {
    coefs <- fit$par
  }

  resid <- fn(coefs)

  df <- length(resid) - length(coefs)
  sigma.sq <- sum(resid^2) / df

  vcov <- tryCatch({
    sigma.sq * solve(hessian(fn, fit$par))
  }, error = function(e) "Failed to compute inverse Hessian.")

  list("coef" = coefs,
       "sigma" = sqrt(sigma.sq),
       "df" = df,
       "data" = dataU,
       "transforms" = transforms,
       "vcov" = vcov,
       "model" = model,
       "shared_asymptote" = as.logical(coefs["m1"] == coefs["m2"]),
       "extraArgs" = as.list(substitute(list(...)))[-1L]   
   )
}
