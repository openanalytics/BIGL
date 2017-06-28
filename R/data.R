#' Partial data with combination experiments of direct-acting antivirals
#'
#' A dataset containing 11 combination experiments of direct-acting antivirals.
#'
#' @name directAntivirals
#' @docType data
#' @format A data frame with 3520 rows and 6 variables:
#' \itemize{
#'   \item experiment: ID of experiment (1-11)
#'   \item cpd1: name of the first compound (4 different compounds)
#'   \item cpd2: name of the second compound (11 different compounds)
#'   \item effect: observed effect (cell count)
#'   \item d1: dose of the first compound
#'   \item d2: dose of the second compound
#' }
NULL

#' Full data with combination experiments of direct-acting antivirals
#'
#' A dataset containing 11 combination experiments of direct-acting antivirals.
#' This dataset is larger than \code{directAntivirals} dataset as it includes
#' concentrations at levels of \code{1e6} which can render plots visually
#' unappealing.
#'
#' @name directAntivirals_ALL
#' @docType data
#' @format A data frame with 4224 rows and 6 variables:
#' \itemize{
#'   \item experiment: ID of experiment (1-11)
#'   \item cpd1: name of the first compound (4 different compounds)
#'   \item cpd2: name of the second compound (11 different compounds)
#'   \item effect: observed effect (cell count)
#'   \item d1: dose of the first compound
#'   \item d2: dose of the second compound
#' }
NULL
