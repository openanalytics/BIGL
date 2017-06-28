#' Highest Single Agent model
#'
#' This function returns response levels for when these are based on
#' Highest Single Agent (HSA) model.
#'
#' @inheritParams generalizedLoewe
hsa <- function(doseInput, parmInput, ...) {

  pars <- parmInput
  increasing <- pars["m1"] > pars["b"] & pars["m2"] > pars["b"]
  decreasing <- pars["m1"] < pars["b"] & pars["m2"] < pars["b"]

  ## If agonist and antagonist, try to determine the leading compound and emit a
  ## warning.
  if (!(increasing | decreasing)) {
    warning("Marginal curves are diverging. HSA might be flawed.")
    lead <- which.max(c(abs(pars["m1"] - pars["b"]),
                        abs(pars["m2"] - pars["b"])))
    if (lead == 1)
      applyFunction <- if (pars["m1"] > pars["b"]) max else min
    else if (lead == 2)
      applyFunction <-  if (pars["m2"] > pars["b"]) max else min
  } else {
    applyFunction <- if (increasing) max else min
  }

  pred1 <- L4(doseInput[, "d1"], pars["h1"], pars["b"], pars["m1"], pars["e1"])
  pred2 <- L4(doseInput[, "d2"], pars["h2"], pars["b"], pars["m2"], pars["e2"])

  apply(cbind(pred1, pred2), 1, applyFunction)

}
