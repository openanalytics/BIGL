#' Bliss Independence Model
#' 
#' This function returns fractional response levels for when these are based on
#' Bliss Independence Model.
#' 
#' @param doseInput Dose-response dataframe containing \code{"d1"} and
#'   \code{"d2"} columns
#' @param parmInput Numeric vector or list with appropriately named
#'   parameter inputs. Typically, it will be coefficients from a
#'   \code{MarginalFit} object.
#' @param ... Further arguments that are currently unused

Blissindependence <- function(doseInput, parmInput, ...){
  
  pars <- parmInput
  #We need the percentage value between 0 and 1 for Bliss estimated effect
  Hilleq <- function(dose,b,logEC50){
    Hilleqread <- 1/(1 + (exp(logEC50)/dose)^(abs(b)))
  }
  
  #calculate prediction mono and rescale to max upper for percentage
  pred1 <- Hilleq(doseInput[, "d1"], pars["h1"], pars["e1"])
  pred1 <- pred1*abs(pars["m1"]-pars["b"])/max(c(abs(pars["m1"]-pars["b"]),abs(pars["m2"]-pars["b"])))
  pred2 <- Hilleq(doseInput[, "d2"], pars["h2"], pars["e2"])
  pred2 <- pred2*abs(pars["m2"]-pars["b"])/max(c(abs(pars["m1"]-pars["b"]),abs(pars["m2"]-pars["b"])))
  
  #Bliss independence combination
  applyfunction <- function(param){
    #prediction value in percentage
    predcombo <- param[1]+param[2]-param[1]*param[2]
    #determine whether curves are rising or lowering
    rise = 1
    if(param[3]>param[4] & param[3] > param[5]){
      rise = -1
    }
    predcombo <- param[3]+rise*max(c(abs(param[4]-param[3]),abs(param[5]-param[3])))*predcombo
    return(predcombo)
  }
  
  apply(cbind(pred1,pred2,pars["b"], pars["m1"],pars["m2"]),1, applyfunction)
}
