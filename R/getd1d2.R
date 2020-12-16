#' A function to get the d1d2 identifier
#' @param dat the data frame containing d1 and d2 entries
#' @return a vector of d1d2 identifiers
getd1d2 = function(dat){
    apply(dat[, c("d1", "d2")], 1, paste, collapse = "_")
}