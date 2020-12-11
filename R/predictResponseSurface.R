#' Predict the entire response surface, so including on-axis points, and return
#' the result as a matrix. For plotting purposes.
#' @inheritParams predictOffAxis
#' @inheritParams fitSurface
predictResponseSurface = function(doseGrid, fitResult, null_model,
                                  transforms = fitResult$transforms){
    fit = fitOffAxis(doseGrid, fitResult, null_model, startvalues = NULL)
    vec = (if(null_model %in% c("loewe")) fit$response else fit)
    names(vec) = getd1d2(doseGrid)
    if (!is.null(transforms)) {
        CompositeT <- with(transforms,
                           function(y, args) PowerT(BiolT(y, args), args))
        vec <- with(transforms, CompositeT(vec, compositeArgs))
    }
    out = matrix(0, length(unique(doseGrid$d1)), length(unique(doseGrid$d2)),
                 dimnames = list(sort(unique(doseGrid$d1)), sort(unique(doseGrid$d2))))
    for(n in names(vec)){
        foo = strsplit(n, split = "_")[[1]]
        out[foo[1], foo[2]] = vec[n]
    }
    out
}
