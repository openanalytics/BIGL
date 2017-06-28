#' Run the BIGL application for demonstrating response surfaces
#'
#' @param ... Pass parameters to \code{\link[shiny]{runApp}}
#' @export
#' @examples
#' \dontrun{
#'   runBIGL()
#' }
runBIGL <- function(...){

  if (requireNamespace("shiny", quietly = TRUE))
    shiny::runApp(appDir = system.file("ui", package = "BIGL"), ...)
  else
    stop("shiny package needs to be installed for the interactive application to work.")

}
