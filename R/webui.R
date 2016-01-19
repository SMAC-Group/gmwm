#' @title Launch Web Modeling UI
#' @description
#' Creates a Web UI Modeling Environment to work with the GMWM objects
#' @examples 
#' \dontrun{
#' # Launch the modeling UI
#' webui()
#' }
#' @export
webui = function(){
  appDir = system.file("webui", package = "gmwm")
  if(!nzchar(appDir)){
    stop("Could not find the webui directory. Try re-installing `gmwm`.", call. = FALSE)
  }
  
  if (requireNamespace("shiny", quietly = TRUE)) {
    message('Hit <escape> in console to stop the application.')
    shiny::runApp(appDir, display.mode = "normal")
  } else {
    stop("In order to use the `webui`, please install the `shiny` r package via: `install.packages('shiny')`")
  }
}