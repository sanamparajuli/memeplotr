#' Launch the memeplotr Shiny application
#'
#' An interactive front end for the package: load a MEME XML file and an
#' optional Newick tree, build the figure by clicking, then copy the generated
#' `ggplot2` call from the "R code" tab into your script. Every control maps
#' onto an argument of the plotting functions, so nothing you can do in the app
#' is unreachable from the command line.
#'
#' @param launch.browser Open the app in a browser. Passed to
#'   [shiny::runApp()].
#' @param port Port to serve on; `NULL` lets Shiny choose.
#' @param ... Further arguments for [shiny::runApp()].
#'
#' @return Invisibly, the result of [shiny::runApp()]. Called for its side
#'   effect of starting the application.
#' @examples
#' if (interactive()) {
#'   run_memeplotr()
#' }
#' @export
run_memeplotr <- function(launch.browser = TRUE, port = NULL, ...) {
  for (pkg in c("shiny", "bslib", "DT")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      .abort(c(sprintf("Package '%s' is needed to run the app.", pkg),
               i = sprintf('install.packages("%s")', pkg)))
    }
  }
  app_dir <- system.file("shiny", package = "memeplotr")
  if (!nzchar(app_dir)) .abort("Could not find the app directory; reinstall memeplotr.")
  shiny::runApp(app_dir, launch.browser = launch.browser, port = port, ...)
}
