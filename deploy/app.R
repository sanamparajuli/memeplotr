## Deployment entry point for shinyapps.io / Posit Connect.
##
## This directory is what you hand to rsconnect::deployApp(). It is three
## lines because the app itself lives in the installed package; rsconnect
## works out that it needs memeplotr from the library(memeplotr) call below.
##
## memeplotr must have been installed with remotes::install_github(), not
## install_local() -- that is what writes RemoteType: github into its
## DESCRIPTION and lets the server install the same package. Check with:
##
##   packageDescription("memeplotr")$RemoteType   # "github"

library(shiny)
library(memeplotr)

## Declared so rsconnect's dependency scan installs them; the app reaches for
## each one conditionally, so it still runs if any is missing.
if (FALSE) {
  library(bslib); library(DT); library(colourpicker)
  library(shinyjs); library(svglite)
}

source(system.file("shiny", "app.R", package = "memeplotr"), local = TRUE)$value
