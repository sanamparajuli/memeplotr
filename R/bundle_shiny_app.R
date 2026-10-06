#' Build a self-contained bundle for shinyapps.io / Posit Connect
#'
#' `memeplotr` is not on CRAN, so a plain `library(memeplotr)` in an `app.R`
#' fails when `rsconnect` resolves dependencies against a repository. This
#' function writes a deployable Shiny directory that carries the package
#' *source* inside it and loads that source at startup with
#' [pkgload::load_all()], so the only packages the server installs are the CRAN
#' ones named in `DESCRIPTION`.
#'
#' The result is an ordinary Shiny app directory: check it with
#' [shiny::runApp()], then upload it with `rsconnect::deployApp()`.
#'
#' @param dir Directory to create.
#' @param pkg_source Package source to embed: an unpacked source directory
#'   (one containing `DESCRIPTION`, `R/` and `inst/`) or a source tarball.
#'   `NULL` (default) searches the working directory for a `memeplotr`
#'   directory, then for the newest `memeplotr_*.tar.gz`.
#' @param overwrite Replace `dir` if it already exists.
#' @param quiet Suppress the message naming what was written.
#'
#' @return The bundle path, invisibly.
#'
#' @details
#' Two substitutions are applied to the packaged `app.R`: `library(memeplotr)`
#' is dropped, since the package is loaded from source instead, and lookups of
#' the form `system.file(package = "memeplotr")` are redirected into the
#' embedded tree. Nothing else is rewritten, so the deployed app and
#' [run_memeplotr()] cannot drift apart.
#'
#' `pkgload` is a runtime dependency *of the bundle*, not of the package.
#'
#' @examples
#' \dontrun{
#' bundle_shiny_app("memeplotr-shinyapp")
#' shiny::runApp("memeplotr-shinyapp")          # check locally first
#' rsconnect::deployApp("memeplotr-shinyapp",
#'                      appName = "memeplotr")  # -> https://<acct>.shinyapps.io/memeplotr/
#' }
#' @export
bundle_shiny_app <- function(dir = "memeplotr-shinyapp", pkg_source = NULL,
                             overwrite = TRUE, quiet = FALSE) {
  pkg_source <- .resolve_pkg_source(pkg_source)
  on.exit(unlink(attr(pkg_source, "tmp"), recursive = TRUE), add = TRUE)

  if (dir.exists(dir)) {
    if (!isTRUE(overwrite)) {
      .abort(c(sprintf("'%s' already exists.", dir),
               i = "Pass `overwrite = TRUE` to replace it."))
    }
    unlink(dir, recursive = TRUE)
  }
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  ## Copy the source verbatim, minus the parts a running app never reads.
  ## vignettes/ and tests/ would otherwise inflate the upload and drag
  ## knitr/testthat into the server's dependency set.
  pkg <- file.path(dir, "memeplotr")
  dir.create(pkg, showWarnings = FALSE)
  keep <- c("DESCRIPTION", "NAMESPACE", "LICENSE", "R", "inst")
  keep <- keep[file.exists(file.path(pkg_source, keep))]
  file.copy(file.path(pkg_source, keep), pkg, recursive = TRUE)

  app_src <- file.path(pkg, "inst", "shiny", "app.R")
  if (!file.exists(app_src)) {
    .abort(c("The package source has no inst/shiny/app.R.",
             i = sprintf("Looked inside '%s'.", pkg_source)))
  }
  app <- readLines(app_src, warn = FALSE)
  app <- app[!grepl("^\\s*library\\(memeplotr\\)\\s*$", app)]
  writeLines(app, file.path(dir, "app_body.R"))
  writeLines(.bundle_launcher(), file.path(dir, "app.R"))

  if (!quiet) {
    message(sprintf(
      "Shiny bundle at '%s': app.R, app_body.R, memeplotr/ (%s). Deploy with rsconnect::deployApp(\"%s\").",
      dir, .desc_version(file.path(pkg, "DESCRIPTION")), dir))
  }
  invisible(dir)
}

## Accept a directory, a tarball, or NULL (search the working directory).
## Returns the path to an unpacked source tree, with any temporary extraction
## directory attached so the caller can clean it up.
.resolve_pkg_source <- function(x) {
  if (is.null(x)) {
    if (dir.exists("memeplotr") && file.exists(file.path("memeplotr", "DESCRIPTION"))) {
      x <- "memeplotr"
    } else {
      tb <- sort(list.files(".", pattern = "^memeplotr_.*\\.tar\\.gz$"),
                 decreasing = TRUE)
      if (!length(tb)) {
        .abort(c("Could not find the memeplotr source.",
                 i = paste("Pass `pkg_source =` pointing at the unpacked",
                           "source directory or a memeplotr_*.tar.gz."),
                 i = sprintf("Searched '%s'.", normalizePath("."))))
      }
      x <- tb[[1L]]
    }
  }
  if (dir.exists(x)) {
    if (!file.exists(file.path(x, "DESCRIPTION"))) {
      .abort(sprintf("'%s' is not a package source directory (no DESCRIPTION).", x))
    }
    return(structure(x, tmp = NULL))
  }
  if (!file.exists(x)) .abort(sprintf("'%s' does not exist.", x))
  tmp <- tempfile("memeplotr-src-")
  dir.create(tmp)
  utils::untar(x, exdir = tmp)
  unpacked <- file.path(tmp, "memeplotr")
  if (!file.exists(file.path(unpacked, "DESCRIPTION"))) {
    .abort(sprintf("'%s' does not contain a memeplotr/ source tree.", x))
  }
  structure(unpacked, tmp = tmp)
}

.desc_version <- function(path) {
  v <- tryCatch(read.dcf(path, "Version")[[1L]], error = function(e) NA_character_)
  if (is.na(v)) "version unknown" else paste("version", v)
}

## Launcher written into the bundle. Held as a character vector rather than an
## inst/ template so a bundle built from one version never picks up another's.
.bundle_launcher <- function() {
  c(
    "## Generated by memeplotr::bundle_shiny_app(). Do not edit by hand --",
    "## regenerate it instead, so the app stays in step with the package.",
    "",
    "if (!requireNamespace(\"pkgload\", quietly = TRUE)) {",
    "  stop(\"This Shiny bundle needs the 'pkgload' package.\", call. = FALSE)",
    "}",
    "pkgload::load_all(\"memeplotr\", export_all = FALSE, helpers = FALSE,",
    "                  attach_testthat = FALSE, warn_conflicts = FALSE,",
    "                  quiet = TRUE)",
    "",
    "## app_body.R asks the installed package for its bundled example files.",
    "## Nothing is installed here, so send those lookups into the bundle.",
    "system.file <- function(..., package = \"base\", lib.loc = NULL,",
    "                        mustWork = FALSE) {",
    "  if (identical(package, \"memeplotr\")) {",
    "    p <- file.path(\"memeplotr\", \"inst\", ...)",
    "    return(if (file.exists(p)) normalizePath(p) else \"\")",
    "  }",
    "  base::system.file(..., package = package, lib.loc = lib.loc,",
    "                    mustWork = mustWork)",
    "}",
    "",
    "## Named here only so rsconnect's dependency scan installs them; the app",
    "## itself reaches for each one conditionally.",
    "if (FALSE) {",
    "  library(bslib); library(DT); library(colourpicker)",
    "  library(shinyjs); library(svglite)",
    "}",
    "",
    "source(\"app_body.R\", local = TRUE)$value"
  )
}
