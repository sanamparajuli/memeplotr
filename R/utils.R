# ---- internal helpers -------------------------------------------------------

#' @keywords internal
#' @noRd
`%||%` <- rlang::`%||%`

#' @keywords internal
#' @noRd
.as_num <- function(x) suppressWarnings(as.numeric(x))

#' @keywords internal
#' @noRd
.as_int <- function(x) suppressWarnings(as.integer(x))

#' @keywords internal
#' @noRd
.blank_to_na <- function(x) {
  x[!nzchar(trimws(x %||% ""))] <- NA_character_
  x
}

#' @keywords internal
#' @noRd
.norm_strand <- function(x) {
  x <- tolower(as.character(x %||% NA))
  out <- rep("*", length(x))
  out[x %in% c("plus", "+", "forward", "fwd", "1")] <- "+"
  out[x %in% c("minus", "-", "reverse", "rev", "-1")] <- "-"
  out[is.na(x)] <- "*"
  out
}

#' @keywords internal
#' @noRd
.abort <- function(msg, ...) rlang::abort(c(msg, ...), call = NULL)

#' @keywords internal
#' @noRd
.need <- function(pkg, what) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    .abort(sprintf("Package '%s' is required for %s. Install it with install.packages('%s').",
                   pkg, what, pkg))
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.consensus_from_pwm <- function(pwm, threshold = 0) {
  if (is.null(pwm) || !nrow(pwm)) return(NA_character_)
  paste(colnames(pwm)[max.col(pwm, ties.method = "first")], collapse = "")
}

#' Shannon information content of PWM rows
#' @keywords internal
#' @noRd
.pwm_ic <- function(pwm, pseudo = 1e-9) {
  p <- pmax(as.matrix(pwm), pseudo)
  p <- p / rowSums(p)
  h <- -rowSums(p * log2(p))
  log2(ncol(p)) - h
}

## millimetres (geom_text `size`) to points (theme `size`); ggplot2's .pt
.pt_per_mm <- 72.27 / 25.4
