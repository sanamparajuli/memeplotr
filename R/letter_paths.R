# ---- vector letter glyphs ---------------------------------------------------
#
# Each glyph is a list of polygons defined on the unit square [0, 1] x [0, 1].
# Building them from primitives (strokes and elliptical ring segments) keeps the
# package font-free and device independent: logo letters are plain ggplot2
# polygons, so they obey fill, colour, alpha and every other aesthetic.

#' @keywords internal
#' @noRd
.gl_rect <- function(x0, y0, x1, y1) {
  data.frame(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1))
}

#' @keywords internal
#' @noRd
.gl_seg <- function(x0, y0, x1, y1, w) {
  dx <- x1 - x0; dy <- y1 - y0
  len <- sqrt(dx^2 + dy^2)
  if (len < 1e-9) return(.gl_rect(x0, y0, x0 + w, y0 + w))
  nx <- -dy / len * w / 2; ny <- dx / len * w / 2
  data.frame(x = c(x0 + nx, x1 + nx, x1 - nx, x0 - nx),
             y = c(y0 + ny, y1 + ny, y1 - ny, y0 - ny))
}

#' @keywords internal
#' @noRd
.gl_ring <- function(cx, cy, rx, ry, w, a0 = 0, a1 = 2 * pi, n = 36) {
  a <- seq(a0, a1, length.out = n)
  rxi <- max(rx - w, 1e-3); ryi <- max(ry - w, 1e-3)
  data.frame(x = c(cx + rx * cos(a), base::rev(cx + rxi * cos(a))),
             y = c(cy + ry * sin(a), base::rev(cy + ryi * sin(a))))
}

#' @keywords internal
#' @noRd
.glyph_defs <- function(w = 0.19) {
  list(
    A = list(.gl_seg(0.06, 0, 0.50, 1, w), .gl_seg(0.94, 0, 0.50, 1, w),
             .gl_rect(0.21, 0.30, 0.79, 0.30 + w * 0.92)),
    B = list(.gl_rect(0, 0, w, 1),
             .gl_ring(w, 0.745, 0.68, 0.255, w * 0.92, -pi / 2, pi / 2),
             .gl_ring(w, 0.255, 0.74, 0.255, w * 0.92, -pi / 2, pi / 2),
             .gl_rect(0, 0.5 - w / 2, 0.55, 0.5 + w / 2)),
    C = list(.gl_ring(0.50, 0.50, 0.50, 0.50, w, 0.32 * pi, 1.68 * pi)),
    D = list(.gl_rect(0, 0, w, 1),
             .gl_ring(w, 0.50, 1 - w, 0.50, w, -pi / 2, pi / 2)),
    E = list(.gl_rect(0, 0, w, 1), .gl_rect(0, 0, 0.95, w),
             .gl_rect(0, 0.5 - w / 2, 0.82, 0.5 + w / 2),
             .gl_rect(0, 1 - w, 0.95, 1)),
    F = list(.gl_rect(0, 0, w, 1),
             .gl_rect(0, 0.5 - w / 2, 0.80, 0.5 + w / 2),
             .gl_rect(0, 1 - w, 0.95, 1)),
    ## the bar attaches at the arc's lower-right end and stops short of the
    ## free upper-right end, leaving the opening that distinguishes G from O
    G = list(.gl_ring(0.50, 0.50, 0.50, 0.50, w, 0.19 * pi, 1.89 * pi),
             .gl_rect(0.55, 0.30, 1, 0.30 + w * 0.92)),
    H = list(.gl_rect(0, 0, w, 1), .gl_rect(1 - w, 0, 1, 1),
             .gl_rect(0, 0.5 - w / 2, 1, 0.5 + w / 2)),
    I = list(.gl_rect(0.5 - w / 2, 0, 0.5 + w / 2, 1),
             .gl_rect(0.14, 0, 0.86, w), .gl_rect(0.14, 1 - w, 0.86, 1)),
    J = list(.gl_rect(0.62, 0.30, 0.81, 1),
             .gl_ring(0.405, 0.30, 0.405, 0.30, w, pi, 2 * pi)),
    K = list(.gl_rect(0, 0, w, 1), .gl_seg(0.12, 0.47, 1, 1, w),
             .gl_seg(0.12, 0.47, 1, 0, w)),
    L = list(.gl_rect(0, 0, w, 1), .gl_rect(0, 0, 0.95, w)),
    M = list(.gl_rect(0, 0, w, 1), .gl_rect(1 - w, 0, 1, 1),
             .gl_seg(0.09, 1, 0.50, 0.26, w), .gl_seg(0.91, 1, 0.50, 0.26, w)),
    N = list(.gl_rect(0, 0, w, 1), .gl_rect(1 - w, 0, 1, 1),
             .gl_seg(0.09, 1, 0.91, 0, w)),
    O = list(.gl_ring(0.50, 0.50, 0.50, 0.50, w)),
    P = list(.gl_rect(0, 0, w, 1),
             .gl_ring(w, 0.745, 1 - w, 0.255, w, -pi / 2, pi / 2),
             .gl_rect(0, 0.49, 0.55, 0.49 + w)),
    Q = list(.gl_ring(0.50, 0.50, 0.50, 0.50, w),
             .gl_seg(0.58, 0.30, 0.98, 0.02, w)),
    R = list(.gl_rect(0, 0, w, 1),
             .gl_ring(w, 0.745, 1 - w, 0.255, w, -pi / 2, pi / 2),
             .gl_rect(0, 0.49, 0.55, 0.49 + w),
             .gl_seg(0.33, 0.52, 1, 0, w)),
    ## the two bowls overlap at the waist: the upper is open bottom-right,
    ## the lower open top-left, which is what makes the stroke continuous
    S = list(.gl_ring(0.50, 0.715, 0.46, 0.285, w, 0, 1.5 * pi),
             .gl_ring(0.50, 0.285, 0.46, 0.285, w, 0.5 * pi, -pi)),
    T = list(.gl_rect(0.5 - w / 2, 0, 0.5 + w / 2, 1), .gl_rect(0, 1 - w, 1, 1)),
    U = list(.gl_ring(0.50, 0.42, 0.50, 0.42, w, pi, 2 * pi),
             .gl_rect(0, 0.42, w, 1), .gl_rect(1 - w, 0.42, 1, 1)),
    V = list(.gl_seg(0.05, 1, 0.50, 0.02, w), .gl_seg(0.95, 1, 0.50, 0.02, w)),
    W = list(.gl_seg(0.02, 1, 0.26, 0.02, w), .gl_seg(0.26, 0.02, 0.50, 0.70, w),
             .gl_seg(0.50, 0.70, 0.74, 0.02, w), .gl_seg(0.74, 0.02, 0.98, 1, w)),
    X = list(.gl_seg(0.05, 0, 0.95, 1, w), .gl_seg(0.05, 1, 0.95, 0, w)),
    Y = list(.gl_seg(0.05, 1, 0.50, 0.46, w), .gl_seg(0.95, 1, 0.50, 0.46, w),
             .gl_rect(0.5 - w / 2, 0, 0.5 + w / 2, 0.52)),
    Z = list(.gl_rect(0, 1 - w, 1, 1), .gl_rect(0, 0, 1, w),
             .gl_seg(0.07, w, 0.93, 1 - w, w)),
    `-` = list(.gl_rect(0.08, 0.44, 0.92, 0.56)),
    `.` = list(.gl_rect(0.38, 0, 0.62, 0.22))
  )
}

#' Vector outlines for logo letters
#'
#' Returns the polygon outline of a character on the unit square. Used
#' internally by [gg_motif_logo()]; exported so that you can draw letters in
#' your own layers or swap in your own glyph set.
#'
#' @param letter A single character. Lower case is upper-cased; unknown
#'   characters fall back to a hollow box.
#' @param stroke Stroke thickness, as a fraction of the glyph box.
#' @param normalise Rescale the finished glyph so that it exactly fills the
#'   unit square. `TRUE` (the default) guarantees a letter stays inside its
#'   cell; `FALSE` keeps the raw construction coordinates.
#'
#' @return A data frame with columns `x`, `y`, `piece` - one row per polygon
#'   vertex, `piece` separating the strokes that make up the glyph.
#' @examples
#' head(letter_outline("W"))
#' @export
letter_outline <- function(letter, stroke = 0.19, normalise = TRUE) {
  letter <- toupper(as.character(letter)[1])
  defs <- .glyph_defs(stroke)
  g <- defs[[letter]]
  if (is.null(g)) {
    g <- list(.gl_rect(0, 0, 1, stroke), .gl_rect(0, 1 - stroke, 1, 1),
              .gl_rect(0, 0, stroke, 1), .gl_rect(1 - stroke, 0, 1, 1))
  }
  out <- dplyr::bind_rows(lapply(seq_along(g), function(i) {
    d <- g[[i]]; d$piece <- i; d
  }))
  ## offsetting the diagonal strokes pushes a few glyphs (K, R, V, W) a few
  ## percent outside the unit square; rescaling to the realised bounding box
  ## keeps every letter exactly inside its cell, so stacked letters never
  ## bleed into their neighbours.
  if (isTRUE(normalise)) {
    rx <- range(out$x); ry <- range(out$y)
    if (diff(rx) > 0) out$x <- (out$x - rx[1L]) / diff(rx)
    if (diff(ry) > 0) out$y <- (out$y - ry[1L]) / diff(ry)
  }
  out
}

#' @keywords internal
#' @noRd
.place_letter <- function(letter, x0, x1, y0, y1, stroke = 0.19) {
  o <- letter_outline(letter, stroke)
  o$x <- x0 + o$x * (x1 - x0)
  o$y <- y0 + o$y * (y1 - y0)
  o
}
