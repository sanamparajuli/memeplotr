# ---- block shapes -----------------------------------------------------------

#' Available motif block shapes
#'
#' @return A character vector of shape names accepted by the `shape` argument
#'   of [geom_motif()] and [gg_motif_map()].
#' @examples
#' motif_shapes()
#' @export
motif_shapes <- function() {
  c("rect", "arrow", "chevron", "hexagon", "roundrect", "capsule", "tag", "notch")
}

#' @keywords internal
#' @noRd
.arc_xy <- function(cx, cy, rx, ry, a0, a1, n = 12) {
  a <- seq(a0, a1, length.out = n)
  list(x = cx + rx * cos(a), y = cy + ry * sin(a))
}

#' Build the polygon outline for a single motif block
#'
#' @param x0,x1 Left and right edge in data units.
#' @param y0,y1 Bottom and top edge in data units.
#' @param shape One of [motif_shapes()].
#' @param dir `1` for a right-pointing block, `-1` for left-pointing.
#' @param head_frac Fraction of the block length used by the point/head.
#' @param round_frac Corner rounding, as a fraction of half the block height.
#' @param rx_per_ry Data-unit aspect correction (x units per y unit) so that
#'   rounded corners look circular on the rendered panel.
#' @param n_arc Points per arc.
#'
#' @return A two-column data frame of polygon vertices.
#' @keywords internal
#' @noRd
.motif_poly <- function(x0, x1, y0, y1, shape = "rect", dir = 1,
                        head_frac = 0.3, round_frac = 0.35,
                        rx_per_ry = 1, n_arc = 10) {
  w <- x1 - x0
  ym <- (y0 + y1) / 2
  hh <- (y1 - y0) / 2
  h <- min(head_frac * w, w * 0.9)
  if (dir < 0) {
    p <- .motif_poly(x0, x1, y0, y1, shape, dir = 1, head_frac, round_frac,
                     rx_per_ry, n_arc)
    p$x <- x0 + x1 - p$x
    return(p[rev(seq_len(nrow(p))), , drop = FALSE])
  }
  out <- switch(
    shape,
    rect = list(x = c(x0, x1, x1, x0), y = c(y0, y0, y1, y1)),
    arrow = list(x = c(x0, x1 - h, x1, x1 - h, x0),
                 y = c(y0, y0, ym, y1, y1)),
    chevron = list(x = c(x0, x1 - h, x1, x1 - h, x0, x0 + h),
                   y = c(y0, y0, ym, y1, y1, ym)),
    notch = list(x = c(x0, x1, x1, x0, x0 + h),
                 y = c(y0, y0, y1, y1, ym)),
    hexagon = list(x = c(x0, x0 + h, x1 - h, x1, x1 - h, x0 + h),
                   y = c(ym, y0, y0, ym, y1, y1)),
    tag = list(x = c(x0, x1 - h, x1, x1, x1 - h, x0),
               y = c(y0, y0, ym - hh * 0.35, ym + hh * 0.35, y1, y1)),
    roundrect = {
      ry <- round_frac * hh
      rx <- min(ry * rx_per_ry, w / 2)
      a1 <- .arc_xy(x1 - rx, y0 + ry, rx, ry, -pi / 2, 0, n_arc)
      a2 <- .arc_xy(x1 - rx, y1 - ry, rx, ry, 0, pi / 2, n_arc)
      a3 <- .arc_xy(x0 + rx, y1 - ry, rx, ry, pi / 2, pi, n_arc)
      a4 <- .arc_xy(x0 + rx, y0 + ry, rx, ry, pi, 3 * pi / 2, n_arc)
      list(x = c(a1$x, a2$x, a3$x, a4$x), y = c(a1$y, a2$y, a3$y, a4$y))
    },
    capsule = {
      rx <- min(hh * rx_per_ry, w / 2)
      a1 <- .arc_xy(x1 - rx, ym, rx, hh, -pi / 2, pi / 2, n_arc * 2)
      a2 <- .arc_xy(x0 + rx, ym, rx, hh, pi / 2, 3 * pi / 2, n_arc * 2)
      list(x = c(a1$x, a2$x), y = c(a1$y, a2$y))
    },
    .abort(sprintf("Unknown shape '%s'. See motif_shapes().", shape))
  )
  data.frame(x = out$x, y = out$y)
}

# ---- Stat -------------------------------------------------------------------

#' @rdname geom_motif
#' @format NULL
#' @usage NULL
#' @export
StatMotif <- ggplot2::ggproto(
  "StatMotif", ggplot2::Stat,
  required_aes = c("xmin", "xmax", "y"),
  compute_panel = function(data, scales, shape = "rect", height = 0.6,
                           head_frac = 0.3, round_frac = 0.35, n_arc = 10,
                           na.rm = FALSE) {
    if (!nrow(data)) return(data)
    n <- nrow(data)
    shp <- if (!is.null(data$shape)) as.character(data$shape) else rep_len(shape, n)
    hgt <- if (!is.null(data$height)) data$height else rep_len(height, n)
    dir <- rep(1, n)
    if (!is.null(data$strand)) {
      s <- as.character(data$strand)
      dir[s %in% c("-", "minus", "reverse")] <- -1
    }
    ## aspect correction so rounded ends are not stretched by the x range
    rx_per_ry <- 1
    xr <- tryCatch(scales$x$dimension(), error = function(e) NULL)
    yr <- tryCatch(scales$y$dimension(), error = function(e) NULL)
    if (!is.null(xr) && !is.null(yr) && diff(yr) > 0) {
      rx_per_ry <- diff(xr) / diff(yr) * 0.35
    }
    keep <- setdiff(names(data), c("x", "y"))
    pieces <- lapply(seq_len(n), function(i) {
      p <- .motif_poly(data$xmin[i], data$xmax[i],
                       data$y[i] - hgt[i] / 2, data$y[i] + hgt[i] / 2,
                       shape = shp[i], dir = dir[i], head_frac = head_frac,
                       round_frac = round_frac, rx_per_ry = rx_per_ry,
                       n_arc = n_arc)
      d <- data[rep(i, nrow(p)), keep, drop = FALSE]
      d$x <- p$x
      d$y <- p$y
      d$group <- i
      d
    })
    out <- do.call(rbind, pieces)
    rownames(out) <- NULL
    out
  }
)

#' Draw motif occurrences as shaped blocks
#'
#' A `ggplot2` layer that renders each motif hit as a polygon whose silhouette
#' encodes orientation. Unlike a plain [ggplot2::geom_rect()], the block can be
#' an arrow, chevron, hexagon, rounded rectangle or capsule, and the direction
#' follows a `strand` aesthetic when one is supplied.
#'
#' @param mapping,data,position,na.rm,show.legend,inherit.aes Standard
#'   `ggplot2` layer arguments.
#' @param stat The statistical transformation; defaults to `"motif"`, which
#'   converts `xmin`/`xmax`/`y` into polygon vertices.
#' @param geom For `stat_motif()`, the geom used to draw the computed
#'   vertices; `"polygon"` unless you have a reason to change it.
#' @param shape Block silhouette, one of [motif_shapes()]. Can also be supplied
#'   per row as a `shape` aesthetic.
#' @param height Block height in `y` data units (the axis is one unit per
#'   sequence, so `0.6` leaves a 40\% gap).
#' @param head_frac Fraction of block length used by the arrow head / point.
#' @param round_frac Corner radius for `shape = "roundrect"`, as a fraction of
#'   half the block height.
#' @param n_arc Number of vertices per arc for curved shapes.
#' @param ... Other arguments passed to the layer, e.g. `colour`, `linewidth`,
#'   `alpha`.
#'
#' @section Aesthetics:
#' Required: `xmin`, `xmax`, `y`. Understood: `fill`, `colour`, `linewidth`,
#' `alpha`, `strand`, `shape`, `height`, `group`.
#'
#' @return A `ggplot2` layer.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' mt <- as_motif_table(res)
#' mt$y <- as.integer(factor(mt$sequence))
#'
#' ggplot2::ggplot(mt) +
#'   geom_motif(ggplot2::aes(xmin = start, xmax = stop, y = y, fill = motif),
#'              shape = "chevron") +
#'   theme_motif()
#' @export
geom_motif <- function(mapping = NULL, data = NULL, stat = "motif",
                       position = "identity", ..., shape = "arrow",
                       height = 0.6, head_frac = 0.3, round_frac = 0.35,
                       n_arc = 10, na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = stat, geom = ggplot2::GeomPolygon,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(shape = shape, height = height, head_frac = head_frac,
                  round_frac = round_frac, n_arc = n_arc, na.rm = na.rm, ...)
  )
}

#' @rdname geom_motif
#' @export
stat_motif <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", ..., shape = "arrow",
                       height = 0.6, head_frac = 0.3, round_frac = 0.35,
                       n_arc = 10, na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE) {
  ggplot2::layer(
    data = data, mapping = mapping, stat = StatMotif, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(shape = shape, height = height, head_frac = head_frac,
                  round_frac = round_frac, n_arc = n_arc, na.rm = na.rm, ...)
  )
}
