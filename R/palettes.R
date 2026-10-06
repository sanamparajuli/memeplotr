# ---- palettes ---------------------------------------------------------------

.memeplotr_palettes <- list(
  okabe_ito = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#D55E00", "#CC79A7", "#999999", "#000000"),
  bright    = c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
                "#AA3377", "#BBBBBB"),
  muted     = c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
                "#DDCC77", "#CC6677", "#882255", "#AA4499"),
  pastel    = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6",
                "#FFFF99", "#B15928", "#1F78B4"),
  set3      = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
                "#CCEBC5", "#FFED6F"),
  earth     = c("#A16928", "#BD925A", "#D6BD8D", "#EDEAC2", "#B5C8B8",
                "#79A7AC", "#2887A1"),
  grey      = c("#1A1A1A", "#4D4D4D", "#808080", "#B3B3B3", "#D9D9D9")
)

#' Colour palettes for motifs
#'
#' Returns `n` colours from one of the built-in qualitative palettes, or from a
#' sequential/continuous palette provided by [grDevices::hcl.colors()]. Palettes
#' shorter than `n` are interpolated rather than recycled, so no two motifs ever
#' receive the same colour by accident.
#'
#' @param n Number of colours required.
#' @param palette Palette name. Built-in qualitative options are
#'   `"okabe_ito"` (colourblind-safe, the default), `"bright"`, `"muted"`,
#'   `"pastel"`, `"set3"`, `"earth"` and `"grey"`. Any name accepted by
#'   [grDevices::hcl.colors()] also works, e.g. `"viridis"`, `"Zissou 1"`,
#'   `"Spectral"`, `"Batlow"`.
#' @param rev Reverse the palette order.
#'
#' @return A character vector of `n` hex colours.
#' @examples
#' motif_palette(6)
#' motif_palette(6, "viridis")
#' memeplotr_palettes()
#' @export
motif_palette <- function(n, palette = "okabe_ito", rev = FALSE) {
  ## names are matched case-insensitively: hcl.pals() capitalises its palettes
  ## ("Viridis", "Spectral") and users reasonably type them in lower case
  hcl <- grDevices::hcl.pals()
  if (length(palette) > 1L) {
    cols <- palette
  } else if (tolower(palette) %in% tolower(names(.memeplotr_palettes))) {
    cols <- .memeplotr_palettes[[match(tolower(palette),
                                       tolower(names(.memeplotr_palettes)))]]
  } else if (tolower(palette) %in% tolower(hcl)) {
    cols <- grDevices::hcl.colors(max(n, 3L), hcl[match(tolower(palette), tolower(hcl))])
  } else {
    .abort(sprintf("Unknown palette '%s'.", palette),
           i = paste0("Built in: ", paste(names(.memeplotr_palettes), collapse = ", "),
                      ". See grDevices::hcl.pals() for the rest."))
  }
  if (rev) cols <- base::rev(cols)
  if (n <= length(cols)) cols[seq_len(n)] else grDevices::colorRampPalette(cols)(n)
}

#' @rdname motif_palette
#' @return `memeplotr_palettes()` returns the names of the built-in palettes.
#' @export
memeplotr_palettes <- function() names(.memeplotr_palettes)

#' @keywords internal
#' @noRd
.discrete_scale <- function(aesthetics, palette, ...) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::discrete_scale(aesthetics = aesthetics, palette = palette, ...)
  } else {
    ggplot2::discrete_scale(aesthetics = aesthetics, scale_name = "motif",
                            palette = palette, ...)
  }
}

#' Motif colour and fill scales
#'
#' Discrete `ggplot2` scales backed by [motif_palette()]. Use `values` to pin
#' specific motifs to specific colours while the rest are taken from the
#' palette.
#'
#' @param palette Palette name, see [motif_palette()].
#' @param values Optional named character vector of colours, names matching
#'   motif labels. Motifs not named here fall back to the palette.
#' @param colours,colors Aliases for `values`, accepted so the scales take the
#'   same argument name as the `gg_motif_*()` plotting functions. The first of
#'   `values`, `colours`, `colors` that is supplied wins.
#' @param rev Reverse the palette order.
#' @param name Legend title.
#' @param ... Passed to [ggplot2::discrete_scale()] / [ggplot2::scale_fill_manual()].
#'
#' @return A `ggplot2` scale, addable with `+`.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_map(res) + scale_fill_motif("Spectral")
#'
#' # pin one motif, leave the rest to the palette
#' gg_motif_map(res) + scale_fill_motif(palette = "grey", values = c(`1` = "firebrick"))
#' gg_motif_map(res) + scale_fill_motif(values = c("MEME-1" = "black"))
#' @export
scale_fill_motif <- function(palette = "okabe_ito", values = NULL, rev = FALSE,
                             name = "Motif", colours = NULL, colors = NULL, ...) {
  .motif_scale("fill", palette, .first_non_null(values, colours, colors),
               rev, name, ...)
}

#' @rdname scale_fill_motif
#' @export
scale_colour_motif <- function(palette = "okabe_ito", values = NULL, rev = FALSE,
                               name = "Motif", colours = NULL, colors = NULL, ...) {
  .motif_scale("colour", palette, .first_non_null(values, colours, colors),
               rev, name, ...)
}

#' @rdname scale_fill_motif
#' @export
scale_color_motif <- scale_colour_motif

#' @keywords internal
#' @noRd
.motif_scale <- function(aes, palette, values, rev, name, ...) {
  if (is.null(values)) {
    return(.discrete_scale(aes, function(n) motif_palette(n, palette, rev),
                           name = name, ...))
  }
  structure(list(values = values, pal = palette, rev = rev, name = name,
                 aes = aes, dots = list(...)),
            class = c("motif_partial_scale", "gg"))
}

#' @export
ggplot_add.motif_partial_scale <- function(object, plot, ...) {
  lv <- .scale_levels(object$aes, plot)
  if (is.null(lv)) lv <- names(object$values)
  unmatched <- setdiff(names(object$values), lv)
  if (length(unmatched)) {
    rlang::warn(paste0(
      "scale_", object$aes, "_motif(): no motif named ",
      paste0("\"", unmatched, "\"", collapse = ", "),
      " in the data. Present: ", paste(lv, collapse = ", "), "."
    ))
  }
  lv <- unique(c(lv, names(object$values)))
  cols <- motif_palette(length(lv), object$pal, object$rev)
  names(cols) <- lv
  cols[names(object$values)] <- unname(object$values)
  manual <- if (identical(object$aes, "fill")) ggplot2::scale_fill_manual else ggplot2::scale_colour_manual
  args <- c(list(values = cols, name = object$name), object$dots)
  ggplot2::ggplot_add(do.call(manual, args), plot, ...)
}

# Discrete levels a scale will see, resolved before the plot is built: evaluate
# each layer's mapping for `aes` against that layer's own data.
#' @keywords internal
#' @noRd
.scale_levels <- function(aes, plot) {
  keys <- if (aes %in% c("colour", "color")) c("colour", "color") else aes
  lv <- NULL
  for (l in plot$layers) {
    ld <- if (inherits(l$data, "waiver")) plot$data else l$data
    if (!is.data.frame(ld) || !nrow(ld)) next
    mp <- l$mapping
    if (isTRUE(l$inherit.aes)) mp <- utils::modifyList(as.list(plot$mapping), as.list(mp))
    q <- NULL
    for (k in keys) if (!is.null(mp[[k]])) { q <- mp[[k]]; break }
    if (is.null(q)) next
    v <- tryCatch(rlang::eval_tidy(q, ld), error = function(e) NULL)
    if (is.null(v) || is.numeric(v)) next
    lv <- union(lv, if (is.factor(v)) levels(droplevels(v)) else sort(unique(as.character(v))))
  }
  if (is.null(lv) && is.data.frame(plot$data) && "motif" %in% names(plot$data)) {
    lv <- levels(factor(plot$data$motif))
  }
  lv
}

#' @keywords internal
#' @noRd
.resolve_motif_colors <- function(levels, palette = "okabe_ito", colors = NULL, rev = FALSE) {
  cols <- motif_palette(length(levels), palette, rev)
  names(cols) <- levels
  if (!is.null(colors)) {
    if (is.null(names(colors))) {
      k <- min(length(colors), length(levels))
      cols[seq_len(k)] <- colors[seq_len(k)]
    } else {
      keep <- intersect(names(colors), levels)
      cols[keep] <- colors[keep]
    }
  }
  cols
}

#' @keywords internal
#' @noRd
.first_non_null <- function(...) {
  for (x in list(...)) if (!is.null(x)) return(x)
  NULL
}
