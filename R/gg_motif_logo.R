# ---- letter colour schemes --------------------------------------------------

.aa_chemistry <- c(
  G = "polar", S = "polar", T = "polar", Y = "polar", C = "polar",
  Q = "neutral", N = "neutral",
  K = "basic", R = "basic", H = "basic",
  D = "acidic", E = "acidic",
  A = "hydrophobic", V = "hydrophobic", L = "hydrophobic", I = "hydrophobic",
  P = "hydrophobic", W = "hydrophobic", F = "hydrophobic", M = "hydrophobic"
)
.chem_cols <- c(polar = "#109648", neutral = "#5E239D", basic = "#255C99",
                acidic = "#D62839", hydrophobic = "#1A1A1A")

## Kyte-Doolittle hydropathy
.aa_kd <- c(A = 1.8, R = -4.5, N = -3.5, D = -3.5, C = 2.5, Q = -3.5, E = -3.5,
            G = -0.4, H = -3.2, I = 4.5, L = 3.8, K = -3.9, M = 1.9, F = 2.8,
            P = -1.6, S = -0.8, T = -0.7, W = -0.9, Y = -1.3, V = 4.2)

.aa_taylor <- c(
  A = "#CCFF00", C = "#FFFF00", D = "#FF0000", E = "#FF0066", F = "#00FF66",
  G = "#FF9900", H = "#0066FF", I = "#66FF00", K = "#6600FF", L = "#33FF00",
  M = "#00FF00", N = "#CC00FF", P = "#FFCC00", Q = "#FF00CC", R = "#0000FF",
  S = "#FF3300", T = "#FF6600", V = "#99FF00", W = "#00CCFF", Y = "#00FFCC"
)

.nt_cols <- c(A = "#3DA64A", C = "#2B6CB0", G = "#DD9B16", T = "#C7362F",
              U = "#C7362F", N = "#9E9E9E")

#' Letter colour schemes for sequence logos
#'
#' @param scheme One of `"auto"`, `"chemistry"`, `"hydrophobicity"`,
#'   `"taylor"`, `"nucleotide"`, `"mono"`.
#' @param letters The alphabet in use, used by `"auto"` to pick between the
#'   nucleotide and protein schemes and to complete the returned vector.
#' @param mono_colour Colour used by `scheme = "mono"`.
#'
#' @return A named character vector of colours, one per letter.
#' @examples
#' logo_colours("chemistry", LETTERS[1:5])
#' logo_colours("nucleotide", c("A", "C", "G", "T"))
#' @export
logo_colours <- function(scheme = c("auto", "chemistry", "hydrophobicity",
                                    "taylor", "nucleotide", "mono"),
                         letters = NULL, mono_colour = "grey20") {
  scheme <- match.arg(scheme)
  letters <- toupper(letters %||% names(.aa_chemistry))
  if (scheme == "auto") {
    scheme <- if (all(letters %in% c("A", "C", "G", "T", "U", "N"))) "nucleotide" else "chemistry"
  }
  out <- switch(
    scheme,
    nucleotide = .nt_cols,
    taylor = .aa_taylor,
    mono = stats::setNames(rep(mono_colour, length(letters)), letters),
    chemistry = stats::setNames(unname(.chem_cols[.aa_chemistry]), names(.aa_chemistry)),
    hydrophobicity = {
      v <- (.aa_kd - min(.aa_kd)) / diff(range(.aa_kd))
      ramp <- grDevices::colorRamp(c("#2166AC", "#F7F7F7", "#B2182B"))
      stats::setNames(grDevices::rgb(ramp(v), maxColorValue = 255), names(.aa_kd))
    }
  )
  miss <- setdiff(letters, names(out))
  if (length(miss)) out <- c(out, stats::setNames(rep("grey60", length(miss)), miss))
  out[letters]
}

#' @rdname logo_colours
#' @export
logo_colors <- logo_colours

# ---- logo data --------------------------------------------------------------

#' Stacked-letter heights for a motif
#'
#' Turns a position probability matrix into the rectangle each letter occupies
#' in a sequence logo. Useful if you want to build the logo layer yourself.
#'
#' @param x A [meme_result], a position probability matrix (positions in rows,
#'   letters in columns), or a list of such matrices.
#' @param motif Motif label or index; required when `x` is a [meme_result] with
#'   more than one motif.
#' @param method `"bits"` for information content (the MEME/WebLogo default) or
#'   `"probability"` for a stack of total height 1.
#' @param small_sample_correction Apply the Schneider small-sample correction
#'   using the motif's site count (`method = "bits"` only).
#' @param min_height Drop letters whose height is below this fraction of the
#'   stack, keeping the figure legible.
#'
#' @return A tibble with `position`, `letter`, `height`, `ymin`, `ymax`.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' head(logo_heights(res, "MEME-3"))
#' @export
logo_heights <- function(x, motif = NULL, method = c("bits", "probability"),
                         small_sample_correction = FALSE, min_height = 0.002) {
  method <- match.arg(method)
  pwm <- .get_pwm(x, motif)
  nsites <- attr(pwm, "nsites") %||% NA_real_
  p <- as.matrix(pwm)
  p[p < 0] <- 0
  rs <- rowSums(p); rs[rs == 0] <- 1
  p <- p / rs
  total <- if (method == "bits") {
    ic <- .pwm_ic(p)
    if (isTRUE(small_sample_correction) && is.finite(nsites) && nsites > 0) {
      ic <- ic - (ncol(p) - 1) / (2 * log(2) * nsites)
    }
    pmax(ic, 0)
  } else {
    rep(1, nrow(p))
  }
  out <- dplyr::bind_rows(lapply(seq_len(nrow(p)), function(i) {
    h <- p[i, ] * total[i]
    keep <- h > min_height * max(total[i], 1e-9)
    if (!any(keep)) return(NULL)
    h <- sort(h[keep])
    tibble::tibble(position = i, letter = names(h), height = unname(h),
                   ymin = unname(c(0, cumsum(h)[-length(h)])),
                   ymax = unname(cumsum(h)))
  }))
  attr(out, "stack_max") <- if (method == "bits") log2(ncol(p)) else 1
  attr(out, "method") <- method
  out
}

#' @keywords internal
#' @noRd
.get_pwm <- function(x, motif = NULL) {
  if (inherits(x, "meme_result")) {
    nm <- names(x$pwm)
    idx <- if (is.null(motif)) {
      if (length(nm) > 1L) .abort("`motif` is required: this result has several motifs.",
                                  i = paste("Available:", paste(nm, collapse = ", ")))
      1L
    } else if (is.numeric(motif)) {
      as.integer(motif)
    } else {
      ## accept the display label, the raw id, or the alternative name
      hit <- match(as.character(motif), x$motifs$motif)
      if (is.na(hit)) hit <- match(as.character(motif), x$motifs$motif_id)
      if (is.na(hit)) hit <- match(as.character(motif), x$motifs$name)
      if (is.na(hit)) hit <- match(as.character(motif), nm)
      if (is.na(hit)) {
        .abort(sprintf("Motif '%s' not found.", motif),
               i = paste("Available:", paste(x$motifs$motif, collapse = ", ")))
      }
      hit
    }
    m <- x$pwm[[idx]]
    attr(m, "nsites") <- x$motifs$nsites[idx]
    attr(m, "label") <- x$motifs$motif[idx]
    return(m)
  }
  if (is.list(x) && !is.data.frame(x)) {
    idx <- if (is.null(motif)) 1L else if (is.numeric(motif)) as.integer(motif) else match(motif, names(x))
    return(x[[idx]])
  }
  as.matrix(x)
}

# ---- logo plot --------------------------------------------------------------

#' Sequence logo for a motif
#'
#' Draws a WebLogo-style stacked-letter logo as plain `ggplot2` polygons. The
#' letters are vector glyphs built by the package itself, so colour schemes,
#' per-letter colours, outlines, transparency and themes are all under your
#' control, and the result composes with [patchwork::wrap_plots()] or
#' [gg_tree_motifs()] like any other plot.
#'
#' @param x A [meme_result], a position probability matrix, or a list of them.
#' @param motif Motif label or index.
#' @param method `"bits"` or `"probability"`, see [logo_heights()].
#' @param scheme Letter colour scheme, see [logo_colours()].
#' @param letter_colours Optional named vector of per-letter colours,
#'   overriding `scheme` for the letters named. `letter_colors` is a synonym.
#' @param letter_colors Synonym for `letter_colours`.
#' @param stroke Glyph stroke thickness, as a fraction of the letter box.
#' @param outline,outline_width Outline colour and width for each glyph.
#'   `outline = NA` (default) draws no outline.
#' @param alpha Glyph opacity.
#' @param gap Horizontal gap between letter stacks, in position units.
#' @param small_sample_correction Schneider small-sample correction.
#' @param reverse_complement Plot the reverse complement (nucleotide motifs
#'   only).
#' @param title,subtitle Plot title and subtitle. `title = NULL` uses the motif
#'   label; `title = NA` draws none.
#' @param y_lab Y axis label; `NULL` picks "Information content (bits)" or
#'   "Probability".
#' @param show_consensus Print the consensus string under the x axis.
#' @param legend Show a letter legend.
#' @param y_limit Upper y limit; defaults to the theoretical maximum
#'   (`log2(alphabet size)` for bits, 1 for probabilities). Used by
#'   [gg_motif_logos()] to put several panels on one scale.
#' @param base_size Base font size for the theme.
#'
#' @return A `ggplot` object.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_logo(res, "MEME-2")
#' gg_motif_logo(res, "MEME-2", scheme = "taylor", outline = "white")
#' gg_motif_logo(res, 1, method = "probability") +
#'   ggplot2::labs(title = "Motif 1, probability stack")
#' @seealso [gg_motif_logos()], [logo_heights()], [logo_colours()]
#' @export
gg_motif_logo <- function(x, motif = NULL, method = c("bits", "probability"),
                          scheme = "auto", letter_colours = NULL,
                          letter_colors = NULL, stroke = 0.19,
                          outline = NA, outline_width = 0.2, alpha = 1,
                          gap = 0.06, small_sample_correction = FALSE,
                          reverse_complement = FALSE,
                          title = NULL, subtitle = NULL, y_lab = NULL,
                          show_consensus = FALSE, legend = FALSE,
                          y_limit = NULL, base_size = 11) {
  method <- match.arg(method)
  letter_colours <- letter_colours %||% letter_colors
  pwm <- .get_pwm(x, motif)
  lab <- attr(pwm, "label") %||% (if (is.character(motif)) motif else "motif")
  if (isTRUE(reverse_complement)) pwm <- .rev_comp_pwm(pwm)

  h <- logo_heights(pwm, method = method,
                    small_sample_correction = small_sample_correction)
  if (!nrow(h)) .abort("No letters above `min_height`; nothing to draw.")
  stack_max <- attr(h, "stack_max")

  polys <- dplyr::bind_rows(lapply(seq_len(nrow(h)), function(i) {
    d <- .place_letter(h$letter[i], h$position[i] - 0.5 + gap / 2,
                       h$position[i] + 0.5 - gap / 2, h$ymin[i], h$ymax[i],
                       stroke = stroke)
    d$letter <- h$letter[i]
    d$grp <- paste(i, d$piece, sep = "_")
    d
  }))

  letters_used <- sort(unique(h$letter))
  cols <- logo_colours(scheme, letters = letters_used)
  if (!is.null(letter_colours)) {
    keep <- intersect(toupper(names(letter_colours)), names(cols))
    cols[keep] <- unname(letter_colours[keep])
  }

  np <- max(h$position)
  if (is.null(y_lab)) {
    y_lab <- if (method == "bits") "Information content (bits)" else "Probability"
  }
  if (is.null(title)) title <- lab
  if (is.na(title)) title <- NULL

  ## both positional scales are built once: adding a second scale_*_continuous
  ## later would work but emits a "scale already present" message
  x_labels <- if (isTRUE(show_consensus)) {
    strsplit(.consensus_from_pwm(as.matrix(pwm)), "")[[1]]
  } else {
    ggplot2::waiver()
  }

  p <- ggplot2::ggplot(polys) +
    ggplot2::geom_polygon(
      ggplot2::aes(x = .data$x, y = .data$y, group = .data$grp,
                   fill = .data$letter),
      colour = outline, linewidth = outline_width, alpha = alpha) +
    ggplot2::scale_fill_manual(values = cols, name = NULL, breaks = letters_used) +
    ggplot2::scale_x_continuous(
      breaks = seq_len(np), labels = x_labels,
      expand = ggplot2::expansion(add = 0.12)) +
    ggplot2::scale_y_continuous(
      limits = c(0, y_limit %||% (if (method == "bits") stack_max else 1)),
      expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(x = "Position", y = y_lab, title = title, subtitle = subtitle) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = if (legend) "right" else "none",
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      plot.title.position = "plot",
      axis.text.x = if (isTRUE(show_consensus)) {
        ggplot2::element_text(family = "mono", size = ggplot2::rel(0.85))
      } else {
        ggplot2::element_text(size = ggplot2::rel(0.75))
      }
    )
  p
}

#' @keywords internal
#' @noRd
.rev_comp_pwm <- function(pwm) {
  m <- as.matrix(pwm)
  cn <- toupper(colnames(m))
  comp <- c(A = "T", T = "A", G = "C", C = "G", U = "A", N = "N")
  if (!all(cn %in% names(comp))) {
    .abort("`reverse_complement = TRUE` only applies to nucleotide motifs.")
  }
  m <- m[nrow(m):1, , drop = FALSE]
  colnames(m) <- unname(comp[cn])
  m[, order(colnames(m)), drop = FALSE]
}

#' Several sequence logos in one figure
#'
#' @param x A [meme_result].
#' @param motifs Motif labels to draw; defaults to all.
#' @param ncol,nrow Grid layout, passed to [patchwork::wrap_plots()].
#' @param shared_y Give every panel the same y axis limit, so panels are
#'   directly comparable.
#' @param ... Passed to [gg_motif_logo()].
#'
#' @return A `patchwork` object.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_logos(res, ncol = 2)
#' @export
gg_motif_logos <- function(x, motifs = NULL, ncol = 1, nrow = NULL,
                           shared_y = TRUE, ...) {
  .need("patchwork")
  if (!inherits(x, "meme_result")) .abort("`x` must be a meme_result.")
  motifs <- motifs %||% x$motifs$motif
  dots <- list(...)
  lim <- if (isTRUE(shared_y)) {
    ## only the arguments logo_heights() understands are forwarded
    hargs <- dots[intersect(names(dots), names(formals(logo_heights))[-(1:2)])]
    max(vapply(motifs, function(m) {
      h <- do.call(logo_heights, c(list(x, m), hargs))
      max(h$ymax, 0)
    }, numeric(1))) * 1.02
  } else {
    NULL
  }
  ## in a stack the panel title has to start at the panel edge, otherwise it
  ## runs into the rotated y axis title of its own panel
  ps <- lapply(motifs, function(m) {
    gg_motif_logo(x, m, y_limit = lim, ...) +
      ggplot2::theme(plot.title.position = "panel")
  })
  patchwork::wrap_plots(ps, ncol = ncol, nrow = nrow)
}
