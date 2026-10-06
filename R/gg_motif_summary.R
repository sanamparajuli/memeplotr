#' @keywords internal
#' @noRd
.prep_tracks <- function(x, source = "sites", order = "input", tree = NULL,
                         seq_order = NULL, keep_empty = FALSE, motifs = NULL,
                         sequences = NULL, pvalue_max = NULL, ...) {
  mt <- .coerce_sites(x, source = source, ...)
  all_seqs <- .all_sequences(x, mt)
  mt <- .filter_sites(mt, motifs = motifs, sequences = sequences,
                      pvalue_max = pvalue_max)
  if (!is.null(sequences)) {
    all_seqs <- all_seqs[all_seqs$sequence %in% sequences, , drop = FALSE]
  }
  if (!nrow(mt)) .abort("No motif hits left after filtering.")
  lv <- .seq_levels(all_seqs, mt, order = order, tree = tree,
                    seq_order = seq_order, keep_empty = keep_empty)
  mlv <- motifs %||% .motif_levels(x, mt)
  mt <- mt[mt$sequence %in% lv, , drop = FALSE]
  mt$motif <- factor(as.character(mt$motif), levels = mlv)
  list(mt = mt, lv = lv, motif_levels = mlv, n = length(lv),
       yof = stats::setNames(length(lv) - seq_along(lv) + 1, lv),
       all_seqs = all_seqs[match(lv, all_seqs$sequence), , drop = FALSE])
}

#' @keywords internal
#' @noRd
.track_y_scale <- function(p, tr, expand_y, labels = TRUE) {
  p +
    ggplot2::scale_y_continuous(
      breaks = unname(tr$yof),
      labels = if (isTRUE(labels)) names(tr$yof) else NULL,
      expand = ggplot2::expansion(add = 0)) +
    ggplot2::coord_cartesian(ylim = c(1 - expand_y, tr$n + expand_y), clip = "off")
}

#' Motif presence or abundance heatmap
#'
#' A sequence-by-motif matrix: useful for spotting clade-specific gains and
#' losses, and for motif sets too dense to read as a composition diagram. With
#' `order = "tree"` the rows line up with [gg_phylo()], so the heatmap can be
#' composed against a tree with [gg_tree_panel()].
#'
#' @param x A [meme_result], `motif_table`, or data frame of hits.
#' @param value `"presence"` for a two-colour present/absent map, `"count"`
#'   for the number of occurrences on a continuous scale.
#' @param order,tree,seq_order,keep_empty Row ordering, as in [gg_motif_map()].
#' @param motifs,sequences,pvalue_max Filters, as in [gg_motif_map()].
#' @param motif_colours Colour present cells by motif rather than by a single
#'   fill. Ignored when `value = "count"`. `motif_colors` is a synonym.
#' @param motif_colors Synonym for `motif_colours`.
#' @param palette Palette used by `motif_colours`, see [motif_palette()].
#' @param present_colour,absent_colour Fills for `value = "presence"`.
#' @param low_colour,high_colour Gradient ends for `value = "count"`.
#' @param label Print the cell value on each tile.
#' @param label_size,label_colour Cell label size and colour.
#' @param border_colour,border_width Tile borders.
#' @param y_labels Draw sequence names on the y axis.
#' @param expand_y Vertical padding in track units; keep it equal to the tree
#'   panel's when composing.
#' @param legend_title,x_lab,title,subtitle Labels.
#' @param ... Passed to [as_motif_table()].
#'
#' @return A `ggplot` object, with the row order in the attribute
#'   `"seq_levels"`.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' nwk <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
#' gg_motif_heatmap(res)
#' gg_motif_heatmap(res, value = "count", label = TRUE)
#' gg_tree_panel(gg_phylo(nwk, tip_labels = FALSE),
#'               gg_motif_heatmap(res, order = "tree", tree = nwk, keep_empty = TRUE))
#' @seealso [gg_motif_map()], [gg_tree_panel()]
#' @export
gg_motif_heatmap <- function(x, value = c("presence", "count"),
                             order = c("input", "tree", "count", "length", "name", "custom"),
                             tree = NULL, seq_order = NULL, keep_empty = FALSE,
                             motifs = NULL, sequences = NULL, pvalue_max = NULL,
                             motif_colours = FALSE, motif_colors = NULL,
                             palette = "okabe_ito",
                             present_colour = "#2C6E91", absent_colour = "grey93",
                             low_colour = "#F2F6F9", high_colour = "#1B4F72",
                             label = FALSE, label_size = 2.7, label_colour = NULL,
                             border_colour = "white", border_width = 0.4,
                             y_labels = TRUE, expand_y = 0.6,
                             legend_title = NULL, x_lab = NULL,
                             title = NULL, subtitle = NULL, ...) {
  value <- match.arg(value)
  order <- match.arg(order)
  motif_colours <- motif_colours %||% motif_colors
  tr <- .prep_tracks(x, order = order, tree = tree, seq_order = seq_order,
                     keep_empty = keep_empty, motifs = motifs,
                     sequences = sequences, pvalue_max = pvalue_max, ...)

  grid <- expand.grid(sequence = tr$lv, motif = tr$motif_levels,
                      stringsAsFactors = FALSE)
  k <- table(factor(tr$mt$sequence, levels = tr$lv),
             factor(as.character(tr$mt$motif), levels = tr$motif_levels))
  grid$n <- as.numeric(k[cbind(grid$sequence, grid$motif)])
  grid$present <- grid$n > 0
  grid$motif <- factor(grid$motif, levels = tr$motif_levels)
  grid$y <- unname(tr$yof[grid$sequence])

  if (value == "count") {
    fill_aes <- ggplot2::aes(fill = .data$n)
    fill_scale <- ggplot2::scale_fill_gradient(
      low = low_colour, high = high_colour,
      name = legend_title %||% "Occurrences")
  } else if (isTRUE(motif_colours)) {
    cols <- .resolve_motif_colors(tr$motif_levels, palette, NULL)
    grid$.f <- factor(ifelse(grid$present, as.character(grid$motif), NA),
                      levels = tr$motif_levels)
    fill_aes <- ggplot2::aes(fill = .data$.f)
    fill_scale <- ggplot2::scale_fill_manual(
      values = cols, na.value = absent_colour, drop = FALSE,
      na.translate = FALSE, name = legend_title %||% "Motif")
  } else {
    fill_aes <- ggplot2::aes(fill = .data$present)
    fill_scale <- ggplot2::scale_fill_manual(
      values = c(`TRUE` = present_colour, `FALSE` = absent_colour),
      labels = c(`TRUE` = "present", `FALSE` = "absent"),
      name = legend_title %||% NULL)
  }

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = .data$motif, y = .data$y)) +
    ggplot2::geom_tile(fill_aes, colour = border_colour, linewidth = border_width) +
    fill_scale
  if (isTRUE(label)) {
    lc <- label_colour %||% "grey20"
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = ifelse(.data$n > 0, format(.data$n), "")),
      size = label_size, colour = lc)
  }
  p <- .track_y_scale(p, tr, expand_y, labels = y_labels) +
    ggplot2::labs(x = x_lab %||% NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_motif(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  attr(p, "seq_levels") <- tr$lv
  p
}

#' Motif counts per sequence
#'
#' A stacked (or dodged) bar chart of how many times each motif occurs in each
#' sequence. Shares the row ordering of [gg_motif_map()], so it can sit beside
#' a tree or a composition diagram.
#'
#' @inheritParams gg_motif_heatmap
#' @param position `"stack"` or `"dodge"`.
#' @param colours,colors Optional named colours for individual motifs.
#' @param alpha,border_colour,border_width Bar styling.
#' @param ... Passed to [as_motif_table()].
#'
#' @return A `ggplot` object, with the row order in the attribute
#'   `"seq_levels"`.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_counts(res)
#' gg_motif_counts(res, position = "dodge", order = "count")
#' @export
gg_motif_counts <- function(x, position = c("stack", "dodge"),
                            order = c("input", "tree", "count", "length", "name", "custom"),
                            tree = NULL, seq_order = NULL, keep_empty = FALSE,
                            motifs = NULL, sequences = NULL, pvalue_max = NULL,
                            palette = "okabe_ito", colours = NULL, colors = NULL,
                            alpha = 1, border_colour = "grey30", border_width = 0.25,
                            y_labels = TRUE, expand_y = 0.6,
                            legend_title = "Motif", x_lab = "Occurrences",
                            title = NULL, subtitle = NULL, ...) {
  position <- match.arg(position)
  order <- match.arg(order)
  tr <- .prep_tracks(x, order = order, tree = tree, seq_order = seq_order,
                     keep_empty = keep_empty, motifs = motifs,
                     sequences = sequences, pvalue_max = pvalue_max, ...)
  d <- as.data.frame(table(sequence = factor(tr$mt$sequence, levels = tr$lv),
                           motif = factor(as.character(tr$mt$motif),
                                          levels = tr$motif_levels)),
                     stringsAsFactors = FALSE, responseName = "n")
  d <- d[d$n > 0, , drop = FALSE]
  d$motif <- factor(d$motif, levels = tr$motif_levels)
  d$y <- unname(tr$yof[d$sequence])
  cols <- .resolve_motif_colors(tr$motif_levels, palette, colours %||% colors)

  ## both axes are continuous, so the orientation has to be declared or
  ## ggplot2 guesses a vertical bar and the widths blow up
  pos <- if (position == "stack") {
    ggplot2::position_stack(reverse = FALSE)
  } else {
    ggplot2::position_dodge2(preserve = "single", padding = 0.1, reverse = TRUE)
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$n, y = .data$y, fill = .data$motif)) +
    ggplot2::geom_col(position = pos, orientation = "y", colour = border_colour,
                      linewidth = border_width, alpha = alpha, width = 0.72) +
    ggplot2::scale_fill_manual(values = cols, name = legend_title, drop = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.04)))
  p <- .track_y_scale(p, tr, expand_y, labels = y_labels) +
    ggplot2::labs(x = x_lab, y = NULL, title = title, subtitle = subtitle) +
    theme_motif(grid = "x")
  attr(p, "seq_levels") <- tr$lv
  p
}

#' Motif co-occurrence
#'
#' How often each pair of motifs is found in the same sequence - the quickest
#' way to see which motifs travel together across a family.
#'
#' @inheritParams gg_motif_heatmap
#' @param measure `"count"` for the number of sequences carrying both motifs,
#'   `"jaccard"` for the intersection over union.
#' @param upper Keep the upper triangle as well as the lower.
#' @param diagonal Keep the diagonal (a motif with itself).
#' @param label Print the value in each cell.
#' @param legend_title Legend title; defaults to the measure name.
#' @param title,subtitle Plot title and subtitle.
#' @param ... Passed to [as_motif_table()].
#'
#' @return A `ggplot` object.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_cooccurrence(res, label = TRUE)
#' gg_motif_cooccurrence(res, measure = "jaccard", upper = TRUE)
#' @export
gg_motif_cooccurrence <- function(x, measure = c("count", "jaccard"),
                                  motifs = NULL, sequences = NULL, pvalue_max = NULL,
                                  upper = FALSE, diagonal = TRUE, label = FALSE,
                                  label_size = 2.8, label_colour = NULL,
                                  low_colour = "#F2F6F9", high_colour = "#1B4F72",
                                  border_colour = "white", border_width = 0.4,
                                  legend_title = NULL, title = NULL, subtitle = NULL,
                                  ...) {
  measure <- match.arg(measure)
  tr <- .prep_tracks(x, motifs = motifs, sequences = sequences,
                     pvalue_max = pvalue_max, ...)
  m <- table(factor(tr$mt$sequence, levels = tr$lv),
             factor(as.character(tr$mt$motif), levels = tr$motif_levels)) > 0
  inter <- t(m) %*% m
  d <- expand.grid(motif_x = tr$motif_levels, motif_y = tr$motif_levels,
                   stringsAsFactors = FALSE)
  d$count <- as.numeric(inter[cbind(d$motif_x, d$motif_y)])
  tot <- colSums(m)
  d$jaccard <- d$count / (tot[d$motif_x] + tot[d$motif_y] - d$count)
  d$jaccard[!is.finite(d$jaccard)] <- 0
  d$value <- d[[measure]]

  ix <- match(d$motif_x, tr$motif_levels); iy <- match(d$motif_y, tr$motif_levels)
  keep <- if (isTRUE(upper)) rep(TRUE, nrow(d)) else ix <= iy
  if (!isTRUE(diagonal)) keep <- keep & ix != iy
  d <- d[keep, , drop = FALSE]
  d$motif_x <- factor(d$motif_x, levels = tr$motif_levels)
  d$motif_y <- factor(d$motif_y, levels = base::rev(tr$motif_levels))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$motif_x, y = .data$motif_y,
                                       fill = .data$value)) +
    ggplot2::geom_tile(colour = border_colour, linewidth = border_width) +
    ggplot2::scale_fill_gradient(
      low = low_colour, high = high_colour,
      name = legend_title %||% if (measure == "count") "Sequences" else "Jaccard")
  if (isTRUE(label)) {
    lab <- if (measure == "count") format(d$value) else sprintf("%.2f", d$value)
    p <- p + ggplot2::geom_text(ggplot2::aes(label = lab),
                                size = label_size,
                                colour = label_colour %||% "grey15")
  }
  p +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    ggplot2::coord_fixed() +
    theme_motif(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Where motifs sit along the sequence
#'
#' Positional distribution of each motif, either in absolute residue
#' coordinates or as a fraction of sequence length, which is the useful view
#' when the sequences differ a lot in length.
#'
#' @inheritParams gg_motif_heatmap
#' @param scale `"absolute"` for residue positions, `"relative"` for
#'   position / sequence length.
#' @param geom `"density"`, `"histogram"`, or `"jitter"` for the raw hits.
#' @param facet Give each motif its own panel.
#' @param bins Number of bins when `geom = "histogram"`.
#' @param alpha,colours,colors,palette Styling, as in [gg_motif_counts()].
#' @param ... Passed to [as_motif_table()].
#'
#' @return A `ggplot` object.
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' gg_motif_positions(res)
#' gg_motif_positions(res, scale = "relative", geom = "histogram", facet = TRUE)
#' @export
gg_motif_positions <- function(x, scale = c("absolute", "relative"),
                               geom = c("density", "histogram", "jitter"),
                               motifs = NULL, sequences = NULL, pvalue_max = NULL,
                               facet = FALSE, bins = 30, alpha = 0.55,
                               palette = "okabe_ito", colours = NULL, colors = NULL,
                               legend_title = "Motif", x_lab = NULL,
                               title = NULL, subtitle = NULL, ...) {
  scale <- match.arg(scale)
  geom <- match.arg(geom)
  tr <- .prep_tracks(x, motifs = motifs, sequences = sequences,
                     pvalue_max = pvalue_max, ...)
  d <- tr$mt
  d$mid <- (d$start + d$stop) / 2
  len <- stats::setNames(tr$all_seqs$seq_length, tr$all_seqs$sequence)
  d$len <- unname(len[d$sequence])
  if (scale == "relative") {
    bad <- !is.finite(d$len) | d$len <= 0
    if (any(bad)) {
      rlang::warn(sprintf("%d hit(s) dropped: sequence length unknown.", sum(bad)))
      d <- d[!bad, , drop = FALSE]
    }
    d$pos <- d$mid / d$len
  } else {
    d$pos <- d$mid
  }
  cols <- .resolve_motif_colors(tr$motif_levels, palette, colours %||% colors)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$pos, fill = .data$motif,
                                       colour = .data$motif))
  p <- p + switch(
    geom,
    ## the default Gaussian rule-of-thumb bandwidth: plug-in selectors fail
    ## outright on motifs with only a couple of hits
    density = ggplot2::geom_density(alpha = alpha, linewidth = 0.4, na.rm = TRUE),
    histogram = ggplot2::geom_histogram(bins = bins, alpha = alpha,
                                        linewidth = 0.2, na.rm = TRUE),
    jitter = ggplot2::geom_jitter(ggplot2::aes(y = .data$motif), height = 0.22,
                                  width = 0, size = 1.6, alpha = max(alpha, 0.7),
                                  na.rm = TRUE)
  )
  p <- p +
    ggplot2::scale_fill_manual(values = cols, name = legend_title, drop = FALSE) +
    ggplot2::scale_colour_manual(values = cols, name = legend_title, drop = FALSE) +
    ggplot2::labs(
      x = x_lab %||% if (scale == "relative") "Relative position" else "Position",
      y = if (geom == "jitter") NULL else if (geom == "density") "Density" else "Hits",
      title = title, subtitle = subtitle) +
    theme_motif(grid = "x")
  if (isTRUE(facet)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$motif), ncol = 2,
                                 scales = "free_y") +
      ggplot2::theme(legend.position = "none")
  }
  p
}
