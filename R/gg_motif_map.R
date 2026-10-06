#' Motif composition (block) diagram
#'
#' Draws one horizontal track per sequence with each motif occurrence as a
#' shaped, coloured block - the MEME "motif locations" figure, but as an
#' ordinary `ggplot2` object that you can keep modifying.
#'
#' @param x A [meme_result] from [read_meme()], a `motif_table` from
#'   [as_motif_table()] / [read_fimo()], or any data frame of hits (passed to
#'   [as_motif_table()] via `...`).
#' @param source For a [meme_result]: `"sites"` (all scanned hits, default) or
#'   `"contributing"` (only the sites MEME used to build the motifs).
#' @param order How to order sequence tracks, top to bottom: `"input"`,
#'   `"tree"` (requires `tree`), `"count"` (most motifs first), `"length"`,
#'   `"name"`, or `"custom"` (uses `seq_order`).
#' @param tree Optional phylogeny: a path to a Newick file or an `ape::phylo`
#'   object. Used only to define tip order when `order = "tree"`; see
#'   [gg_tree_motifs()] to draw the tree alongside.
#' @param seq_order Character vector of sequence names giving an explicit
#'   top-to-bottom order. Implies `order = "custom"`.
#' @param motifs Optional character vector of motif labels to keep (and the
#'   order they appear in the legend).
#' @param sequences Optional character vector of sequences to keep.
#' @param pvalue_max Optional p-value cut-off applied to the hits.
#' @param meta Optional data frame with a `sequence` column plus any extra
#'   columns you want available for `facet_by` or downstream layers.
#' @param facet_by Optional column name (in `meta` or in the hit table) to
#'   facet the tracks by, e.g. a clade or species column.
#' @param shape Block silhouette, one of [motif_shapes()].
#' @param strand_shape If `TRUE`, blocks point in the direction given by the
#'   `strand` column. Ignored for protein motifs, which have no strand.
#' @param height Block height, in track units (1 unit = one sequence).
#' @param head_frac,round_frac Shape tuning, see [geom_motif()].
#' @param palette Palette name, see [motif_palette()].
#' @param colors Optional colours overriding the palette. Either a named vector
#'   (names matching motif labels) or an unnamed vector used in motif order.
#'   `colours` is accepted as a synonym.
#' @param colours Synonym for `colors`.
#' @param alpha Block fill opacity.
#' @param border_colour,border_width Block outline colour and width. Use
#'   `border_colour = NA` for no outline.
#' @param backbone Draw the full-length sequence line behind the blocks.
#' @param backbone_colour,backbone_linewidth Backbone appearance.
#' @param label Label blocks with the motif name (`TRUE`), with nothing
#'   (`FALSE`, default), or with any column name in the hit table.
#' @param label_size,label_colour Block label appearance.
#' @param length_label Print each sequence's length at the right-hand end.
#' @param legend_title Legend title.
#' @param x_lab,title,subtitle Axis label and titles. `x_lab = NULL` picks
#'   "Position (aa)" or "Position (nt)" from the alphabet.
#' @param expand_y Extra space above and below the outermost tracks, in track
#'   units. Keep the default when combining with a tree.
#' @param ... Passed to [as_motif_table()] when `x` is a plain data frame.
#'
#' @param drop_untreed With `order = "tree"`, drop sequences that are absent
#'   from the tree instead of appending them below the tips. [gg_tree_motifs()]
#'   always drops them: an extra row here that the tree panel does not have
#'   would shift every row out of register with its branch and label.
#' @param keep_empty With `order = "tree"`, keep a track for tree tips that
#'   have no motif hits. Needed when the plot is aligned against a tree panel,
#'   and set for you by [gg_tree_motifs()].
#' @return A `ggplot` object. The attribute `"seq_levels"` holds the
#'   top-to-bottom sequence order actually used.
#'
#' @examples
#' res <- read_meme(system.file("extdata", "example_meme.xml", package = "memeplotr"))
#' nwk <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
#'
#' gg_motif_map(res)
#'
#' # everything is a normal ggplot afterwards
#' gg_motif_map(res, shape = "hexagon", palette = "Spectral", order = "count") +
#'   ggplot2::labs(title = "MYB family motif architecture") +
#'   ggplot2::theme(legend.position = "right")
#'
#' # pin individual motifs to your own colours
#' gg_motif_map(res, colors = c("MEME-1" = "#111111", "MEME-3" = "#E64B35"))
#'
#' @seealso [gg_tree_motifs()], [geom_motif()], [scale_fill_motif()]
#' @export
gg_motif_map <- function(x, source = c("sites", "contributing"),
                         order = c("input", "tree", "count", "length", "name", "custom"),
                         tree = NULL, seq_order = NULL, keep_empty = FALSE,
                         drop_untreed = FALSE, motifs = NULL,
                         sequences = NULL, pvalue_max = NULL, meta = NULL,
                         facet_by = NULL,
                         shape = "arrow", strand_shape = FALSE, height = 0.6,
                         head_frac = 0.3, round_frac = 0.35,
                         palette = "okabe_ito", colors = NULL, colours = NULL,
                         alpha = 1, border_colour = "grey25", border_width = 0.25,
                         backbone = TRUE, backbone_colour = "grey75",
                         backbone_linewidth = 0.8,
                         label = FALSE, label_size = 2.4, label_colour = "grey10",
                         length_label = FALSE, legend_title = "Motif",
                         x_lab = NULL, title = NULL, subtitle = NULL,
                         expand_y = 0.6, ...) {
  source <- match.arg(source)
  order <- match.arg(order)
  colors <- colors %||% colours

  alphabet <- if (inherits(x, "meme_result")) x$alphabet$type else NA_character_
  mt <- .coerce_sites(x, source = source, ...)
  all_seqs <- .all_sequences(x, mt)

  mt <- .filter_sites(mt, motifs = motifs, sequences = sequences,
                      pvalue_max = pvalue_max)
  if (!is.null(sequences)) all_seqs <- all_seqs[all_seqs$sequence %in% sequences, , drop = FALSE]
  if (!nrow(mt)) .abort("No motif hits left after filtering.")

  lv <- .seq_levels(all_seqs, mt, order = order, tree = tree,
                    seq_order = seq_order, keep_empty = keep_empty,
                    drop_untreed = drop_untreed)
  ## hits for sequences no longer in `lv` would otherwise become NA rows
  mt <- mt[mt$sequence %in% lv, , drop = FALSE]
  if (!nrow(mt)) .abort("No motif hits left for the sequences being plotted.")
  all_seqs <- all_seqs[match(lv, all_seqs$sequence), , drop = FALSE]
  ## with keep_empty = TRUE some levels have no row in all_seqs; they still get
  ## a track, which is what keeps the rows aligned with a tree panel
  all_seqs$sequence <- lv
  n <- length(lv)
  yof <- stats::setNames(n - seq_len(n) + 1, lv)

  mt <- mt[mt$sequence %in% lv, , drop = FALSE]
  mt$y <- unname(yof[mt$sequence])
  all_seqs$y <- unname(yof[all_seqs$sequence])

  mt$motif <- factor(as.character(mt$motif),
                     levels = motifs %||% .motif_levels(x, mt))
  cols <- .resolve_motif_colors(levels(mt$motif), palette, colors)

  if (!is.null(meta)) {
    mt <- dplyr::left_join(mt, meta, by = "sequence")
    all_seqs <- dplyr::left_join(all_seqs, meta, by = "sequence")
  }

  aes_blocks <- ggplot2::aes(xmin = .data$start, xmax = .data$stop,
                             y = .data$y, fill = .data$motif)
  if (isTRUE(strand_shape)) {
    aes_blocks <- utils::modifyList(aes_blocks, ggplot2::aes(strand = .data$strand))
  }

  p <- ggplot2::ggplot()
  bb <- all_seqs[is.finite(all_seqs$seq_length), , drop = FALSE]
  if (isTRUE(backbone) && nrow(bb)) {
    p <- p + ggplot2::geom_segment(
      data = bb,
      ggplot2::aes(x = 1, xend = .data$seq_length, y = .data$y, yend = .data$y),
      colour = backbone_colour, linewidth = backbone_linewidth,
      lineend = "round", inherit.aes = FALSE)
  }
  p <- p + geom_motif(
    data = mt, mapping = aes_blocks, shape = shape, height = height,
    head_frac = head_frac, round_frac = round_frac,
    colour = border_colour, linewidth = border_width, alpha = alpha)

  lab_col <- if (isTRUE(label)) "motif" else if (is.character(label)) label else NULL
  if (!is.null(lab_col)) {
    if (!lab_col %in% names(mt)) .abort(sprintf("`label` column '%s' not in the hit table.", lab_col))
    mt$.lab <- as.character(mt[[lab_col]])
    mt$.xm <- (mt$start + mt$stop) / 2
    p <- p + ggplot2::geom_text(
      data = mt, ggplot2::aes(x = .data$.xm, y = .data$y, label = .data$.lab),
      size = label_size, colour = label_colour, inherit.aes = FALSE)
  }
  if (isTRUE(length_label) && nrow(bb)) {
    p <- p + ggplot2::geom_text(
      data = bb,
      ggplot2::aes(x = .data$seq_length, y = .data$y,
                   label = format(.data$seq_length, big.mark = ",")),
      hjust = -0.25, size = 2.6, colour = "grey35", inherit.aes = FALSE)
  }

  if (is.null(x_lab)) {
    x_lab <- switch(alphabet %||% "", protein = "Position (aa)",
                    dna = "Position (nt)", rna = "Position (nt)", "Position")
  }

  p <- p +
    ggplot2::scale_fill_manual(values = cols, name = legend_title, drop = FALSE) +
    ggplot2::scale_y_continuous(breaks = unname(yof), labels = names(yof),
                                expand = ggplot2::expansion(add = 0)) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.01, if (length_label) 0.08 else 0.02))) +
    ggplot2::coord_cartesian(ylim = c(1 - expand_y, n + expand_y), clip = "off") +
    ggplot2::labs(x = x_lab, y = NULL, title = title, subtitle = subtitle) +
    theme_motif()

  if (!is.null(facet_by)) {
    if (!facet_by %in% names(mt)) .abort(sprintf("`facet_by` column '%s' not found.", facet_by))
    p <- p + ggplot2::facet_grid(rows = ggplot2::vars(.data[[facet_by]]),
                                 scales = "free_y", space = "free_y", switch = "y")
  }
  attr(p, "seq_levels") <- lv
  p
}

# ---- shared internals -------------------------------------------------------

#' @keywords internal
#' @noRd
.coerce_sites <- function(x, source = "sites", ...) {
  if (inherits(x, "motif_table")) return(x)
  if (inherits(x, "meme_result")) return(as_motif_table(x, source = source))
  if (is.data.frame(x)) return(as_motif_table(x, ...))
  .abort("`x` must be a meme_result, a motif_table, or a data frame of motif hits.")
}

#' @keywords internal
#' @noRd
.all_sequences <- function(x, mt) {
  if (inherits(x, "meme_result") && nrow(x$sequences)) {
    return(tibble::tibble(sequence = x$sequences$sequence,
                          seq_length = as.numeric(x$sequences$length)))
  }
  u <- unique(mt$sequence)
  len <- vapply(u, function(s) {
    v <- mt$seq_length[mt$sequence == s]
    if (all(is.na(v))) max(mt$stop[mt$sequence == s], na.rm = TRUE) else max(v, na.rm = TRUE)
  }, numeric(1))
  tibble::tibble(sequence = u, seq_length = unname(len))
}

#' @keywords internal
#' @noRd
.motif_levels <- function(x, mt) {
  if (inherits(x, "meme_result")) {
    lv <- x$motifs$motif
    return(lv[lv %in% unique(as.character(mt$motif))])
  }
  if (is.factor(mt$motif)) levels(droplevels(mt$motif)) else unique(as.character(mt$motif))
}

#' @keywords internal
#' @noRd
.filter_sites <- function(mt, motifs = NULL, sequences = NULL, pvalue_max = NULL) {
  if (!is.null(motifs)) mt <- mt[as.character(mt$motif) %in% motifs, , drop = FALSE]
  if (!is.null(sequences)) mt <- mt[mt$sequence %in% sequences, , drop = FALSE]
  if (!is.null(pvalue_max)) mt <- mt[is.na(mt$pvalue) | mt$pvalue <= pvalue_max, , drop = FALSE]
  mt
}

#' @keywords internal
#' @noRd
.seq_levels <- function(all_seqs, mt, order = "input", tree = NULL,
                        seq_order = NULL, keep_empty = FALSE,
                        drop_untreed = FALSE) {
  if (!is.null(seq_order)) order <- "custom"
  u <- all_seqs$sequence
  lv <- switch(
    order,
    input = u,
    name = sort(u),
    length = u[base::order(all_seqs$seq_length, decreasing = TRUE)],
    count = {
      k <- table(factor(mt$sequence, levels = u))
      u[base::order(as.numeric(k), decreasing = TRUE)]
    },
    custom = {
      miss <- setdiff(seq_order, u)
      if (length(miss)) {
        rlang::warn(paste0("Dropping unknown names in `seq_order`: ",
                           paste(utils::head(miss, 5), collapse = ", ")))
      }
      c(intersect(seq_order, u), setdiff(u, seq_order))
    },
    tree = {
      if (is.null(tree)) .abort("`order = \"tree\"` needs a `tree`.")
      phy <- read_tree(tree)
      tips <- .tip_order(phy)
      miss <- setdiff(u, tips)
      if (length(miss)) {
        ## `drop_untreed` is what gg_tree_motifs() needs: an extra row here but
        ## not in the tree panel gives the two panels different row counts, so
        ## every label drifts off the track it names.
        where <- if (isTRUE(drop_untreed)) "dropped (absent from the tree panel's rows)" else "appended at the bottom"
        rlang::warn(c(sprintf("%d sequence(s) absent from the tree; %s.", length(miss), where),
                      i = paste(utils::head(miss, 5), collapse = ", "),
                      i = "Tip labels must match the MEME sequence names exactly."))
      }
      if (isTRUE(drop_untreed)) miss <- character()
      if (isTRUE(keep_empty)) c(tips, miss) else c(intersect(tips, u), miss)
    }
  )
  unique(lv)
}
