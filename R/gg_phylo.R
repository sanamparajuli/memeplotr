#' Read a phylogeny
#'
#' Thin wrapper around [ape::read.tree()] that also accepts an existing `phylo`
#' object, a Newick string, or a file path, so that every memeplotr function
#' can take `tree =` in whichever form you have.
#'
#' @param x An `ape::phylo` object, a path to a Newick file, or a Newick string
#'   (recognised by the trailing `;`).
#' @return An `ape::phylo` object.
#' @examples
#' read_tree(system.file("extdata", "example_tree.nwk", package = "memeplotr"))
#' read_tree("((a:1,b:1):1,c:2);")
#' @export
read_tree <- function(x) {
  if (inherits(x, "phylo")) return(x)
  if (inherits(x, "multiPhylo")) return(x[[1]])
  if (!is.character(x) || length(x) != 1L) .abort("`tree` must be a phylo object, a file path, or a Newick string.")
  phy <- if (grepl(";\\s*$", x)) ape::read.tree(text = x) else ape::read.tree(x)
  if (is.null(phy)) .abort("Could not parse the tree.")
  if (inherits(phy, "multiPhylo")) phy <- phy[[1]]
  phy
}

#' @keywords internal
#' @noRd
.tree_layout <- function(phy, ladderize = TRUE, right = TRUE) {
  if (isTRUE(ladderize)) phy <- ape::ladderize(phy, right = right)
  n <- length(phy$tip.label)
  N <- n + phy$Nnode
  has_bl <- !is.null(phy$edge.length)
  xx <- if (has_bl) ape::node.depth.edgelength(phy) else {
    d <- ape::node.depth(phy, method = 1)
    max(d) - d
  }
  parent <- phy$edge[, 1]
  child <- phy$edge[, 2]
  kids <- split(child, parent)
  root <- setdiff(parent, child)[1]

  ## preorder descent gives the conventional top-to-bottom tip order
  ord <- integer(0)
  stack <- root
  while (length(stack)) {
    nd <- stack[1]; stack <- stack[-1]
    ch <- kids[[as.character(nd)]]
    if (is.null(ch)) ord <- c(ord, nd) else stack <- c(ch, stack)
  }
  yy <- rep(NA_real_, N)
  yy[ord] <- n - seq_along(ord) + 1          # first tip at the top
  ## internal nodes sit midway between their extreme children; resolve by
  ## sweeping until every node has a value (equivalent to a postorder pass)
  internals <- setdiff(seq_len(N), seq_len(n))
  for (k in seq_len(N)) {
    todo <- internals[is.na(yy[internals])]
    if (!length(todo)) break
    for (nd in todo) {
      ch <- kids[[as.character(nd)]]
      if (!is.null(ch) && !anyNA(yy[ch])) yy[nd] <- mean(range(yy[ch]))
    }
  }

  nodes <- tibble::tibble(
    node = seq_len(N), x = xx, y = yy,
    is_tip = seq_len(N) <= n,
    label = c(phy$tip.label, phy$node.label %||% rep(NA_character_, phy$Nnode))
  )
  h_seg <- tibble::tibble(x = xx[parent], xend = xx[child],
                          y = yy[child], yend = yy[child])
  v_seg <- dplyr::bind_rows(lapply(unique(parent), function(p) {
    ch <- kids[[as.character(p)]]
    tibble::tibble(x = xx[p], xend = xx[p], y = min(yy[ch]), yend = max(yy[ch]))
  }))
  list(phy = phy, nodes = nodes, h = h_seg, v = v_seg, n_tip = n,
       max_x = max(xx, na.rm = TRUE), has_bl = has_bl,
       tip_order = phy$tip.label[ord])
}

#' @keywords internal
#' @noRd
.tip_order <- function(phy, ladderize = TRUE, right = TRUE) {
  .tree_layout(phy, ladderize = ladderize, right = right)$tip_order
}

#' Draw a phylogeny with ggplot2
#'
#' A rectangular phylogram or cladogram drawn entirely with `ggplot2` segments,
#' so no Bioconductor dependency is required and every element remains
#' modifiable. Tips are laid out one unit apart with the first tip at the top,
#' matching the `y` scale used by [gg_motif_map()] - which is what makes
#' [gg_tree_motifs()] align exactly.
#'
#' @param tree A `phylo` object, Newick file path, or Newick string.
#' @param ladderize,right Ladderize the tree before layout, and in which
#'   direction.
#' @param cladogram Ignore branch lengths and align all tips at the right.
#' @param tip_labels Draw tip labels.
#' @param tip_label_size,tip_label_colour,tip_label_face Tip label appearance.
#' @param tip_label_offset Horizontal offset of tip labels, as a fraction of
#'   tree depth. Only used when `tip_label_position = "panel"`.
#' @param tip_label_position Where tip labels live. `"panel"` draws them inside
#'   the plotting area next to each tip, which is compact but can overflow the
#'   panel when names are long. `"axis"` draws them as right-hand axis text, so
#'   grid reserves their full width outside the panel and they can never run
#'   into an adjacent panel - use this for long sequence names and for
#'   [gg_tree_motifs()].
#' @param align_tips Extend dotted guide lines from each tip to the label
#'   column.
#' @param node_labels Draw internal node labels (e.g. bootstrap support) when
#'   present.
#' @param node_label_size,node_label_colour Node label appearance.
#' @param branch_colour,branch_linewidth Branch appearance. `branch_colour` may
#'   also be a column name in `meta` to colour terminal branches by a trait.
#' @param meta Optional data frame with a `sequence` (or `label`) column used
#'   to colour tips/branches.
#' @param tip_points Draw a point at each tip.
#' @param tip_point_size Size of tip points.
#' @param scale_bar Draw a branch-length scale bar.
#' @param expand_y Vertical expansion, in tip units; keep in sync with
#'   [gg_motif_map()] when combining plots.
#' @param xlim_mult Extra horizontal room for labels, as a multiple of depth.
#'
#' @param title,subtitle Title and subtitle for the panel.
#' @return A `ggplot` object with the attribute `"tip_order"`.
#' @examples
#' nwk <- system.file("extdata", "example_tree.nwk", package = "memeplotr")
#' gg_phylo(nwk)
#' gg_phylo(nwk, cladogram = TRUE, tip_points = TRUE, branch_colour = "steelblue")
#' @export
gg_phylo <- function(tree, ladderize = TRUE, right = TRUE, cladogram = FALSE,
                     tip_labels = TRUE, tip_label_size = 3,
                     tip_label_colour = "grey10", tip_label_face = "plain",
                     tip_label_offset = 0.02,
                     tip_label_position = c("panel", "axis"),
                     align_tips = FALSE,
                     node_labels = FALSE, node_label_size = 2.3,
                     node_label_colour = "grey45",
                     branch_colour = "grey25", branch_linewidth = 0.5,
                     meta = NULL, tip_points = FALSE, tip_point_size = 1.8,
                     scale_bar = FALSE, title = NULL, subtitle = NULL,
                     expand_y = 0.6, xlim_mult = 0.35) {
  tip_label_position <- match.arg(tip_label_position)
  on_axis <- isTRUE(tip_labels) && tip_label_position == "axis"
  phy <- read_tree(tree)
  if (isTRUE(cladogram)) phy$edge.length <- NULL
  lay <- .tree_layout(phy, ladderize = ladderize, right = right)
  tips <- lay$nodes[lay$nodes$is_tip, , drop = FALSE]
  depth <- max(lay$max_x, 1e-9)
  off <- tip_label_offset * depth

  p <- ggplot2::ggplot()
  if (isTRUE(align_tips)) {
    p <- p + ggplot2::geom_segment(
      data = tips,
      ggplot2::aes(x = .data$x, xend = depth, y = .data$y, yend = .data$y),
      linetype = "dotted", colour = "grey80", linewidth = 0.3)
  }
  p <- p +
    ggplot2::geom_segment(data = lay$v,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      colour = branch_colour, linewidth = branch_linewidth, lineend = "round") +
    ggplot2::geom_segment(data = lay$h,
      ggplot2::aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$yend),
      colour = branch_colour, linewidth = branch_linewidth, lineend = "round")

  if (isTRUE(tip_points)) {
    td <- tips
    if (!is.null(meta)) {
      key <- intersect(c("sequence", "label", "tip"), names(meta))[1]
      if (!is.na(key)) {
        td <- dplyr::left_join(td, dplyr::rename(meta, label = !!key), by = "label")
      }
    }
    p <- p + ggplot2::geom_point(data = td,
      ggplot2::aes(x = .data$x, y = .data$y), size = tip_point_size,
      colour = branch_colour)
  }
  if (isTRUE(tip_labels) && !on_axis) {
    xend <- if (isTRUE(align_tips)) depth else NULL
    td <- tips
    td$.lx <- if (is.null(xend)) td$x + off else xend + off
    p <- p + ggplot2::geom_text(data = td,
      ggplot2::aes(x = .data$.lx, y = .data$y, label = .data$label),
      hjust = 0, size = tip_label_size, colour = tip_label_colour,
      fontface = tip_label_face)
  }
  if (isTRUE(node_labels)) {
    nd <- lay$nodes[!lay$nodes$is_tip & !is.na(lay$nodes$label), , drop = FALSE]
    if (nrow(nd)) {
      p <- p + ggplot2::geom_text(data = nd,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = node_label_size, colour = node_label_colour,
        hjust = 1.15, vjust = -0.4)
    }
  }
  if (isTRUE(scale_bar) && lay$has_bl) {
    bl <- signif(depth / 5, 1)
    p <- p + ggplot2::annotate("segment", x = 0, xend = bl,
                               y = 0.15, yend = 0.15, linewidth = 0.5) +
      ggplot2::annotate("text", x = bl / 2, y = 0.15, label = bl,
                        vjust = -0.6, size = 2.6)
  }

  ## y limits are set on the coord (not the scale) so that they match
  ## gg_motif_map() exactly and patchwork can align the two panels
  ## with labels on the axis, grid reserves their width outside the panel, so
  ## no amount of name length can run into a neighbouring panel; in-panel
  ## labels instead need `xlim_mult` room and are drawn with clip = "off"
  ysc <- if (on_axis) {
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0),
                                position = "right",
                                breaks = tips$y, labels = tips$label)
  } else {
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(add = 0))
  }
  p <- p + ysc +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(
        mult = c(0.02, if (isTRUE(tip_labels) && !on_axis) xlim_mult else 0.02))) +
    ggplot2::coord_cartesian(ylim = c(1 - expand_y, lay$n_tip + expand_y),
                             clip = "off") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    theme_tree_motif() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.title.position = "plot")
  if (on_axis) {
    p <- p + ggplot2::theme(
      axis.text.y.right = ggplot2::element_text(
        size = tip_label_size * .pt_per_mm, colour = tip_label_colour,
        face = tip_label_face, hjust = 0,
        margin = ggplot2::margin(l = 2.2)))
  }
  attr(p, "tip_order") <- lay$tip_order
  attr(p, "n_tip") <- lay$n_tip
  p
}
